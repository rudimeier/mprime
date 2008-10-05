; Copyright 2001-2007 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements part of a discrete weighted transform to
; quickly multiply two numbers.
;
; This code handles the last 8 levels of two pass FFTs that use the
; SSE2 instructions.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xmult.mac
INCLUDE xpass2.mac
INCLUDE xnormal.mac

_TEXT SEGMENT

;;
;; Routines to do the normalization after a multiply
;;

;; When doing zero-padded FFTs, the multiplied 7 words around the halfway point
;; must be subtracted from the bottom of the FFT.  This must be done before
;; normalization multiplies the FFT data by k.  This macro does that.

xsub_7_words MACRO
	LOCAL	nozpad, zlp
	cmp	THIS_BLOCK, 8		;; Have we subtracted all 7 words?
	jge	short nozpad		;; Yes, skip this code
	mov	eax, cache_line_multiplier ;; Load loop counter
	mov	edx, THIS_BLOCK
	mov	rcx, rsi		;; Copy source ptr (we preserve rsi)
	mov	rdi, zpad_addr		;; Addr of first zpad element
zlp:	movsd	xmm0, Q [rcx]		;; Load FFT word
	movsd	xmm1, Q [rdi][rdx*8]	;; Load ZPAD data
	mulsd	xmm1, XMM_NORM012_FF	;; Scale by FFTLEN/2
	subsd	xmm0, xmm1
	addsd	xmm7, xmm1		;; Adjust sumout
	movsd	Q [rcx], xmm0		;; Store FFT word
	lea	rcx, [rcx+64]		;; Bump pointers
	inc	rdx
	dec	eax			;; Iterate 2*clm (up to 8) times
	jnz	short zlp		;; Loop if necessary
nozpad:
	ENDM

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

saved_rsi	EQU	PPTR [rsp+first_local]
loopcount1	EQU	DPTR [rsp+first_local+SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+SZPTR+4]
loopcount3	EQU	DPTR [rsp+first_local+SZPTR+8]

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	noadd, setlp, ilp0, ilp1, ilexit, done
	PROCFP	lab
	int_prolog SZPTR+12,0,0
zero	mov	zero_fft, 1		;; Set flag saying zero upper half
	movapd	xmm7, XMM_SUMOUT	;; Load SUMOUT
	movapd	xmm6, XMM_MAXERR	;; Load maximum error
no zero	mov	edx, ADDIN_ROW		;; Is this the time to do our addin?
no zero	cmp	edx, THIS_BLOCK
no zero	jne	short noadd		;; Jump if addin does not occur now
no zero	mov	edi, ADDIN_OFFSET	;; Get address to add value into
no zero	movsd	xmm0, Q [rsi][rdi]	;; Get the value
no zero	addsd	xmm0, ADDIN_VALUE	;; Add in the requested value
no zero	movsd	Q [rsi][rdi], xmm0	;; Save the new value
no zero	subsd	xmm7, ADDIN_VALUE	;; Do not include addin in sumout
noadd:	mov	saved_rsi, rsi		;; Save for xtop_carry_adjust
	mov	rbx, norm_ptr2		;; Load column multipliers ptr
ttp	mov	eax, cache_line_multiplier ;; Load inner loop counter
	lea	rdi, XMM_COL_MULTS	;; Load col mult scratch area
setlp:	xnorm_2d_setup ttp
ttp	lea	rdi, [rdi+512]		;; Next scratch area section
ttp	lea	rbx, [rbx+32]		;; Next column multiplier
ttp	sub	al, 1			;; Each cache line has its own col mult
ttp	jnz	setlp
ttp	mov	norm_ptr2, rbx		;; Save column multipliers ptr

	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2, eax		;; Save loop counter
	mov	loopcount3, 0		;; Clear outermost loop counter
	sub	rax, rax		;; Clear big/lit flags
	sub	rcx, rcx
ttp	mov	al, [rdi+0]		;; Load big vs. little flags
ttp	mov	cl, [rdi+1]		;; Load big vs. little flags
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1, ebx		;; Save loop counter
	lea	rbx, XMM_COL_MULTS	;; Load col mult scratch area
	xprefetcht1 [rdx+128]		;; Prefetch group multiplier
ilp1:	xprefetchw [rsi+64]
	xnorm_2d ttp, zero, echk, const	;; Normalize 8 values
	lea	rsi, [rsi+64]		;; Next cache line
ttp	lea	rbx, [rbx+512]		;; Next column multipliers
ttp	lea	rdi, [rdi+4]		;; Next big/little flags
	sub	loopcount1, 1		;; Test loop counter
	jnz	ilp1			;; Loop til done
	add	rsi, normblkdst		;; Skip gap in blkdst or clmblkdst
	lea	rbp, [rbp+64]		;; Next set of carries
ttp	lea	rdx, [rdx+128]		;; Next set of 8 group multipliers
	sub	loopcount2, 1		;; Test loop counter
	jz	ilexit			;; Jump when loop complete
	add	loopcount3, 80000000h/4 ;; 8 iterations
	jnc	ilp0
	add	rsi, normblkdst8	;; Add 128 every 8 clmblkdsts
	jmp	ilp0			;; Iterate
ilexit:	movapd	XMM_SUMOUT, xmm7	;; Save SUMOUT
	movapd	XMM_MAXERR, xmm6	;; Save maximum error
ttp	mov	norm_ptr1, rdi		;; Save big/little flags array ptr

	; Handle adjusting the carry out of the topmost FFT word

	mov	eax, THIS_BLOCK		;; Check for processing last block
	cmp	eax, LAST_PASS1_BLOCK
	jne	done			;; Jump if not last block
	mov	rsi, saved_rsi		;; Restore FFT data ptr
	xnorm_top_carry			;; Adjust carry if k > 1

done:	int_epilog SZPTR+12,0,0
	ENDPP	lab
	ENDM

loopcount1z	EQU	DPTR [rsp+first_local]
loopcount2z	EQU	DPTR [rsp+first_local+4]
loopcount3z	EQU	DPTR [rsp+first_local+8]

zpnorm	MACRO	lab, ttp, echk, const
	LOCAL	setlp, ilp0, ilp1, ilexit
	PROCFP	lab
	int_prolog 12,0,0
const	mov	const_fft, 1		;; Set flag saying mul-by-const
	movapd	xmm7, XMM_SUMOUT	;; Load SUMOUT
	movapd	xmm6, XMM_MAXERR	;; Load maximum error
	xsub_7_words

	mov	rbx, norm_ptr2		;; Load column multipliers ptr
ttp	mov	eax, cache_line_multiplier ;; Load inner loop counter
	lea	rdi, XMM_COL_MULTS	;; Load col mult scratch area
setlp:	xnorm_2d_setup ttp
ttp	lea	rdi, [rdi+512]		;; Next scratch area section
ttp	lea	rbx, [rbx+32]		;; Next column multiplier
ttp	sub	al, 1			;; Each cache line has its own col mult
ttp	jnz	setlp
ttp	mov	norm_ptr2, rbx		;; Save column multipliers ptr

	mov	rdx, norm_grp_mults	;; Addr of the group multipliers
	mov	rbp, carries		;; Addr of the carries
	mov	rdi, norm_ptr1		;; Load big/little flags array ptr
	mov	eax, addcount1		;; Load loop counter
	mov	loopcount2z, eax	;; Save loop counter
	mov	loopcount3z, 0		;; Clear outermost loop counter
	sub	rax, rax		;; Clear big/lit flags
ttp	mov	al, [rdi+0]		;; Load big vs. little flags
ilp0:	mov	ebx, cache_line_multiplier ;; Load inner loop counter
	mov	loopcount1z, ebx	;; Save loop counter
	lea	rbx, XMM_COL_MULTS	;; Load col mult scratch area
	xprefetcht1 [rdx+128]		;; Prefetch group multiplier
ilp1:	xprefetchw [rsi+64]
	xnorm_2d_zpad ttp, echk, const	;; Normalize 8 values
	lea	rsi, [rsi+64]		;; Next cache line
ttp	lea	rbx, [rbx+512]		;; Next column multipliers
ttp	lea	rdi, [rdi+4]		;; Next big/little flags
	sub	loopcount1z, 1		;; Test loop counter
	jnz	ilp1			;; Loop til done
	add	rsi, normblkdst		;; Skip gap in blkdst or clmblkdst
	lea	rbp, [rbp+64]		;; Next set of carries
ttp	lea	rdx, [rdx+128]		;; Next set of 8 group multipliers
	sub	loopcount2z, 1		;; Test loop counter
	jz	ilexit			;; Jump when loop complete
	add	loopcount3z, 80000000h/4 ;; 8 iterations
	jnc	ilp0
	add	rsi, normblkdst8	;; Add 128 every 8 clmblkdsts
	jmp	ilp0			;; Iterate
ilexit:	movapd	XMM_SUMOUT, xmm7	;; Save SUMOUT
	movapd	XMM_MAXERR, xmm6	;; Save maximum error
ttp	mov	norm_ptr1, rdi		;; Save big/little flags array ptr
	int_epilog 12,0,0
	ENDPP	lab
	ENDM

; The 16 different normalization routines.  One for each combination of
; rational/irrational, zeroing/no zeroing, error check/no error check, and
; mul by const/no mul by const.

PREFETCHING = 1

	inorm	xr2, noexec, noexec, noexec, noexec
	inorm	xr2e, noexec, noexec, exec, noexec
	inorm	xr2c, noexec, noexec, noexec, exec
	inorm	xr2ec, noexec, noexec, exec, exec
	inorm	xr2z, noexec, exec, noexec, noexec
	inorm	xr2ze, noexec, exec, exec, noexec
	inorm	xi2, exec, noexec, noexec, noexec
	inorm	xi2e, exec, noexec, exec, noexec
	inorm	xi2c, exec, noexec, noexec, exec
	inorm	xi2ec, exec, noexec, exec, exec
	inorm	xi2z, exec, exec, noexec, noexec
	inorm	xi2ze, exec, exec, exec, noexec

	zpnorm	xr2zp, noexec, noexec, noexec
	zpnorm	xr2zpe, noexec, exec, noexec
	zpnorm	xr2zpc, noexec, noexec, exec
	zpnorm	xr2zpec, noexec, exec, exec
	zpnorm	xi2zp, exec, noexec, noexec
	zpnorm	xi2zpe, exec, exec, noexec
	zpnorm	xi2zpc, exec, noexec, exec
	zpnorm	xi2zpec, exec, exec, exec

_TEXT	ENDS
END

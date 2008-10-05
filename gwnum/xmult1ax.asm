; Copyright 2001-2007 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
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
INCLUDE memory.mac
INCLUDE xnormal.mac

_TEXT SEGMENT

;;
;; Add two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCF	gwxaddq1
	ad_prolog 0,0,rsi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	eax, addcount1		; Load loop counter
uaddlp:	movapd	xmm0, [rdx]		; Load second number
	addpd	xmm0, [rcx]		; Add in first number
	movapd	xmm1, [rdx+16]		; Load second number
	addpd	xmm1, [rcx+16]		; Add in first number
	movapd	xmm2, [rdx+32]		; Load second number
	addpd	xmm2, [rcx+32]		; Add in first number
	movapd	xmm3, [rdx+48]		; Load second number
	addpd	xmm3, [rcx+48]		; Add in first number
	movapd	[rsi], xmm0		; Save result
	movapd	[rsi+16], xmm1		; Save result
	movapd	[rsi+32], xmm2		; Save result
	movapd	[rsi+48], xmm3		; Save result
	lea	rcx, [rcx+64]		; Next source
	lea	rdx, [rdx+64]		; Next source
	lea	rsi, [rsi+64]		; Next dest
	sub	eax, 1			; Check loop counter
	jnz	short uaddlp		; Loop if necessary
	ad_epilog 0,0,rsi
gwxaddq1 ENDP

;;
;; Add two numbers with carry propogation
;;

saved_dest_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_col_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+2*SZPTR+0]
loopcount2	EQU	DPTR [rsp+first_local+2*SZPTR+4]

PROCF	gwxadd1
	ad_prolog 2*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	movapd	xmm2, XMM_BIGVAL	; Start process with no carry
	movapd	xmm3, xmm2
	mov	eax, normcount1		; Load loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
nadd0:	mov	loopcount1, eax		; Save loop counter
	and	eax, 07FFh		; Grab 11 bits of the counter
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for xnorm_op_1d_mid
	mov	saved_col_ptr, rbp	; remember rbp for xnorm_op_1d_mid
	sub	rax, rax		; Clear big/lit flags
	sub	rbx, rbx
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	naddlp			; Yes, use two-to-phi multipliers
raddlp:	xnorm_op_1d addpd, noexec	; Add and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	raddlp			; Loop til done
	jmp	nadddn			; Loop til done
naddlp:	xnorm_op_1d addpd, exec		; Add and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	naddlp			; Loop til done
nadddn:	mov	rbx, saved_col_ptr	; Restore multipliers pointer
	mov	rax, saved_dest_ptr	; Restore dest pointer
	xnorm_op_1d_mid_cleanup		; Rotate carries and add in carries
	mov	eax, loopcount1		; Restore loop counter
	shr	eax, 11			; Get next loop amount
	jnz	nadd0
	mov	rsi, DESTARG		; Address of result
	mov	rbx, norm_col_mults	; Address of the multipliers
	xnorm_op_1d_cleanup
	ad_epilog 2*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxadd1 ENDP

;;
;; Subtract two numbers without carry propogation.  Caller can use this for
;; consecutive add or subtract operations.  However, the last operation
;; before a multiply must use the routine that will normalize data.
;;

PROCF	gwxsubq1
	ad_prolog 0,0,rsi
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	eax, addcount1		; Load loop counter
usublp:	movapd	xmm0, [rdx]		; Load second number
	subpd	xmm0, [rcx]		; Subtract first number
	movapd	xmm1, [rdx+16]		; Load second number
	subpd	xmm1, [rcx+16]		; Subtract first number
	movapd	xmm2, [rdx+32]		; Load second number
	subpd	xmm2, [rcx+32]		; Subtract first number
	movapd	xmm3, [rdx+48]		; Load second number
	subpd	xmm3, [rcx+48]		; Subtract first number
	movapd	[rsi], xmm0		; Save result
	movapd	[rsi+16], xmm1		; Save result
	movapd	[rsi+32], xmm2		; Save result
	movapd	[rsi+48], xmm3		; Save result
	lea	rcx, [rcx+64]		; Next source
	lea	rdx, [rdx+64]		; Next source
	lea	rsi, [rsi+64]		; Next dest
	sub	eax, 1			; Check loop counter
	jnz	short usublp		; Loop if necessary
	ad_epilog 0,0,rsi
gwxsubq1 ENDP

;;
;; Subtract two numbers with carry propogation
;;

saved_dest_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_col_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+2*SZPTR+0]
loopcount2	EQU	DPTR [rsp+first_local+2*SZPTR+4]

PROCF	gwxsub1
	ad_prolog 2*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	movapd	xmm2, XMM_BIGVAL	; Start process with no carry
	movapd	xmm3, xmm2
	mov	eax, normcount1		; Load loop counter
	mov	rbp, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
nsub0:	mov	loopcount1, eax		; Save loop counter
	and	eax, 07FFh		; Grab 11 bits of the counter
	mov	loopcount2, eax
	mov	saved_dest_ptr, rsi	; remember rsi for xnorm_op_1d_mid
	mov	saved_col_ptr, rbp	; remember rbp for xnorm_op_1d_mid
	sub	rax, rax		; Clear big/lit flag
	sub	rbx, rbx
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	nsublp			; Yes, use two-to-phi multipliers
rsublp:	xnorm_op_1d subpd, noexec	; Subtract and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	rsublp			; Loop til done
	jmp	nsubdn			; Jump if done
nsublp:	xnorm_op_1d subpd, exec		; Subtract and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	nsublp			; Loop til done
nsubdn:	mov	rbx, saved_col_ptr	; Restore multipliers pointer
	mov	rax, saved_dest_ptr	; Restore dest pointer
	xnorm_op_1d_mid_cleanup		; Rotate carries and add in carries
	mov	eax, loopcount1		; Restore loop counter
	shr	eax, 11			; Get next loop amount
	jnz	nsub0
	mov	rsi, DESTARG		; Address of result
	mov	rbx, norm_col_mults	; Address of the multipliers
	xnorm_op_1d_cleanup
	ad_epilog 2*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxsub1 ENDP

;;
;; Add and subtract two numbers without carry propogation.
;;

PROCF	gwxaddsubq1
	ad_prolog 0,0,rbp,rsi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination #1
	mov	rbp, DEST2ARG	  	; Address of destination #2
	mov	eax, addcount1		; Load loop counter
uaddsublp:
	movapd	xmm0, [rcx]		; Load first number
	movapd	xmm1, xmm0		; Dup first number
	addpd	xmm0, [rdx]		; Add in second number
	subpd	xmm1, [rdx]		; Subtract out second number
	movapd	xmm2, [rcx+16]		; Load first number
	movapd	xmm3, xmm2		; Dup first number
	addpd	xmm2, [rdx+16]		; Add in second number
	subpd	xmm3, [rdx+16]		; Subtract out second number
	movapd	xmm4, [rcx+32]		; Load first number
	movapd	xmm5, xmm4		; Dup first number
	addpd	xmm4, [rdx+32]		; Add in second number
	subpd	xmm5, [rdx+32]		; Subtract out second number
	movapd	xmm6, [rcx+48]		; Load first number
	movapd	xmm7, xmm6		; Dup first number
	addpd	xmm6, [rdx+48]		; Add in second number
	subpd	xmm7, [rdx+48]		; Subtract out second number
	movapd	[rsi], xmm0		; Save result
	movapd	[rbp], xmm1		; Save result
	movapd	[rsi+16], xmm2		; Save result
	movapd	[rbp+16], xmm3		; Save result
	movapd	[rsi+32], xmm4		; Save result
	movapd	[rbp+32], xmm5		; Save result
	movapd	[rsi+48], xmm6		; Save result
	movapd	[rbp+48], xmm7		; Save result
	lea	rcx, [rcx+64]		; Next source
	lea	rdx, [rdx+64]		; Next source
	lea	rsi, [rsi+64]		; Next dest
	lea	rbp, [rbp+64]		; Next dest
	sub	eax, 1			; Check loop counter
	jnz	uaddsublp		; Loop if necessary
	ad_epilog 0,0,rbp,rsi,xmm6,xmm7
gwxaddsubq1 ENDP

;;
;; Add and subtract two numbers with carry propogation
;;

saved_dest1_ptr	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_dest2_ptr	EQU	PPTR [rsp+first_local+1*SZPTR]
saved_col_ptr	EQU	PPTR [rsp+first_local+2*SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+3*SZPTR+0]
loopcount2	EQU	DPTR [rsp+first_local+3*SZPTR+4]

PROCF	gwxaddsub1
	ad_prolog 3*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
	mov	rcx, SRCARG		; Address of first number
	mov	rdx, SRC2ARG		; Address of second number
	mov	rsi, DESTARG		; Address of destination
	mov	rbp, DEST2ARG	  	; Address of destination #2
	movapd	xmm2, XMM_BIGVAL	; Start process with no carry
	movapd	xmm3, xmm2
	movapd	xmm6, xmm2
	movapd	xmm7, xmm2
	mov	eax, normcount1		; Load loop counter
	mov	rbx, norm_col_mults	; Address of the multipliers
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
naddsub0:
	mov	loopcount1, eax		; Save loop counter
	and	eax, 07FFh		; Grab 11 bits of the counter
	mov	loopcount2, eax
	mov	saved_col_ptr, rbx	; remember rbx for xnorm_op_1d_mid
	mov	saved_dest2_ptr, rbp	; remember dest #2 for xnorm_op_1d_mid
	mov	saved_dest1_ptr, rsi	; remember dest #1 for xnorm_op_1d_mid
	sub	rax, rax		; Clear big/lit flag
	cmp	RATIONAL_FFT, 0		; Test for irrational FFTs
	je	naddsublp		; Yes, use two-to-phi multipliers
raddsublp:
	xnorm_addsub_1d noexec		; Add/sub and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	raddsublp		; Loop til done
	jmp	naddsubdn		; Jump if done
naddsublp:
	xnorm_addsub_1d exec		; Add/sub and normalize 8 values
	sub	loopcount2, 1		; Decrement loop counter
	jnz	naddsublp		; Loop til done
naddsubdn:
	xchg	rbx, saved_col_ptr	; Save/Restore multipliers pointer
					; Rotate carries and add in carries
	xnorm_addsub_1d_mid_cleanup saved_dest1_ptr, saved_dest2_ptr
	mov	rbx, saved_col_ptr	; Restore multipliers pointer
	mov	eax, loopcount1		; Restore loop counter
	shr	eax, 11			; Get next loop amount
	jnz	naddsub0

	mov	rsi, DESTARG		; Address of result #1
	mov	rbp, DEST2ARG		; Address of result #2
	mov	rbx, norm_col_mults	; Address of the multipliers
	xnorm_addsub_1d_cleanup
	ad_epilog 3*SZPTR+8,0,rbx,rbp,rsi,rdi,xmm6,xmm7
gwxaddsub1 ENDP

;;
;; Copy one number and zero some low order words.
;;

PROCF	gwxcopyzero1
	ad_prolog 0,0,rsi,rdi
	mov	rsi, SRCARG		; Address of first number
	mov	rdi, DESTARG		; Address of destination
	sub	ecx, ecx		; Offset to compare to COPYZERO
	mov	eax, addcount1		; Load loop counter
cz1:	xcopyzero			; Copy/zero 8 values
	lea	rsi, [rsi+64]		; Next source
	lea	rdi, [edi+64]		; Next dest
	lea	rcx, [rcx+64]		; Next compare offset
	sub	eax, 1			; Test loop counter
	jnz	cz1			; Loop if necessary
	ad_epilog 0,0,rsi,rdi
gwxcopyzero1 ENDP

;;
;; Routines to do the normalization after a multiply
;;

; Macro to loop through all the FFT values and apply the proper normalization
; routine.

saved_reg1	EQU	PPTR [rsp+first_local+0*SZPTR]
saved_reg2	EQU	PPTR [rsp+first_local+1*SZPTR]
saved_reg3	EQU	PPTR [rsp+first_local+2*SZPTR]

inorm	MACRO	lab, ttp, zero, echk, const
	LOCAL	ilp0, ilp1
	PROCF	lab
	int_prolog 3*SZPTR,0,0
	mov	rsi, DESTARG		;; Addr of multiplied number
no zero	mov	edi, ADDIN_OFFSET	;; Get address to add value into
no zero	movsd	xmm0, Q [rsi][rdi]	;; Get the value
no zero	addsd	xmm0, ADDIN_VALUE	;; Add in the requested value
no zero	movsd	Q [rsi][rdi], xmm0	;; Save the new value
no zero	subsd	xmm7, ADDIN_VALUE	;; Do not include addin in sumout
	movapd	xmm2, XMM_BIGVAL	;; Start process with no carry
	movapd	xmm3, xmm2
	movlpd	xmm6, MAXERR		;; Current maximum error
	movhpd	xmm6, MAXERR
	mov	rbp, norm_col_mults	;; Addr of the multipliers
	mov	rdi, norm_biglit_array	;; Addr of the big/little flags array
	sub	rax, rax		;; Clear big/lit flags
	sub	rcx, rcx
	mov	ebx, normcount1		;; Load loop counter
ilp0:	mov	saved_reg1, rbx		;; Save loop counter
	and	ebx, 07FFh		;; Grab 11 bits of the counter
	mov	saved_reg2, rdi		;; remember edi for xnorm012_1d_mid
	mov	saved_reg3, rsi		;; remember esi for xnorm012_1d_mid
	mov	rdx, rbp		;; remember ebp for xnorm012_1d_mid
ilp1:	xnorm_1d ttp, zero, echk, const ;; Normalize 8 values
	lea	rsi, [rsi+64]		;; Next cache line
ttp	lea	rbp, [rbp+128]		;; Next set of 8 multipliers
ttp	lea	rdi, [rdi+4]		;; Next big/little flags
	sub	ebx, 1			;; Test loop counter
	jnz	ilp1			;; Loop til done
	mov	rbx, saved_reg3		;; Restore FFT data addr
	xchg	rdi, saved_reg2		;; Restore big/lit pointer
	xnorm012_1d_mid zero		;; Rotate carries and add in carries
	mov	rdi, saved_reg2		;; Restore big/lit pointer
	mov	rbx, saved_reg1		;; Restore loop counter
	shr	ebx, 11			;; Get next loop amount
	jnz	ilp0
zero	jmp	zdn			;; Go to zero upper half end code
no zero	jmp	idn			;; Go to normal end code
&lab	ENDP
	ENDM

zpnorm	MACRO	lab, ttp, echk, const
	LOCAL	ilp0, ilp1
	PROCF	lab
	int_prolog 3*SZPTR,0,0
	mov	rsi, DESTARG		;; Addr of multiplied number
	movapd	xmm2, XMM_BIGVAL	;; Start process with no carry
	subpd	xmm3, xmm3
	movlpd	xmm6, MAXERR		;; Current maximum error
	movhpd	xmm6, MAXERR
	mov	rbp, norm_col_mults	;; Addr of the multipliers
	mov	rdi, norm_biglit_array	;; Addr of the big/little flags array
	sub	rax, rax		;; Clear big/lit flag
	mov	ebx, normcount1		;; Load loop counter
ilp0:	mov	saved_reg1, rbx		;; Save loop counter
	and	ebx, 07FFh		;; Grab 11 bits of the counter
	mov	saved_reg2, rdi		;; remember edi for xnorm012_1d_mid
	mov	saved_reg3, rsi		;; remember esi for xnorm012_1d_mid
	mov	rdx, rbp		;; remember ebp for xnorm012_1d_mid
ilp1:	xnorm_1d_zpad ttp, echk, const	;; Normalize 8 values
	lea	rsi, [rsi+64]		;; Next cache line
ttp	lea	rbp, [rbp+128]		;; Next set of 8 multipliers
ttp	lea	rdi, [rdi+4]		;; Next big/little flags
	sub	ebx, 1			;; Test loop counter
	jnz	ilp1			;; Loop til done
	mov	rbx, saved_reg3		;; Restore FFT data addr
	xchg	rdi, saved_reg2		;; Restore big/lit pointer
	xnorm012_1d_mid_zpad const	;; Rotate carries and add in carries
	mov	rdi, saved_reg2		;; Restore big/lit pointer
	mov	rbx, saved_reg1		;; Restore loop counter
	shr	ebx, 11			;; Get next loop amount
	jnz	ilp0
const	jmp	zpcdn			;; Go to zero padded FFT end code
no const jmp	zpdn			;; Go to zero padded FFT end code
&lab	ENDP
	ENDM

; The 16 different normalization routines.  One for each combination of
; rational/irrational, zeroing/no zeroing, error check/no error check, and
; mul by const/no mul by const.

	inorm	xr1, noexec, noexec, noexec, noexec
	inorm	xr1e, noexec, noexec, exec, noexec
	inorm	xr1c, noexec, noexec, noexec, exec
	inorm	xr1ec, noexec, noexec, exec, exec
	inorm	xr1z, noexec, exec, noexec, noexec
	inorm	xr1ze, noexec, exec, exec, noexec
	inorm	xi1, exec, noexec, noexec, noexec
	inorm	xi1e, exec, noexec, exec, noexec
	inorm	xi1c, exec, noexec, noexec, exec
	inorm	xi1ec, exec, noexec, exec, exec
	inorm	xi1z, exec, exec, noexec, noexec
	inorm	xi1ze, exec, exec, exec, noexec

	zpnorm	xr1zp, noexec, noexec, noexec
	zpnorm	xr1zpe, noexec, exec, noexec
	zpnorm	xr1zpc, noexec, noexec, exec
	zpnorm	xr1zpec, noexec, exec, exec
	zpnorm	xi1zp, exec, noexec, noexec
	zpnorm	xi1zpe, exec, exec, noexec
	zpnorm	xi1zpc, exec, noexec, exec
	zpnorm	xi1zpec, exec, exec, exec

; Common code to finish off the one-pass FFTs normalization.  The
; Windows 64-bit ABI frowns on us jumping from one procedure into another.
; However, my reading of the spec is that as long as the two procedures have
; identical prologs then stack unwinding for exception handling will work OK.

PROCF	__common_xnorm1_end_code

	;; Dummy prolog to match normalization code
	int_prolog 3*SZPTR,0,0

; Finish off the normalization process by adding any carry to first values.
; Handle both the with and without two-to-phi array cases.

zpcdn:	mov	rsi, DESTARG		; Address of squared number
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	mov	rbp, norm_col_mults	; Restart the column multipliers
	sub	rax, rax
	xnorm012_1d_zpad exec		; Add in carries
	jmp	cmnend			; All done, go cleanup

zpdn:	mov	rsi, DESTARG		; Address of squared number
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	mov	rbp, norm_col_mults	; Restart the column multipliers
	sub	rax, rax
	xnorm012_1d_zpad noexec		; Add in carries
	jmp	cmnend			; All done, go cleanup

zdn:	mov	rsi, DESTARG		; Address of squared number
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	mov	rbp, norm_col_mults	; Restart the column multipliers
	sub	rax, rax
	xnorm012_1d exec		; Add in carries
	jmp	cmnend			; All done, go cleanup

idn:	mov	rsi, DESTARG		; Address of squared number
	xnorm_top_carry_1d		; Adjust top carry when k > 1
	mov	rdi, norm_biglit_array	; Address of the big/little flags array
	mov	rbp, norm_col_mults	; Restart the column multipliers
	sub	rax, rax
	xnorm012_1d noexec		; Add in carries

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

cmnend:	movapd	XMM_TMP1, xmm7		; Add together the two partial sumouts
	addsd	xmm7, XMM_TMP1+8
	mulsd	xmm7, ttmp_ff_inv
	movsd	Q [rsi-24], xmm7	; Save sum of FFT outputs
	movapd	XMM_TMP1, xmm6		; Compute new maximum error
	maxsd	xmm6, XMM_TMP1+8
	movsd	MAXERR, xmm6

; Return

	int_epilog 3*SZPTR,0,0
__common_xnorm1_end_code ENDP

_TEXT	ENDS
END

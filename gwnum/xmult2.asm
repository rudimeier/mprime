; Copyright 2001-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; FFT sizes between than 5K and 128K doubles are supported.
; This code does two passes, 8 levels on the second pass.
;
; You will not stand a chance of understanding any of this code without
; thoroughly familiarizing yourself with fast fourier transforms.  This
; code was adapted from an algorithm described in Richard Crandall's article
; on Discrete Weighted Transforms and Large-Integer Arithmetic.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xarch.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE xnormal.mac

PUBLIC	xgw_finish_fft
PUBLIC	xgw_carries
PUBLIC	xgw_finish_mult
PUBLIC	pass1_aux_entry_point_return

_TEXT SEGMENT

;;*************************************************************************
;; Routine for auxillary threads to call to start processing pass 1 blocks
;;*************************************************************************

; pass1_aux_entry_point ()
; Entry point for auxillary threads to do process blocks in pass 1
; Windows 32-bit and Linux 32-bit
;	Parameter asm_data = [esp+4]
; Windows 64-bit
;	Parameter asm_data = rcx
; Linux 64-bit
;	Parameter asm_data = rdi

PROCFL pass1_aux_entry_point
	ad_prolog 0,1,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

;; Jump to common pass 1 code.  We must jump rather than call because the
;; pass 1 code operates with a push_amt of zero.

	mov	rax, THREAD_WORK_ROUTINE ; Pass 1 entry point
	jmp	rax			; Go process data blocks
pass1_aux_entry_point ENDP

;;*************************************************************************
;; Routine for auxillary threads to call to start processing pass 2 blocks
;;*************************************************************************

; pass2_aux_entry_point ()
; Entry point for auxillary threads to process blocks in pass 2
; Windows 32-bit and Linux 32-bit
;	Parameter asm_data = [esp+4]
; Windows 64-bit
;	Parameter asm_data = rcx
; Linux 64-bit
;	Parameter asm_data = rdi

PROCFL	pass2_aux_entry_point
	ad_prolog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

;; Call intermediate routine to mimic main thread's prologs.  We must do this
;; so that common epilog code performs properly.

	call	internal_pass2_aux_entry_point

;; Return to C code

	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15
pass2_aux_entry_point ENDP

; Intermediate routine to call to set up an identical prolog as the
; main thread's processing of pass 2 data blocks.

PROCF	internal_pass2_aux_entry_point
	int_prolog 0,0,1

;; Jump to common code handling pass 2 data blocks.

	mov	rax, THREAD_WORK_ROUTINE ; Pass 2 entry point
	jmp	rax			; Go process data blocks

internal_pass2_aux_entry_point ENDP


;;*****************************************
;; Routine for finishing off a two-pass FFT
;;*****************************************

; Split the accumulated carries into two carries - a high carry and a
; low carry.  Handle both the with and without two-to-phi array cases.
; Add these carries back into the FFT data.

loopcount1	EQU	DPTR [rsp+first_local]
loopcount2	EQU	DPTR [rsp+first_local+4]

PROCF	xgw_carries
	int_prolog 8,0,0
	cmp	ZERO_PADDED_FFT, 0	; Special case the zero padded FFT case
	jne	xgw_carries_zpad
	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Compute addr after carries
	shl	rbx, 6
	add	rbx, rsi
	xnorm012_2d_part1
	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load 3 section counts
	mov	loopcount1, eax		; Save for later
ilp0:	mov	eax, loopcount1		; Get list of counts
	mov	ebx, eax		; Form count for this section
	and	rbx, 07FFh
	jz	spldn			; No rows to do.  We're all done!
	mov	loopcount2, ebx		; Save count of carry rows this section
	shr	eax, 11			; Move counts list along
	mov	loopcount1, eax
	shl	rbx, 6			; Compute addr of the last carries row
	add	rbx, rsi
	xnorm012_2d_part2
	cmp	B_IS_2, 0		; Is b = 2?
	jne	ilp1			; Yes, do simpler roundings
nb2ilp1:mov	rbx, norm_col_mults	; Addr of the column multipliers
	xnorm012_2d noexec		; Split carries for one cache line
	mov	ebx, cache_line_multiplier; Cache lines in each pass1 loop
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	cmp	RATIONAL_FFT, 0		; Don't bump these two pointers
	jne	short nb2iskip		; for rational FFTs
	bump	rdx, 128		; Next group multiplier
	lea	rdi, [rdi+rbx*4]	; Next big/little flags pointer
nb2iskip:sub	loopcount2, 1		; Test loop counter
	jnz	nb2ilp1			; Next carry row in section
	jmp	ilp0			; Next section
ilp1:	mov	rbx, norm_col_mults	; Addr of the column multipliers
	xnorm012_2d exec		; Split carries for one cache line
	mov	ebx, cache_line_multiplier; Cache lines in each pass1 loop
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	cmp	RATIONAL_FFT, 0		; Don't bump these two pointers
	jne	short iskip		; for rational FFTs
	bump	rdx, 128		; Next group multiplier
	lea	rdi, [rdi+rbx*4]	; Next big/little flags pointer
iskip:	sub	loopcount2, 1		; Test loop counter
	jnz	ilp1			; Next carry row in section
	jmp	ilp0			; Next section
spldn:	mov	zero_fft, 0		; Clear zero-high-words-fft flag
	jmp	cdn			; Jump to common exit code

xgw_carries_zpad:
	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Compute addr after carries
	shl	rbx, 6
	add	rbx, rsi
	mov	rsi, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rbp, norm_col_mults	; Addr of the group multipliers
	xnorm012_2d_zpad_part1
	mov	rsi, carries		; Addr of the carries
	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load 3 section counts
	mov	loopcount1, eax		; Save for later
zlp0:	mov	eax, loopcount1		; Get list of counts
	mov	ebx, eax		; Form count for this section
	and	rbx, 07FFh
	jz	zpldn			; No rows to do.  We're all done!
	mov	loopcount2, ebx		; Save count of carry rows this section
	shr	eax, 11			; Move counts list along
	mov	loopcount1, eax
	shl	rbx, 6			; Compute addr of the last carries row
	add	rbx, rsi
	xnorm012_2d_zpad_part2
	movapd	XMM_TMP1, xmm6		; xnorm012_2d_zpad needs to use this register
zlp1:	mov	rbx, norm_col_mults	; Addr of the column multipliers
	xnorm012_2d_zpad		; Split carries for one cache line
c2d:	mov	ebx, cache_line_multiplier; Cache lines in each pass1 loop
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	cmp	RATIONAL_FFT, 0		; Don't bump these two pointers
	jne	short zskip		; for rational FFTs
	bump	rdx, 128		; Next group multiplier
	lea	rdi, [rdi+rbx*4]	; Next big/little flags pointer
zskip:	sub	loopcount2, 1		; Test loop counter
	jnz	zlp1			; Next carry row in section
	movapd	xmm6, XMM_TMP1		; Restore saved register
	jmp	zlp0			; Next section
zpldn:	mov	const_fft, 0		; Clear mul-by-const-fft flag

cdn:	int_epilog 8,0,0
xgw_carries ENDP


; Common code to finish off the two-pass FFTs.  The Windows 64-bit ABI
; frowns on us jumping from one procedure into another.
; However, my reading of the spec is that as long as the two procedures have
; identical prologs then stack unwinding for exception handling will work OK.
; Of course, this code won't be linked into a 64-bit Windows executable,
; but we include the dummy prolog to be consistent.

PROCF	__common_2pass_xfft_exit_code

	;; Create a dummy prolog
	ad_prolog 0,1,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

; Common code to finish up ffts

xgw_finish_fft:
	mov	rsi, DESTARG		; Restore source pointer
	mov	DWORD PTR [rsi-28], 3	; Set has-been-FFTed flags
	xfft_1_ret

; Common code to finish up multiplies

; Finish the multiply

xgw_finish_mult:

; Set FFT-started flag

	mov	rsi, DESTARG		; Addr of FFT data
	mov	eax, POSTFFT		; Set FFT started flag
	mov	DWORD PTR [rsi-28], eax

; Normalize SUMOUT value by multiplying by 1 / (fftlen/2).

	movsd	xmm7, XMM_SUMOUT	; Add together the two partial sumouts
	addsd	xmm7, XMM_SUMOUT+8
	mulsd	xmm7, ttmp_ff_inv
	movsd	Q [rsi-24], xmm7	; Save sum of FFT outputs
	movsd	xmm6, XMM_MAXERR	; Compute new maximum error
	maxsd	xmm6, XMM_MAXERR+8
	movsd	MAXERR, xmm6
	;; Fall through to return code

; Common routine for pass1 auxillary threads to exit

pass1_aux_entry_point_return:
	ad_epilog 0,1,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

__common_2pass_xfft_exit_code ENDP

_TEXT	ENDS
END

; Copyright 2001-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; These routine implement some common cleanup code for r4dwpn (r4delay with partial normalization) FFTs
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

PUBLIC	xgw_carries_wpn
PUBLIC	xgw_finish_mult_wpn

_TEXT SEGMENT

;;*****************************************
;; Routine for finishing off a r4dwpn FFT
;;*****************************************

; Split the accumulated carries into two carries - a high carry and a
; low carry.  Handle both the with and without two-to-phi array cases.
; Add these carries back into the FFT data.

biglit_incr	EQU	PPTR [rsp+first_local]
grp_incr	EQU	PPTR [rsp+first_local+SZPTR]
loopcount1	EQU	DPTR [rsp+first_local+2*SZPTR]
loopcount2	EQU	DPTR [rsp+first_local+2*SZPTR+4]

PROCF	xgw_carries_wpn
	int_prolog 2*SZPTR+8,0,0
	mov	eax, cache_line_multiplier; Cache lines in each pass1 loop
	shl	rax, 1			; Compute biglit increment
	mov	edx, 4*XMM_GMD		; Compute grp increment
	cmp	RATIONAL_FFT, 0		; Don't bump these two pointers
	je	short iskip		; for rational FFTs
	sub	rax, rax		; Zero biglit_incr
	sub	rdx, rdx		; Zero grp_incr
iskip:	mov	biglit_incr, rax	; Save computed increments
	mov	grp_incr, rdx
	mov	rsi, carries		; Addr of the carries
	mov	ebx, addcount1		; Load block count
	shl	rbx, 6
	add	rbx, rsi
	cmp	ZERO_PADDED_FFT, 0	; Special case the zero padded FFT case
	jne	xgw_carries_zpad
	xnorm012_wpn_part1
	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
ilp0:	mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
ilp1:	cmp	B_IS_2, 0		; Is b = 2?
	jne	b2			; Yes, do simpler roundings
	xnorm012_wpn noexec		; Split carries for one cache line
	jmp	nb2			; Rejoin common code
b2:	xnorm012_wpn exec		; Split carries for one cache line
nb2:	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	ilp1
	add	rdx, grp_incr		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	ilp0			; Next carry row
idn:	mov	zero_fft, 0		; Clear zero-high-words-fft flag
	jmp	cdn			; Jump to common exit code

xgw_carries_zpad:
	mov	rsi, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	xnorm012_wpn_zpad_part1
	mov	rsi, carries		; Addr of the carries
	mov	rbp, DESTARG		; Addr of the FFT data
	mov	rdi, norm_biglit_array	; Addr of the big/little flags array
	mov	rdx, norm_grp_mults	; Addr of the group multipliers
	mov	eax, count3		; Load count of grp multipliers
	mov	loopcount2, eax		; Save count
zlp0:	mov	eax, count2		; Load wpn count
	mov	loopcount1, eax		; Save count
zlp1:	xnorm012_wpn_zpad		; Split carries for one cache line
	bump	rsi, 64			; Next carries pointer
	add	rbp, pass1blkdst	; Next FFT data pointer
	add	rdi, biglit_incr	; Next big/little flags pointer
	sub	loopcount1, 1		; Test loop counter
	jnz	zlp1
	add	rdx, grp_incr		; Next group multiplier
	sub	loopcount2, 1		; Test loop counter
	jnz	zlp0			; Next carry row
	mov	const_fft, 0		; Clear mul-by-const-fft flag

cdn:	int_epilog 2*SZPTR+8,0,0
xgw_carries_wpn ENDP


; Common code to finish off r4dwpn FFTs.  The Windows 64-bit ABI
; frowns on us jumping from one procedure into another.
; However, my reading of the spec is that as long as the two procedures have
; identical prologs then stack unwinding for exception handling will work OK.
; Of course, this code won't be linked into a 64-bit Windows executable,
; but we include the dummy prolog to be consistent.

PROCF	__common_wpn_xfft_exit_code

	;; Create a dummy prolog
	ad_prolog 0,1,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

; Common code to finish up multiplies

; Finish the multiply

xgw_finish_mult_wpn:

; Set FFT-started flag

	mov	rsi, DESTARG		; Addr of FFT data
	mov	eax, POSTFFT		; Set FFT started flag
	mov	DWORD PTR [rsi-28], eax

; Calculate SUMOUT and MAXERR values

	movsd	xmm7, XMM_SUMOUT	; Add together the two partial sumouts
	addsd	xmm7, XMM_SUMOUT+8
	movsd	Q [rsi-24], xmm7	; Save sum of FFT outputs
	movsd	xmm6, XMM_MAXERR	; Compute new maximum error
	maxsd	xmm6, XMM_MAXERR+8
	movsd	MAXERR, xmm6

; Return

	ad_epilog 0,1,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

__common_wpn_xfft_exit_code ENDP

_TEXT	ENDS
END

; Copyright 2001-2007 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements part of a discrete weighted transform to
; quickly multiply two numbers.
;
; This code handles the last 8-13 levels of two pass FFTs that use the
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
INCLUDE	xlucas.mac
INCLUDE xpass2.mac

_TEXT SEGMENT

PREFETCHING = 1

;;*****************************************************
;; Macros to do the code common to all pass 2 routines
;;*****************************************************

;; Entry point for real FFTs: do one real and count1 complex blocks, also
;; entry point for all-complex FFTs: do count1 complex blocks

xpass2_entry MACRO complex_start, get_nxt
	int_prolog 0,0,1
	MOVOFFSET rax, get_nxt		;; Auxillary thread entry point
	mov	THREAD_WORK_ROUTINE, rax ;; save addr for C code to call
	c_call	PASS2_WAKE_UP_THREADS	;; C callback routine
	cmp	ALL_COMPLEX_FFT, 1	;; Test if there is a real-data block
	je	complex_start		;; Jump to process all-complex blocks
	ENDM

;; Call C code to get next block to process.  Loop to process the block
;; or return when work completes.

xpass2_loop MACRO complex_loop

;; GET_NEXT_BLOCK returns TRUE if we are done.  Otherwise, it calculates
;; the next blocks to process and prefetch.

	c_call	PASS2_GET_NEXT_BLOCK	;; C callback routine
	and	rax, rax		;; Test return code
	jz	complex_loop		;; If false, process another block

;; All done

	int_epilog 0,0,1
	ENDM

;;*************************************************************************
;; Routine for auxillary threads to call to start processing pass 2 blocks
;;*************************************************************************

IFNDEF AMD

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

ENDIF

;;*****************************************************
;; Routine to do the last 8 levels in a two-pass FFT
;;*****************************************************

PROCFP	xpass2_8_levels
	xpass2_entry p8lp, p8get	; Set up block counts and pointers
	xpass2_8_levels_real		; Do the real data block
	jmp	p8get			; Go process complex blocks
p8lp:	xpass2_8_levels_complex
p8get:	xpass2_loop p8lp
ENDPP	xpass2_8_levels

;;*****************************************************
;; Routine to do the last 10 levels in a two-pass FFT
;;*****************************************************

PROCFP	xpass2_10_levels
	xpass2_entry p10lp, p10get	; Set up block counts and pointers
	xpass2_10_levels_real		; Do the real data block
	jmp	p10get			; Go process complex blocks
p10lp:	xpass2_10_levels_complex
p10get:	xpass2_loop p10lp
ENDPP	xpass2_10_levels


;;*****************************************************
;; Routine to do the last 11 levels in a two-pass FFT
;;*****************************************************

PROCFP	xpass2_11_levels
	xpass2_entry p11lp, p11get	; Set up block counts and pointers
	xpass2_11_levels_real		; Do the real data block
	jmp	p11get			; Go process complex blocks
p11lp:	xpass2_11_levels_complex
p11get:	xpass2_loop p11lp
ENDPP	xpass2_11_levels

;;*****************************************************
;; Routine to do the last 12 levels in a two-pass FFT
;;*****************************************************

PROCFP	xpass2_12_levels
	xpass2_entry p12lp, p12get	; Set up block counts and pointers
	xpass2_12_levels_real		; Do the real data block
	jmp	p12get			; Go process complex blocks
p12lp:	xpass2_12_levels_complex
p12get:	xpass2_loop p12lp
ENDPP	xpass2_12_levels


;;*****************************************************
;; Routine to do the last 13 levels in a two-pass FFT
;;*****************************************************

PROCFP	xpass2_13_levels
	xpass2_entry p13lp, p13get	; Set up block counts and pointers
	xpass2_13_levels_real		; Do the real data block
	jmp	p13get			; Go process complex blocks
p13lp:	xpass2_13_levels_complex
p13get:	xpass2_loop p13lp
ENDPP	xpass2_13_levels


PREFETCHING = 0

;;*************************************************************************
;; Routine to do the last 11 levels in a two-pass FFT without prefetching
;;*************************************************************************

PROCFP	xpass2_11_levels_np
	xpass2_entry np11lp, np11get	; Set up block counts and pointers
	xpass2_11_levels_real		; Do the real data block
	jmp	np11get			; Go process complex blocks
np11lp:	xpass2_11_levels_complex
np11get:xpass2_loop np11lp
ENDPP	xpass2_11_levels_np


;;*************************************************************************
;; Routine to do the last 12 levels in a two-pass FFT without prefetching
;;*************************************************************************

PROCFP	xpass2_12_levels_np
	xpass2_entry np12lp, np12get	; Set up block counts and pointers
	xpass2_12_levels_real		; Do the real data block
	jmp	np12get			; Go process complex blocks
np12lp:	xpass2_12_levels_complex
np12get:xpass2_loop np12lp
ENDPP	xpass2_12_levels_np

;;*************************************************************************
;; Routine to do the last 13 levels in a two-pass FFT without prefetching
;;*************************************************************************

PROCFP	xpass2_13_levels_np
	xpass2_entry np13lp, np13get	; Set up block counts and pointers
	xpass2_13_levels_real		; Do the real data block
	jmp	np13get			; Go process complex blocks
np13lp:	xpass2_13_levels_complex
np13get:xpass2_loop np13lp
ENDPP	xpass2_13_levels_np


_TEXT	ENDS
END

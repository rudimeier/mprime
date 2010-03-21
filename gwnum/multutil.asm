; Copyright 2001-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine used to be part of mult.asm.  However, linux x86-64 objcopy had
; difficulty with "mov rax, OFFSET variable" when variable was declared in the
; same file.
;

	TITLE   multutil

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE xmult.mac

VERSION_NUMBER = 2514		;; Version 25.14

;; Make prctabs public.  Linux x86-64 objcopy handles loading offsets of
;; variables only if they are public variables in one file and defined
;; extrn wherever they are used.

IFNDEF X86_64
PUBLIC	prctab1, prctab2
ENDIF
PUBLIC	xprctab1, xprctab2, xprctab2a

  ;; Routines located in other asm files

IFNDEF X86_64
EXTRN	gwadd1:PROC
EXTRN	gwaddq1:PROC
EXTRN	gwsub1:PROC
EXTRN	gwsubq1:PROC
EXTRN	gwaddsub1:PROC
EXTRN	gwaddsubq1:PROC
EXTRN	gwcopyzero1:PROC
EXTRN	gwmuls1:PROC
EXTRN	gwadd2:PROC
EXTRN	gwaddq2:PROC
EXTRN	gwsub2:PROC
EXTRN	gwsubq2:PROC
EXTRN	gwaddsub2:PROC
EXTRN	gwaddsubq2:PROC
EXTRN	gwcopyzero2:PROC
EXTRN	gwmuls2:PROC
ENDIF

EXTRN	gwxadd1:PROC
EXTRN	gwxaddq1:PROC
EXTRN	gwxsub1:PROC
EXTRN	gwxsubq1:PROC
EXTRN	gwxaddsub1:PROC
EXTRN	gwxaddsubq1:PROC
EXTRN	gwxcopyzero1:PROC
EXTRN	gwxadds1:PROC
EXTRN	gwxmuls1:PROC
EXTRN	gwxadd2:PROC
EXTRN	gwxaddq2:PROC
EXTRN	gwxsub2:PROC
EXTRN	gwxsubq2:PROC
EXTRN	gwxaddsub2:PROC
EXTRN	gwxaddsubq2:PROC
EXTRN	gwxcopyzero2:PROC
EXTRN	gwxadds2:PROC
EXTRN	gwxmuls2:PROC

;; List of possible suffixes for the normalization routines

exnorm	MACRO num, prefix, suffix
	exnorm1	&prefix, &num&z&suffix
	exnorm1	&prefix, &num&ze&suffix
	ENDM
exnorm1	MACRO prefix, suffix
	EXTRN	&prefix&r&suffix:PROC
	EXTRN	&prefix&i&suffix:PROC
	ENDM

IFNDEF X86_64
exnorm 1
exnorm 2
ENDIF
exnorm 1, x
exnorm 2, x
exnorm 2, x, AMD

; Jmptable addresses for gwinfo1 to return

IFNDEF X86_64
EXTRN	jmptable:DWORD
EXTRN	jmptablep:DWORD
ENDIF
EXTRN	xjmptable:DWORD
EXTRN	xjmptablep:DWORD

;
; More global variables needed in FFT routines
;

IFNDEF X86_64
_GWDATA SEGMENT PAGE PUBLIC 'DATA'
ELSE
_GWDATA SEGMENT PAGE
ENDIF

IFNDEF X86_64
prctab1	DD	OFFSET gwadd1, OFFSET gwaddq1, OFFSET gwsub1, OFFSET gwsubq1
	DD	OFFSET gwaddsub1, OFFSET gwaddsubq1, OFFSET gwcopyzero1
	DD	OFFSET gwaddq1, OFFSET gwsubq1, OFFSET gwaddsubq1
	DD	0, OFFSET gwmuls1 ;; Buggy, won't propogate carry over multiple words. pad on zeropads
	DD	OFFSET r1z, OFFSET r1ze, OFFSET i1z, OFFSET i1ze
prctab2	DD	OFFSET gwadd2, OFFSET gwaddq2, OFFSET gwsub2, OFFSET gwsubq2
	DD	OFFSET gwaddsub2, OFFSET gwaddsubq2, OFFSET gwcopyzero2
	DD	OFFSET gwaddq2, OFFSET gwsubq2, OFFSET gwaddsubq2
	DD	0, OFFSET gwmuls2 ;; Buggy, won't propogate carry over multiple words. pad on zeropads
	DD	OFFSET r2z, OFFSET r2ze, OFFSET i2z, OFFSET i2ze
ENDIF

xprctab1 DP	OFFSET gwxadd1, OFFSET gwxaddq1, OFFSET gwxsub1
 	DP	OFFSET gwxsubq1, OFFSET gwxaddsub1, OFFSET gwxaddsubq1
	DP	OFFSET gwxcopyzero1
	DP	OFFSET gwxaddq1, OFFSET gwxsubq1, OFFSET gwxaddsubq1
	DP	OFFSET gwxadds1, OFFSET gwxmuls1
	DP	OFFSET xr1z, OFFSET xr1ze, OFFSET xi1z, OFFSET xi1ze
xprctab2 DP	OFFSET gwxadd2, OFFSET gwxaddq2, OFFSET gwxsub2
	DP	OFFSET gwxsubq2, OFFSET gwxaddsub2, OFFSET gwxaddsubq2
	DP	OFFSET gwxcopyzero2
	DP	OFFSET gwxaddq2, OFFSET gwxsubq2, OFFSET gwxaddsubq2
	DP	OFFSET gwxadds2, OFFSET gwxmuls2
	DP	OFFSET xr2z, OFFSET xr2ze, OFFSET xi2z, OFFSET xi2ze
xprctab2a DP	OFFSET gwxadd2, OFFSET gwxaddq2, OFFSET gwxsub2
	DP	OFFSET gwxsubq2, OFFSET gwxaddsub2, OFFSET gwxaddsubq2
	DP	OFFSET gwxcopyzero2
	DP	OFFSET gwxaddq2, OFFSET gwxsubq2, OFFSET gwxaddsubq2
	DP	OFFSET gwxadds2, OFFSET gwxmuls2
	DP	OFFSET xr2zAMD, OFFSET xr2zeAMD, OFFSET xi2zAMD, OFFSET xi2zeAMD

	;; Align so that other GWDATA areas are also aligned on a cache line
	align 128
_GWDATA ENDS


;; FFT setup routines

_TEXT SEGMENT

; gwinfo1 (resptr)
;	Return address of jmp tables for C code to examine
; Windows 32-bit (_gwinfo1)
; Linux 32-bit (gwinfo1)
;	Parameter resptr = [esp+4]
; Windows 64-bit (gwinfo1) - leaf routine, no unwind info necessary
;	Parameter resptr = rcx
; Linux 64-bit (gwinfo1)
;	Parameter resptr = rdi

PROCL	gwinfo1
	IFNDEF X86_64
	mov	ecx, [esp+4]		; Address of data struct to return info
	ENDIF
	IFDEF LINUX64
	mov	rcx, rdi		; Address of data struct to return info
	ENDIF
	mov	rax, OFFSET xjmptable	; P4 mersenne mod FFTs
	mov	[rcx+0*SZPTR], rax
	mov	rax, OFFSET xjmptablep	; P4 2^N+1 mod FFTs
	mov	[rcx+1*SZPTR], rax
	IFNDEF X86_64
	mov	rax, OFFSET jmptable	; x86 mersenne mod FFTs
	mov	[rcx+2*SZPTR], rax
	mov	rax, OFFSET jmptablep	; x86 2^N+1 mod FFTs
	mov	[rcx+3*SZPTR], rax
	ENDIF
	mov	eax, VERSION_NUMBER
	mov	[rcx+4*SZPTR], eax
	ret
gwinfo1 ENDP

_TEXT	ENDS
END

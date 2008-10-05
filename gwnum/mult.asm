; Copyright 1995-2008 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements the setup, common routines, and global variables
; for the various discrete-weighted transforms
;

	TITLE   gwdata

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE xmult.mac

;; Make jmptables public.  Linux x86-64 objcopy handles loading offsets of
;; variables only if they are public variables in one file and defined
;; extrn wherever they are used.

IFNDEF X86_64
PUBLIC	jmptable, jmptablep
ENDIF
PUBLIC	xjmptable, xjmptablep

;; Proctables tables declared in another file to workaround linux x86-64
;; objcopy bugs

IFNDEF X86_64
EXTRN	prctab1:DWORD
EXTRN	prctab2:DWORD
ENDIF
EXTRN	xprctab1:DWORD
EXTRN	xprctab2:DWORD
EXTRN	xprctab2a:DWORD

;; List of possible suffixes for the FFTs

SFX_486		EQU	1	;; 486/Pentium version - no suffix

SFX_PPRO	EQU	2	;; Pentium Pro / P2 optimized
SFX_PPRO_CLM2	EQU	4

SFX_P3		EQU	8	;; P3 and Athlon optimized
SFX_P3_CLM2	EQU	10h

SFX_P4		EQU	20h	;; P4 optimized (SSE2, no suffix)
SFX_P4_CLM0	EQU	40h
SFX_P4_CLM1	EQU	80h
SFX_P4_CLM2	EQU	100h
SFX_P4_CLM4	EQU	200h
SFX_P4_CLM8	EQU	400h
SFX_P4_CLM10	EQU	SFX_P4_CLM1 + SFX_P4_CLM0
SFX_P4_CLM20	EQU	SFX_P4_CLM2 + SFX_P4_CLM0
SFX_P4_CLM21	EQU	SFX_P4_CLM2 + SFX_P4_CLM1
SFX_P4_CLM210	EQU	SFX_P4_CLM2 + SFX_P4_CLM1 + SFX_P4_CLM0
SFX_P4_CLM41	EQU	SFX_P4_CLM4 + SFX_P4_CLM1
SFX_P4_CLM42	EQU	SFX_P4_CLM4 + SFX_P4_CLM2
SFX_P4_CLM420	EQU	SFX_P4_CLM4 + SFX_P4_CLM2 + SFX_P4_CLM0
SFX_P4_CLM421	EQU	SFX_P4_CLM4 + SFX_P4_CLM2 + SFX_P4_CLM1
SFX_P4_CLM4210	EQU	SFX_P4_CLM4 + SFX_P4_CLM2 + SFX_P4_CLM1 + SFX_P4_CLM0
SFX_P4_CLM84	EQU	SFX_P4_CLM8 + SFX_P4_CLM4
SFX_P4_CLM8421	EQU	SFX_P4_CLM8 + SFX_P4_CLM4 + SFX_P4_CLM2 + SFX_P4_CLM1

SFX_AMD64	EQU	800h	;; AMD64 optimized (SSE2, AMD suffix)
SFX_AMD64_CLM0	EQU	1000h
SFX_AMD64_CLM1	EQU	2000h
SFX_AMD64_CLM2	EQU	4000h
SFX_AMD64_CLM4	EQU	8000h
SFX_AMD64_CLM8	EQU	10000h
SFX_AMD64_CLM10 EQU	SFX_AMD64_CLM1 + SFX_AMD64_CLM0
SFX_AMD64_CLM20 EQU	SFX_AMD64_CLM2 + SFX_AMD64_CLM0
SFX_AMD64_CLM21 EQU	SFX_AMD64_CLM2 + SFX_AMD64_CLM1
SFX_AMD64_CLM210 EQU	SFX_AMD64_CLM2 + SFX_AMD64_CLM1 + SFX_AMD64_CLM0
SFX_AMD64_CLM41	EQU	SFX_AMD64_CLM4 + SFX_AMD64_CLM1
SFX_AMD64_CLM42	EQU	SFX_AMD64_CLM4 + SFX_AMD64_CLM2
SFX_AMD64_CLM421 EQU	SFX_AMD64_CLM4 + SFX_AMD64_CLM2 + SFX_AMD64_CLM1
SFX_AMD64_CLM4210 EQU	SFX_AMD64_CLM4 + SFX_AMD64_CLM2 + SFX_AMD64_CLM1 + SFX_AMD64_CLM0
SFX_AMD64_CLM84 EQU	SFX_AMD64_CLM8 + SFX_AMD64_CLM4
SFX_AMD64_CLM8421 EQU	SFX_AMD64_CLM8 + SFX_AMD64_CLM4 + SFX_AMD64_CLM2 + SFX_AMD64_CLM1

exfft	MACRO fft_length, x, levels
	exfft1	fft_length, _1, x, 0, 1, 2, 4, 8, levels
	exfft1	fft_length, _2, x, 0, 1, 2, 4, 8, levels
	exfft1	fft_length, _3, x, 0, 1, 2, 4, 8, levels
	exfft1	fft_length, _4, x, 0, 1, 2, 4, 8, levels
	ENDM
exfft1	MACRO fft_length, suffix, x, c0, c1, c2, c4, c8, levels
	IFNDEF X86_64
	IF x AND SFX_486
	EXTRN	fft&fft_length&suffix:PROC
	ENDIF
	IF x AND SFX_PPRO
	EXTRN	fft&fft_length&suffix&PPRO:PROC
	ENDIF
	IF x AND SFX_PPRO_CLM2
	EXTRN	fft&fft_length&c2&suffix&PPRO:PROC
	ENDIF
	IF x AND SFX_P3
	EXTRN	fft&fft_length&suffix&P3:PROC
	ENDIF
	IF x AND SFX_P3_CLM2
	EXTRN	fft&fft_length&c2&suffix&P3:PROC
	ENDIF
	ENDIF
	IF x AND SFX_P4
	EXTRN	xfft&fft_length&suffix:PROC
	ENDIF
	IF x AND SFX_P4_CLM0
	EXTRN	xfft&fft_length&c0&levels&suffix:PROC
	ENDIF
	IF x AND SFX_P4_CLM1
	EXTRN	xfft&fft_length&c1&levels&suffix:PROC
	ENDIF
	IF x AND SFX_P4_CLM2
	EXTRN	xfft&fft_length&c2&levels&suffix:PROC
	ENDIF
	IF x AND SFX_P4_CLM4
	EXTRN	xfft&fft_length&c4&levels&suffix:PROC
	ENDIF
	IF x AND SFX_P4_CLM8
	EXTRN	xfft&fft_length&c8&levels&suffix:PROC
	ENDIF
	IF x AND SFX_AMD64
	EXTRN	xfft&fft_length&suffix&AMD:PROC
	ENDIF
	IF x AND SFX_AMD64_CLM0
	EXTRN	xfft&fft_length&c0&levels&suffix&AMD:PROC
	ENDIF
	IF x AND SFX_AMD64_CLM1
	EXTRN	xfft&fft_length&c1&levels&suffix&AMD:PROC
	ENDIF
	IF x AND SFX_AMD64_CLM2
	EXTRN	xfft&fft_length&c2&levels&suffix&AMD:PROC
	ENDIF
	IF x AND SFX_AMD64_CLM4
	EXTRN	xfft&fft_length&c4&levels&suffix&AMD:PROC
	ENDIF
	IF x AND SFX_AMD64_CLM8
	EXTRN	xfft&fft_length&c8&levels&suffix&AMD:PROC
	ENDIF
	ENDM

	exfft	32, SFX_PPRO + SFX_P4
	exfft	40, SFX_PPRO
	exfft	48, SFX_PPRO + SFX_P4
	exfft	56, SFX_PPRO
	exfft	64, SFX_PPRO + SFX_P4
	exfft	80, SFX_PPRO + SFX_P4
	exfft	96, SFX_PPRO + SFX_P4
	exfft	112, SFX_PPRO + SFX_P4
	exfft	128, SFX_PPRO + SFX_P4
	exfft	160, SFX_PPRO + SFX_P4
	exfft	192, SFX_PPRO + SFX_P4
	exfft	224, SFX_PPRO + SFX_P4
	exfft	256, SFX_PPRO + SFX_P4
	exfft	320, SFX_PPRO + SFX_P4
	exfft	384, SFX_PPRO + SFX_P4
	exfft	448, SFX_PPRO + SFX_P4
	exfft	512, SFX_PPRO + SFX_P4
	exfft	640, SFX_PPRO + SFX_P4
	exfft	768, SFX_PPRO + SFX_P4
	exfft	896, SFX_PPRO + SFX_P4
	exfft	1024, SFX_PPRO + SFX_P4
	exfft	1280, SFX_PPRO + SFX_P4
	exfft	1536, SFX_PPRO + SFX_P4
	exfft	1792, SFX_PPRO + SFX_P4
	exfft	2048, SFX_PPRO + SFX_P4
	exfft	2560, SFX_PPRO + SFX_P4
	exfft	3072, SFX_PPRO + SFX_P4
	exfft	3584, SFX_PPRO + SFX_P4
	exfft	4096, SFX_PPRO + SFX_P4
	exfft	5120, SFX_PPRO + SFX_P4
	exfft	6144, SFX_PPRO + SFX_P4
	exfft	7168, SFX_PPRO + SFX_P4
	exfft	8192, SFX_PPRO + SFX_P4
	exfft	10K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	12K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	14K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	16K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	20K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	24K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	28K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	32K, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	40K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	48K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	56K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	64K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	80K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	96K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	112K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	128K, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	160K, SFX_PPRO + SFX_P3
allfft	exfft	160K, SFX_P4 + SFX_AMD64
allfft	exfft	160K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	160K, SFX_P4_CLM2 + SFX_AMD64_CLM2, 10
	exfft	192K, SFX_PPRO + SFX_P3
allfft	exfft	192K, SFX_P4 + SFX_AMD64
allfft	exfft	192K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	192K, SFX_P4_CLM41 + SFX_AMD64_CLM2, 10
	exfft	224K, SFX_PPRO + SFX_P3
allfft	exfft	224K, SFX_P4 + SFX_AMD64
allfft	exfft	224K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	224K, SFX_P4_CLM42 + SFX_AMD64_CLM2, 10
	exfft	256K, SFX_PPRO + SFX_P3
allfft	exfft	256K, SFX_P4 + SFX_AMD64
allfft	exfft	256K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	256K, SFX_P4_CLM42 + SFX_AMD64_CLM2, 10
	exfft	320K, SFX_PPRO + SFX_P3
allfft	exfft	320K, SFX_P4 + SFX_AMD64
allfft	exfft	320K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
allfft	exfft	320K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 11
	exfft	320K, SFX_P4_CLM41 + SFX_AMD64_CLM1, 10
	exfft	384K, SFX_PPRO + SFX_P3
allfft	exfft	384K, SFX_P4 + SFX_AMD64
allfft	exfft	384K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
allfft	exfft	384K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 11
	exfft	384K, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	448K, SFX_PPRO + SFX_P3
allfft	exfft	448K, SFX_P4 + SFX_AMD64
allfft	exfft	448K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
allfft	exfft	448K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 11
	exfft	448K, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	512K, SFX_PPRO + SFX_P3
allfft	exfft	512K, SFX_P4 + SFX_AMD64
allfft	exfft	512K, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 10
allfft	exfft	512K, SFX_P4_CLM421 + SFX_AMD64_CLM4, 11
	exfft	512K, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	640K, SFX_PPRO + SFX_P3
allfft	exfft	640K, SFX_P4_CLM1, 8
allfft	exfft	640K, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 10
allfft	exfft	640K, SFX_P4_CLM8421 + SFX_AMD64_CLM84, 11
allfft	exfft	640K, SFX_P4_CLM8421 + SFX_AMD64_CLM84, 11
allfft	exfft	640K9, SFX_P4_CLM1, 11
allfft	exfft	640K, SFX_P4_CLM42 + SFX_AMD64_CLM8421, 12
	exfft	640K, SFX_P4_CLM41, 11
	exfft	640K, SFX_P4_CLM2 + SFX_AMD64_CLM2, 12
	exfft	768K, SFX_PPRO + SFX_P3
allfft	exfft	768K, SFX_P4_CLM1, 8
allfft	exfft	768K, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 10
allfft	exfft	768K, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 11
allfft	exfft	768K9, SFX_P4_CLM10, 11
allfft	exfft	768K, SFX_P4_CLM42 + SFX_AMD64_CLM421, 12
	exfft	768K, SFX_P4_CLM0, 10
	exfft	768K, SFX_P4_CLM42, 11
	exfft	768K, SFX_P4_CLM2 + SFX_AMD64_CLM2, 12
	exfft	896K, SFX_PPRO + SFX_P3
allfft	exfft	896K, SFX_P4_CLM1, 8
allfft	exfft	896K, SFX_P4_CLM210 + SFX_AMD64_CLM2, 10
allfft	exfft	896K, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 11
allfft	exfft	896K9, SFX_P4_CLM10, 11
allfft	exfft	896K, SFX_P4_CLM42 + SFX_AMD64_CLM421, 12
	exfft	896K, SFX_P4_CLM420, 11
	exfft	896K, SFX_P4_CLM2 + SFX_AMD64_CLM21, 12
	exfft	1024K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	1024K, SFX_P4_CLM1, 8
allfft	exfft	1024K, SFX_P4_CLM210 + SFX_AMD64_CLM2, 10
allfft	exfft	1024K, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 11
allfft	exfft	1024K9, SFX_P4_CLM10, 11
allfft	exfft	1024K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
	exfft	1024K, SFX_P4_CLM420, 11
	exfft	1024K, SFX_P4_CLM1 + SFX_AMD64_CLM21, 12
	exfft	1280K, SFX_PPRO + SFX_P3
allfft	exfft	1280K, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	1280K, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 11
allfft	exfft	1280K9, SFX_P4_CLM10, 11
allfft	exfft	1280K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	1280K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	1280K, SFX_P4_CLM20, 11
	exfft	1280K, SFX_P4_CLM4 + SFX_AMD64_CLM1, 12
	exfft	1536K, SFX_PPRO + SFX_P3
allfft	exfft	1536K, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	1536K, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 11
allfft	exfft	1536K9, SFX_P4_CLM10, 11
allfft	exfft	1536K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	1536K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	1536K, SFX_P4_CLM10, 11
	exfft	1536K, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	1792K, SFX_PPRO + SFX_P3
allfft	exfft	1792K, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	1792K, SFX_P4_CLM4210 + SFX_AMD64_CLM4210, 11
allfft	exfft	1792K9, SFX_P4_CLM10, 11
allfft	exfft	1792K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	1792K9, SFX_P4_CLM10, 12
allfft	exfft	1792K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	1792K, SFX_P4_CLM10, 11
	exfft	1792K, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	2048K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	2048K, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	2048K, SFX_P4_CLM4210 + SFX_AMD64_CLM4210, 11
allfft	exfft	2048K9, SFX_P4_CLM10, 11
allfft	exfft	2048K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	2048K9, SFX_P4_CLM10, 12
allfft	exfft	2048K, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	2048K, SFX_P4_CLM0, 11
	exfft	2048K, SFX_P4_CLM421 + SFX_AMD64_CLM1, 12
	exfft	2560K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	2560K, SFX_P4_CLM0, 10
allfft	exfft	2560K, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	2560K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	2560K, SFX_P4_CLM4 + SFX_AMD64_CLM421, 13
	exfft	2560K, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	2560K, SFX_AMD64_CLM1, 13
	exfft	3072K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	3072K, SFX_P4_CLM0, 10
allfft	exfft	3072K, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	3072K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	3072K, SFX_P4_CLM4 + SFX_AMD64_CLM421, 13
	exfft	3072K, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	3072K, SFX_AMD64_CLM1, 13
	exfft	3584K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	3584K, SFX_P4_CLM0, 10
allfft	exfft	3584K, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	3584K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	3584K, SFX_P4_CLM4 + SFX_AMD64_CLM421, 13
	exfft	3584K, SFX_P4_CLM421 + SFX_AMD64_CLM1, 12
	exfft	3584K, SFX_AMD64_CLM1, 13
	exfft	4096K, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	4096K, SFX_P4_CLM0, 10
allfft	exfft	4096K, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	4096K, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	4096K, SFX_P4_CLM42 + SFX_AMD64_CLM421, 13
	exfft	4096K, SFX_P4_CLM421 + SFX_AMD64_CLM1, 12
	exfft	4096K, SFX_AMD64_CLM1, 13
allfft	exfft	5M, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	5M, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	5M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	5M, SFX_P4_CLM41 + SFX_AMD64_CLM41, 12
allfft	exfft	6M, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	6M, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	6M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	6M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 12
allfft	exfft	7M, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	7M, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	7M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	7M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 12
allfft	exfft	8M, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	8M, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	8M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	8M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 12
allfft	exfft	10M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	10M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	10M, SFX_P4_CLM41 + SFX_AMD64_CLM41, 13
allfft	exfft	12M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	12M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	12M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 13
allfft	exfft	14M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	14M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	14M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 13
allfft	exfft	16M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	16M, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	16M, SFX_P4_CLM421 + SFX_AMD64_CLM41, 13
allfft	exfft	20M, SFX_P4_CLM210 + SFX_AMD64_CLM4210, 13
	exfft	20M, SFX_P4_CLM21 + SFX_AMD64_CLM20, 13
allfft	exfft	24M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 13
	exfft	24M, SFX_P4_CLM210 + SFX_AMD64_CLM20, 13
allfft	exfft	28M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 13
	exfft	28M, SFX_P4_CLM210 + SFX_AMD64_CLM20, 13
allfft	exfft	32M, SFX_P4_CLM210 + SFX_AMD64_CLM210, 13
	exfft	32M, SFX_P4_CLM210 + SFX_AMD64_CLM10, 13

	exfft	32p, SFX_PPRO + SFX_P4
	exfft	48p, SFX_PPRO + SFX_P4
	exfft	64p, SFX_PPRO + SFX_P4
	exfft	96p, SFX_PPRO + SFX_P4
	exfft	128p, SFX_PPRO + SFX_P4
	exfft	192p, SFX_PPRO + SFX_P4
	exfft	256p, SFX_PPRO + SFX_P4
	exfft	384p, SFX_PPRO + SFX_P4
	exfft	512p, SFX_PPRO + SFX_P4
	exfft	768p, SFX_PPRO + SFX_P4
	exfft	1024p, SFX_PPRO + SFX_P4
	exfft	1536p, SFX_PPRO + SFX_P4
	exfft	2048p, SFX_PPRO + SFX_P4
	exfft	3072p, SFX_PPRO + SFX_P4
	exfft	4096p, SFX_PPRO + SFX_P4
	exfft	6144p, SFX_PPRO + SFX_P4
	exfft	8192p, SFX_PPRO + SFX_P4
	exfft	12Kp, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	16Kp, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	24Kp, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	32Kp, SFX_PPRO + SFX_P4 + SFX_AMD64
	exfft	48Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	64Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	96Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	128Kp, SFX_PPRO + SFX_P3 + SFX_P4 + SFX_AMD64
	exfft	192Kp, SFX_PPRO + SFX_P3
allfft	exfft	192Kp, SFX_P4 + SFX_AMD64
allfft	exfft	192Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	192Kp, SFX_P4_CLM1 + SFX_AMD64_CLM2, 10
	exfft	256Kp, SFX_PPRO + SFX_P3
allfft	exfft	256Kp, SFX_P4 + SFX_AMD64
allfft	exfft	256Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
	exfft	256Kp, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	384Kp, SFX_PPRO + SFX_P3
allfft	exfft	384Kp, SFX_P4 + SFX_AMD64
allfft	exfft	384Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 10
allfft	exfft	384Kp, SFX_P4_CLM4 + SFX_AMD64_CLM4, 11
	exfft	384Kp, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	512Kp, SFX_PPRO + SFX_P3
allfft	exfft	512Kp, SFX_P4 + SFX_AMD64
allfft	exfft	512Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 10
allfft	exfft	512Kp, SFX_P4_CLM4 + SFX_AMD64_CLM4, 11
	exfft	512Kp, SFX_P4_CLM42 + SFX_AMD64_CLM1, 10
	exfft	768Kp, SFX_PPRO + SFX_P3
allfft	exfft	768Kp, SFX_P4_CLM1, 8
allfft	exfft	768Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 10
allfft	exfft	768Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 11
allfft	exfft	768Kp9, SFX_P4_CLM10, 11
allfft	exfft	768Kp, SFX_P4_CLM42 + SFX_AMD64_CLM421, 12
	exfft	768Kp, SFX_P4_CLM420, 11
	exfft	768Kp, SFX_P4_CLM2 + SFX_AMD64_CLM2, 12
	exfft	1024Kp, SFX_PPRO + SFX_P3
allfft	exfft	1024Kp, SFX_P4_CLM1, 8
allfft	exfft	1024Kp, SFX_P4_CLM210 + SFX_AMD64_CLM2, 10
allfft	exfft	1024Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM42, 11
allfft	exfft	1024Kp9, SFX_P4_CLM10, 11
allfft	exfft	1024Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
	exfft	1024Kp, SFX_P4_CLM420, 11
	exfft	1024Kp, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	1536Kp, SFX_PPRO + SFX_P3
allfft	exfft	1536Kp, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	1536Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 11
allfft	exfft	1536Kp9, SFX_P4_CLM10, 11
allfft	exfft	1536Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	1536Kp, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	1536Kp, SFX_P4_CLM10, 11
	exfft	1536Kp, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	2048Kp, SFX_PPRO + SFX_P3
allfft	exfft	2048Kp, SFX_P4_CLM10 + SFX_AMD64_CLM10, 10
allfft	exfft	2048Kp, SFX_P4_CLM4210 + SFX_AMD64_CLM4210, 11
allfft	exfft	2048Kp9, SFX_P4_CLM10, 11
allfft	exfft	2048Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	2048Kp9, SFX_P4_CLM10, 12
allfft	exfft	2048Kp, SFX_P4_CLM4 + SFX_AMD64_CLM4, 13
	exfft	2048Kp, SFX_P4_CLM0, 11
	exfft	2048Kp, SFX_P4_CLM41 + SFX_AMD64_CLM1, 12
	exfft	3072Kp, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	3072Kp, SFX_P4_CLM0, 10
allfft	exfft	3072Kp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	3072Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	3072Kp, SFX_P4_CLM4 + SFX_AMD64_CLM421, 13
	exfft	3072Kp, SFX_P4_CLM421 + SFX_AMD64_CLM1, 12
	exfft	3072Kp, SFX_AMD64_CLM1, 13
	exfft	4096Kp, SFX_PPRO_CLM2 + SFX_P3_CLM2
allfft	exfft	4096Kp, SFX_P4_CLM0, 10
allfft	exfft	4096Kp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 11
allfft	exfft	4096Kp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 12
allfft	exfft	4096Kp, SFX_P4_CLM42 + SFX_AMD64_CLM421, 13
	exfft	4096Kp, SFX_P4_CLM421 + SFX_AMD64_CLM1, 12
	exfft	4096Kp, SFX_AMD64_CLM1, 13
allfft	exfft	6Mp, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	6Mp, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	6Mp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	6Mp, SFX_P4_CLM421 + SFX_AMD64_CLM41, 12
allfft	exfft	8Mp, SFX_P4_CLM10 + SFX_AMD64_CLM10, 11
allfft	exfft	8Mp, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 12
allfft	exfft	8Mp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	8Mp, SFX_P4_CLM421 + SFX_AMD64_CLM41, 12
allfft	exfft	12Mp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	12Mp, SFX_P4_CLM421 + SFX_AMD64_CLM421, 13
	exfft	12Mp, SFX_P4_CLM21 + SFX_AMD64_CLM41, 13
allfft	exfft	16Mp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 12
allfft	exfft	16Mp, SFX_P4_CLM4210 + SFX_AMD64_CLM421, 13
	exfft	16Mp, SFX_P4_CLM21 + SFX_AMD64_CLM41, 13
allfft	exfft	24Mp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 13
	exfft	24Mp, SFX_P4_CLM210 + SFX_AMD64_CLM20, 13
allfft	exfft	32Mp, SFX_P4_CLM210 + SFX_AMD64_CLM210, 13
	exfft	32Mp, SFX_P4_CLM210 + SFX_AMD64_CLM10, 13

;
; Global variables needed by FFT setup code
;

IFNDEF X86_64
_GWDATA SEGMENT PAGE PUBLIC 'DATA'
ELSE
_GWDATA SEGMENT PAGE
ENDIF

; Jmptable definitions

RPF	EQU	80000000h	; FFT requires prefetch capability
RPFW	EQU	40000000h	; FFT requires prefetchw (3DNow!) capability

IFNDEF	X86_64
jmptable DD	755,	32,	0.0000036,	672
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft32_1PPRO, OFFSET fft32_2PPRO
	DD			OFFSET fft32_3PPRO, OFFSET fft32_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			3, 1, 1, 1, 0
	DD	939,	40,	0.0000057,	948
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft40_1PPRO, OFFSET fft40_2PPRO
	DD			OFFSET fft40_3PPRO, OFFSET fft40_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			4, 1, 1, 1, 0
	DD	1113,	48,	0.0000065,	1128
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48_1PPRO, OFFSET fft48_2PPRO
	DD			OFFSET fft48_3PPRO, OFFSET fft48_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			5, 1, 1, 1, 0
	DD	1303,	56,	0.0000084,	1356
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft56_1PPRO, OFFSET fft56_2PPRO
	DD			OFFSET fft56_3PPRO, OFFSET fft56_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			6, 1, 1, 1, 0
	DD	1499,	64,	0.0000083,	1392
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64_1PPRO, OFFSET fft64_2PPRO
	DD			OFFSET fft64_3PPRO, OFFSET fft64_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			7, 1, 1, 1, 0
	DD	1857,	80,	0.0000121,	1848
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft80_1PPRO, OFFSET fft80_2PPRO
	DD			OFFSET fft80_3PPRO, OFFSET fft80_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			9, 2, 1, 1, 0
	DD	2211,	96,	0.0000141,	2208
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96_1PPRO, OFFSET fft96_2PPRO
	DD			OFFSET fft96_3PPRO, OFFSET fft96_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			11, 2, 1, 1, 0
	DD	2585,	112,	0.0000179,	2616
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft112_1PPRO, OFFSET fft112_2PPRO
	DD			OFFSET fft112_3PPRO, OFFSET fft112_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			13, 3, 1, 1, 0
	DD	2953,	128,	0.0000178,	2832
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128_1PPRO, OFFSET fft128_2PPRO
	DD			OFFSET fft128_3PPRO, OFFSET fft128_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			15, 3, 1, 1, 0
	DD	3663,	160,	0.0000296,	3840
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft160_1PPRO, OFFSET fft160_2PPRO
	DD			OFFSET fft160_3PPRO, OFFSET fft160_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			19, 4, 2, 1, 0
	DD	4359,	192,	0.000035,	4608
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft192_1PPRO, OFFSET fft192_2PPRO
	DD			OFFSET fft192_3PPRO, OFFSET fft192_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			23, 5, 2, 1, 0
	DD	5093,	224,	0.000045,	5424
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft224_1PPRO, OFFSET fft224_2PPRO
	DD			OFFSET fft224_3PPRO, OFFSET fft224_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			27, 6, 3, 1, 0
	DD	5833,	256,	0.000045,	5712
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft256_1PPRO, OFFSET fft256_2PPRO
	DD			OFFSET fft256_3PPRO, OFFSET fft256_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			31, 7, 3, 1, 0
	DD	7243,	320,	0.000062,	7680
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft320_1PPRO, OFFSET fft320_2PPRO
	DD			OFFSET fft320_3PPRO, OFFSET fft320_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			39, 9, 2, 1, 0
	DD	8639,	384,	0.000075,	9216
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft384_1PPRO, OFFSET fft384_2PPRO
	DD			OFFSET fft384_3PPRO, OFFSET fft384_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			47, 11, 2, 1, 0
	DD	10085,	448,	0.000093,	10800
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft448_1PPRO, OFFSET fft448_2PPRO
	DD			OFFSET fft448_3PPRO, OFFSET fft448_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			55, 13, 3, 1, 0
	DD	11537,	512,	0.000097,	11472
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft512_1PPRO, OFFSET fft512_2PPRO
	DD			OFFSET fft512_3PPRO, OFFSET fft512_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			63, 15, 3, 1, 0
	DD	14301,	640,	0.000140,	15552
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft640_1PPRO, OFFSET fft640_2PPRO
	DD			OFFSET fft640_3PPRO, OFFSET fft640_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			79, 19, 4, 2, 0
	DD	17047,	768,	0.000167,	18672
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft768_1PPRO, OFFSET fft768_2PPRO
	DD			OFFSET fft768_3PPRO, OFFSET fft768_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			95, 23, 5, 2, 0
	DD	19881,	896,	0.000210,	21840
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft896_1PPRO, OFFSET fft896_2PPRO
	DD			OFFSET fft896_3PPRO, OFFSET fft896_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			111, 27, 6, 3, 0
	DD	22799,	1024,	0.000218,	22992
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1024_1PPRO, OFFSET fft1024_2PPRO
	DD			OFFSET fft1024_3PPRO, OFFSET fft1024_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			127, 31, 7, 3, 0
	DD	28295,	1280,	0.000302,	31152
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1280_1PPRO, OFFSET fft1280_2PPRO
	DD			OFFSET fft1280_3PPRO, OFFSET fft1280_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			159, 39, 9, 2, 0
	DD	33761,	1536,	0.000365,	37392
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1536_1PPRO, OFFSET fft1536_2PPRO
	DD			OFFSET fft1536_3PPRO, OFFSET fft1536_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			191, 47, 11, 2, 0
	DD	39411,	1792,	0.000456,	43680
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1792_1PPRO, OFFSET fft1792_2PPRO
	DD			OFFSET fft1792_3PPRO, OFFSET fft1792_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			223, 55, 13, 3, 0
	DD	45061,	2048,	0.000490,	46032
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft2048_1PPRO, OFFSET fft2048_2PPRO
	DD			OFFSET fft2048_3PPRO, OFFSET fft2048_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			255, 63, 15, 3, 0
	DD	55825,	2560,	0.000708,	62544
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft2560_1PPRO, OFFSET fft2560_2PPRO
	DD			OFFSET fft2560_3PPRO, OFFSET fft2560_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			319, 79, 19, 4, 2, 0
	DD	66519,	3072,	0.000851,	75072
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft3072_1PPRO, OFFSET fft3072_2PPRO
	DD			OFFSET fft3072_3PPRO, OFFSET fft3072_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			383, 95, 23, 5, 2, 0
	DD	77599,	3584,	0.00107,	87648
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft3584_1PPRO, OFFSET fft3584_2PPRO
	DD			OFFSET fft3584_3PPRO, OFFSET fft3584_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			447, 111, 27, 6, 3, 0
	DD	89047,	4096,	0.00113,	92112
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft4096_1PPRO, OFFSET fft4096_2PPRO
	DD			OFFSET fft4096_3PPRO, OFFSET fft4096_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			511, 127, 31, 7, 3, 0
	DD	110400,	5120,	0.00152,	24848
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft5120_1PPRO, OFFSET fft5120_2PPRO
	DD			OFFSET fft5120_3PPRO, OFFSET fft5120_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			10, 1
	DD			9, 1, 1, 1, 0
	DD	131100,	6144,	0.00191,	28016
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft6144_1PPRO, OFFSET fft6144_2PPRO
	DD			OFFSET fft6144_3PPRO, OFFSET fft6144_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			12, 1
	DD			11, 1, 1, 1, 0
	DD	152800,	7168,	0.00226,	31232
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft7168_1PPRO, OFFSET fft7168_2PPRO
	DD			OFFSET fft7168_3PPRO, OFFSET fft7168_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			14, 1
	DD			13, 1, 1, 1, 0
	DD	175300,	8192,	0.00242,	34400
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft8192_1PPRO, OFFSET fft8192_2PPRO
	DD			OFFSET fft8192_3PPRO, OFFSET fft8192_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			16, 1
	DD			15, 1, 1, 1, 0
	DD	217700,	10240,	0.00333,	40880
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft10K_1PPRO, OFFSET fft10K_2PPRO
	DD			OFFSET fft10K_3PPRO, OFFSET fft10K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			20, 1
	DD			19, 1, 1, 1, 0
	DD	258200,	12288,	0.00397,	47264
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft12K_1PPRO, OFFSET fft12K_2PPRO
	DD			OFFSET fft12K_3PPRO, OFFSET fft12K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			24, 1
	DD			23, 1, 1, 1, 0
	DD	301400,	14336,	0.00488,	53696
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft14K_1PPRO, OFFSET fft14K_2PPRO
	DD			OFFSET fft14K_3PPRO, OFFSET fft14K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			28, 1
	DD			27, 1, 1, 1, 0
	DD	346100,	16384,	0.00522,	59936
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft16K_1PPRO, OFFSET fft16K_2PPRO
	DD			OFFSET fft16K_3PPRO, OFFSET fft16K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			32, 1
	DD			31, 1, 1, 1, 0
	DD	430300,	20480,	0.00692,	72800
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft20K_1PPRO, OFFSET fft20K_2PPRO
	DD			OFFSET fft20K_3PPRO, OFFSET fft20K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			40, 1
	DD			39, 1, 1, 1, 0
	DD	511600,	24576,	0.00826,	85568
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft24K_1PPRO, OFFSET fft24K_2PPRO
	DD			OFFSET fft24K_3PPRO, OFFSET fft24K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			48, 1
	DD			47, 1, 1, 1, 0
	DD	596100,	28672,	0.0101,		98384
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft28K_1PPRO, OFFSET fft28K_2PPRO
	DD			OFFSET fft28K_3PPRO, OFFSET fft28K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			56, 1
	DD			55, 1, 1, 1, 0
	DD	683700,	32768,	0.0109,		111008
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft32K_1PPRO, OFFSET fft32K_2PPRO
	DD			OFFSET fft32K_3PPRO, OFFSET fft32K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			64, 1
	DD			63, 1, 1, 1, 0
	DD	848800,	40960,	0.0151,		136832
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft40K_1P3, OFFSET fft40K_2P3
	DD			OFFSET fft40K_3P3, OFFSET fft40K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			80, 1
	DD			79, 1, 1, 1, 0
	DD	848800,	40960,	0.0151,		136832
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft40K_1PPRO, OFFSET fft40K_2PPRO
	DD			OFFSET fft40K_3PPRO, OFFSET fft40K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			80, 1
	DD			79, 1, 1, 1, 0
	DD	1009000, 49152,	0.0184,		162416
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48K_1P3, OFFSET fft48K_2P3
	DD			OFFSET fft48K_3P3, OFFSET fft48K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			96, 1
	DD			95, 1, 1, 1, 0
	DD	1009000, 49152,	0.0184,		162416
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48K_1PPRO, OFFSET fft48K_2PPRO
	DD			OFFSET fft48K_3PPRO, OFFSET fft48K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			96, 1
	DD			95, 1, 1, 1, 0
	DD	1177000, 57344,	0.0227,		188048
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft56K_1P3, OFFSET fft56K_2P3
	DD			OFFSET fft56K_3P3, OFFSET fft56K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			112, 1
	DD			111, 1, 1, 1, 0
	DD	1177000, 57344,	0.0227,		188048
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft56K_1PPRO, OFFSET fft56K_2PPRO
	DD			OFFSET fft56K_3PPRO, OFFSET fft56K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			112, 1
	DD			111, 1, 1, 1, 0
	DD	1350000, 65536,	0.0252,		213152
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64K_1P3, OFFSET fft64K_2P3
	DD			OFFSET fft64K_3P3, OFFSET fft64K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			128, 1
	DD			127, 1, 1, 1, 0
	DD	1350000, 65536,	0.0252,		213152
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64K_1PPRO, OFFSET fft64K_2PPRO
	DD			OFFSET fft64K_3PPRO, OFFSET fft64K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			128, 1
	DD			127, 1, 1, 1, 0
	DD	1678000, 81920, 0.0360,		264752
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft80K_1P3, OFFSET fft80K_2P3
	DD			OFFSET fft80K_3P3, OFFSET fft80K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	1678000, 81920, 0.0360,		264752
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft80K_1PPRO, OFFSET fft80K_2PPRO
	DD			OFFSET fft80K_3PPRO, OFFSET fft80K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	1994000, 98304, 0.0445,		297536
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96K_1P3, OFFSET fft96K_2P3
	DD			OFFSET fft96K_3P3, OFFSET fft96K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			48, 1
	DD			47, 1, 1, 1, 0
	DD	1994000, 98304, 0.0445,		297536
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96K_1PPRO, OFFSET fft96K_2PPRO
	DD			OFFSET fft96K_3PPRO, OFFSET fft96K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			48, 1
	DD			47, 1, 1, 1, 0
	DD	2324000, 114688, 0.0548,	341072
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft112K_1P3, OFFSET fft112K_2P3
	DD			OFFSET fft112K_3P3, OFFSET fft112K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			56, 1
	DD			55, 1, 1, 1, 0
	DD	2324000, 114688, 0.0548,	341072
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft112K_1PPRO, OFFSET fft112K_2PPRO
	DD			OFFSET fft112K_3PPRO, OFFSET fft112K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			56, 1
	DD			55, 1, 1, 1, 0
	DD	2664000, 131072, 0.0604,	384416
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128K_1P3, OFFSET fft128K_2P3
	DD			OFFSET fft128K_3P3, OFFSET fft128K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			64, 1
	DD			63, 1, 1, 1, 0
	DD	2664000, 131072, 0.0604,	384416
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128K_1PPRO, OFFSET fft128K_2PPRO
	DD			OFFSET fft128K_3PPRO, OFFSET fft128K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			64, 1
	DD			63, 1, 1, 1, 0
	DD	3310000, 163840, 0.0830,	471680
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			5120		;; scratch area size
	DD			OFFSET fft160K_1P3, OFFSET fft160K_2P3
	DD			OFFSET fft160K_3P3, OFFSET fft160K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			80, 1
	DD			79, 1, 1, 1, 0
	DD	3310000, 163840, 0.0830,	471680
	DD			2		;; Flags, min_l2_cache, clm
	DD			5120		;; scratch area size
	DD			OFFSET fft160K_1PPRO, OFFSET fft160K_2PPRO
	DD			OFFSET fft160K_3PPRO, OFFSET fft160K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			80, 1
	DD			79, 1, 1, 1, 0
	DD	3933000, 196608, 0.0982,	558704
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			OFFSET fft192K_1P3, OFFSET fft192K_2P3
	DD			OFFSET fft192K_3P3, OFFSET fft192K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			96, 1
	DD			95, 1, 1, 1, 0
	DD	3933000, 196608, 0.0982,	558704
	DD			2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			OFFSET fft192K_1PPRO, OFFSET fft192K_2PPRO
	DD			OFFSET fft192K_3PPRO, OFFSET fft192K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			96, 1
	DD			95, 1, 1, 1, 0
	DD	4593000, 229376, 0.1193,	645776
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			7168		;; scratch area size
	DD			OFFSET fft224K_1P3, OFFSET fft224K_2P3
	DD			OFFSET fft224K_3P3, OFFSET fft224K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			112, 1
	DD			111, 1, 1, 1, 0
	DD	4593000, 229376, 0.1193,	645776
	DD			2		;; Flags, min_l2_cache, clm
	DD			7168		;; scratch area size
	DD			OFFSET fft224K_1PPRO, OFFSET fft224K_2PPRO
	DD			OFFSET fft224K_3PPRO, OFFSET fft224K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			112, 1
	DD			111, 1, 1, 1, 0
	DD	5264000, 262144, 0.1316,	732320
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			OFFSET fft256K_1P3, OFFSET fft256K_2P3
	DD			OFFSET fft256K_3P3, OFFSET fft256K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			128, 1
	DD			127, 1, 1, 1, 0
	DD	5264000, 262144, 0.1316,	732320
	DD			2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			OFFSET fft256K_1PPRO, OFFSET fft256K_2PPRO
	DD			OFFSET fft256K_3PPRO, OFFSET fft256K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			128, 1
	DD			127, 1, 1, 1, 0
	DD	6545000, 327680, 0.1726,	906800
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			OFFSET fft320K_1P3, OFFSET fft320K_2P3
	DD			OFFSET fft320K_3P3, OFFSET fft320K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	6545000, 327680, 0.1726,	906800
	DD			2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			OFFSET fft320K_1PPRO, OFFSET fft320K_2PPRO
	DD			OFFSET fft320K_3PPRO, OFFSET fft320K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	7772000, 393216, 0.2107,	1080848
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft384K_1P3, OFFSET fft384K_2P3
	DD			OFFSET fft384K_3P3, OFFSET fft384K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			192, 1
	DD			191, 1, 1, 1, 0
	DD	7772000, 393216, 0.2107,	1080848
	DD			2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft384K_1PPRO, OFFSET fft384K_2PPRO
	DD			OFFSET fft384K_3PPRO, OFFSET fft384K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			192, 1
	DD			191, 1, 1, 1, 0
	DD	9071000, 458752, 0.2520,	1254944
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			OFFSET fft448K_1P3, OFFSET fft448K_2P3
	DD			OFFSET fft448K_3P3, OFFSET fft448K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			224, 1
	DD			223, 1, 1, 1, 0
	DD	9071000, 458752, 0.2520,	1254944
	DD			2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			OFFSET fft448K_1PPRO, OFFSET fft448K_2PPRO
	DD			OFFSET fft448K_3PPRO, OFFSET fft448K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			224, 1
	DD			223, 1, 1, 1, 0
	DD	10380000, 524288, 0.2808,	1428128
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft512K_1P3, OFFSET fft512K_2P3
	DD			OFFSET fft512K_3P3, OFFSET fft512K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			256, 1
	DD			255, 1, 1, 1, 0
	DD	10380000, 524288, 0.2808,	1428128
	DD			2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft512K_1PPRO, OFFSET fft512K_2PPRO
	DD			OFFSET fft512K_3PPRO, OFFSET fft512K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			256, 1
	DD			255, 1, 1, 1, 0
	DD	12890000, 655360, 0.372,	1777232
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			OFFSET fft640K_1P3, OFFSET fft640K_2P3
	DD			OFFSET fft640K_3P3, OFFSET fft640K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			320, 1
	DD			319, 1, 1, 1, 0
	DD	12890000, 655360, 0.372,	1777232
	DD			2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			OFFSET fft640K_1PPRO, OFFSET fft640K_2PPRO
	DD			OFFSET fft640K_3PPRO, OFFSET fft640K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			320, 1
	DD			319, 1, 1, 1, 0
	DD	15310000, 786432, 0.453,	2125376
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft768K_1P3, OFFSET fft768K_2P3
	DD			OFFSET fft768K_3P3, OFFSET fft768K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			384, 1
	DD			383, 1, 1, 1, 0
	DD	15310000, 786432, 0.453,	2125376
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft768K_1PPRO, OFFSET fft768K_2PPRO
	DD			OFFSET fft768K_3PPRO, OFFSET fft768K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			384, 1
	DD			383, 1, 1, 1, 0
	DD	17890000, 917504, 0.536,	2473568
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			OFFSET fft896K_1P3, OFFSET fft896K_2P3
	DD			OFFSET fft896K_3P3, OFFSET fft896K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			448, 1
	DD			447, 1, 1, 1, 0
	DD	17890000, 917504, 0.536,	2473568
	DD			2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			OFFSET fft896K_1PPRO, OFFSET fft896K_2PPRO
	DD			OFFSET fft896K_3PPRO, OFFSET fft896K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			448, 1
	DD			447, 1, 1, 1, 0
	DD	20460000, 1048576, 0.600,	2819744
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft1024K2_1P3, OFFSET fft1024K2_2P3
	DD			OFFSET fft1024K2_3P3, OFFSET fft1024K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			512, 1
	DD			511, 1, 1, 1, 0
	DD	20460000, 1048576, 0.600,	2819744
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft1024K2_1PPRO, OFFSET fft1024K2_2PPRO
	DD			OFFSET fft1024K2_3PPRO, OFFSET fft1024K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			512, 1
	DD			511, 1, 1, 1, 0
	DD	25390000, 1310720, 0.776,	3474992
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			OFFSET fft1280K_1P3, OFFSET fft1280K_2P3
	DD			OFFSET fft1280K_3P3, OFFSET fft1280K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	25390000, 1310720, 0.776,	3474992
	DD			2		;; Flags, min_l2_cache, clm
	DD			10240+4*64	;; scratch area size
	DD			OFFSET fft1280K_1PPRO, OFFSET fft1280K_2PPRO
	DD			OFFSET fft1280K_3PPRO, OFFSET fft1280K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			160, 1
	DD			159, 1, 1, 1, 0
	DD	30190000, 1572864, 0.934,	4140560
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft1536K_1P3, OFFSET fft1536K_2P3
	DD			OFFSET fft1536K_3P3, OFFSET fft1536K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			192, 1
	DD			191, 1, 1, 1, 0
	DD	30190000, 1572864, 0.934,	4140560
	DD			2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft1536K_1PPRO, OFFSET fft1536K_2PPRO
	DD			OFFSET fft1536K_3PPRO, OFFSET fft1536K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			192, 1
	DD			191, 1, 1, 1, 0
	DD	35200000, 1835008, 1.113,	4806176
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			OFFSET fft1792K_1P3, OFFSET fft1792K_2P3
	DD			OFFSET fft1792K_3P3, OFFSET fft1792K_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			224, 1
	DD			223, 1, 1, 1, 0
	DD	35200000, 1835008, 1.113,	4806176
	DD			2		;; Flags, min_l2_cache, clm
	DD			14336+6*64	;; scratch area size
	DD			OFFSET fft1792K_1PPRO, OFFSET fft1792K_2PPRO
	DD			OFFSET fft1792K_3PPRO, OFFSET fft1792K_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			224, 1
	DD			223, 1, 1, 1, 0
	DD	40300000, 2097152, 1.226,	5470880
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft2048K2_1P3, OFFSET fft2048K2_2P3
	DD			OFFSET fft2048K2_3P3, OFFSET fft2048K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, 1
	DD			255, 1, 1, 1, 0
	DD	40300000, 2097152, 1.226,	5470880
	DD			2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft2048K2_1PPRO, OFFSET fft2048K2_2PPRO
	DD			OFFSET fft2048K2_3PPRO, OFFSET fft2048K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, 1
	DD			255, 1, 1, 1, 0
	DD	50020000, 2621440, 1.636,	6803024
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			OFFSET fft2560K2_1P3, OFFSET fft2560K2_2P3
	DD			OFFSET fft2560K2_3P3, OFFSET fft2560K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			320, 1
	DD			319, 1, 1, 1, 0
	DD	50020000, 2621440, 1.636,	6803024
	DD			2		;; Flags, min_l2_cache, clm
	DD			20480+9*64	;; scratch area size
	DD			OFFSET fft2560K2_1PPRO, OFFSET fft2560K2_2PPRO
	DD			OFFSET fft2560K2_3PPRO, OFFSET fft2560K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			320, 1
	DD			319, 1, 1, 1, 0
	DD	59510000, 3145728, 1.990,	8134208
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft3072K2_1P3, OFFSET fft3072K2_2P3
	DD			OFFSET fft3072K2_3P3, OFFSET fft3072K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, 1
	DD			383, 1, 1, 1, 0
	DD	59510000, 3145728, 1.990,	8134208
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft3072K2_1PPRO, OFFSET fft3072K2_2PPRO
	DD			OFFSET fft3072K2_3PPRO, OFFSET fft3072K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, 1
	DD			383, 1, 1, 1, 0
	DD	69360000, 3670016, 2.380,	9465440
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			OFFSET fft3584K2_1P3, OFFSET fft3584K2_2P3
	DD			OFFSET fft3584K2_3P3, OFFSET fft3584K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			448, 1
	DD			447, 1, 1, 1, 0
	DD	69360000, 3670016, 2.380,	9465440
	DD			2		;; Flags, min_l2_cache, clm
	DD			28672+13*64	;; scratch area size
	DD			OFFSET fft3584K2_1PPRO, OFFSET fft3584K2_2PPRO
	DD			OFFSET fft3584K2_3PPRO, OFFSET fft3584K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			448, 1
	DD			447, 1, 1, 1, 0
	DD	79370000, 4194304, 2.604,	10794656
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft4096K2_1P3, OFFSET fft4096K2_2P3
	DD			OFFSET fft4096K2_3P3, OFFSET fft4096K2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, 1
	DD			511, 1, 1, 1, 0
	DD	79370000, 4194304, 2.604,	10794656
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft4096K2_1PPRO, OFFSET fft4096K2_2PPRO
	DD			OFFSET fft4096K2_3PPRO, OFFSET fft4096K2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, 1
	DD			511, 1, 1, 1, 0
	DD	0
jmptablep DD	755,	32,	0.000004,	976
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft32p_1PPRO, OFFSET fft32p_2PPRO
	DD			OFFSET fft32p_3PPRO, OFFSET fft32p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			4, 1, 1, 1, 0
	DD	1111,	48,	0.000007,	1608
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48p_1PPRO, OFFSET fft48p_2PPRO
	DD			OFFSET fft48p_3PPRO, OFFSET fft48p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			6, 3, 1, 1, 0
	DD	1485,	64,	0.000010,	1952
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64p_1PPRO, OFFSET fft64p_2PPRO
	DD			OFFSET fft64p_3PPRO, OFFSET fft64p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			8, 4, 1, 1, 0
	DD	2199,	96,	0.0000141,	3072
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96p_1PPRO, OFFSET fft96p_2PPRO
	DD			OFFSET fft96p_3PPRO, OFFSET fft96p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			12, 3, 1, 1, 0
	DD	2947,	128,	0.000021,	3904
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128p_1PPRO, OFFSET fft128p_2PPRO
	DD			OFFSET fft128p_3PPRO, OFFSET fft128p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			16, 4, 1, 1, 1, 0
	DD	4345,	192,	0.000035,	6288
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft192p_1PPRO, OFFSET fft192p_2PPRO
	DD			OFFSET fft192p_3PPRO, OFFSET fft192p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			24, 6, 3, 1, 0
	DD	5817,	256,	0.000051,	7808
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft256p_1PPRO, OFFSET fft256p_2PPRO
	DD			OFFSET fft256p_3PPRO, OFFSET fft256p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			32, 8, 4, 1, 1, 0
	DD	8607,	384,	0.000075,	12432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft384p_1PPRO, OFFSET fft384p_2PPRO
	DD			OFFSET fft384p_3PPRO, OFFSET fft384p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			48, 12, 3, 1, 0
	DD	11515,	512,	0.000106,	15616
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft512p_1PPRO, OFFSET fft512p_2PPRO
	DD			OFFSET fft512p_3PPRO, OFFSET fft512p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			64, 16, 4, 1, 1, 0
	DD	17001,	768,	0.000167,	25008
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft768p_1PPRO, OFFSET fft768p_2PPRO
	DD			OFFSET fft768p_3PPRO, OFFSET fft768p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			96, 24, 6, 3, 1, 0
	DD	22701,	1024,	0.000249,	31232
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1024p_1PPRO, OFFSET fft1024p_2PPRO
	DD			OFFSET fft1024p_3PPRO, OFFSET fft1024p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			128, 32, 8, 4, 1, 0
	DD	33569,	1536,	0.000365,	49872
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft1536p_1PPRO, OFFSET fft1536p_2PPRO
	DD			OFFSET fft1536p_3PPRO, OFFSET fft1536p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			192, 48, 12, 3, 0
	DD	44951,	2048,	0.000582,	62464
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft2048p_1PPRO, OFFSET fft2048p_2PPRO
	DD			OFFSET fft2048p_3PPRO, OFFSET fft2048p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			256, 64, 16, 4, 1, 0
	DD	66319,	3072,	0.000851,	99888
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft3072p_1PPRO, OFFSET fft3072p_2PPRO
	DD			OFFSET fft3072p_3PPRO, OFFSET fft3072p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			384, 96, 24, 6, 3, 0
	DD	88747,	4096,	0.00135,	124928
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft4096p_1PPRO, OFFSET fft4096p_2PPRO
	DD			OFFSET fft4096p_3PPRO, OFFSET fft4096p_4PPRO
	DD			OFFSET prctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1, 1
	DD			512, 128, 32, 8, 4, 0
	DD	130600,	6144,	0.00191,	26512
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft6144p_1PPRO, OFFSET fft6144p_2PPRO
	DD			OFFSET fft6144p_3PPRO, OFFSET fft6144p_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			12, 1
	DD			12, 1, 1, 1, 0
	DD	174000,	8192,	0.00284,	32960
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft8192p_1PPRO, OFFSET fft8192p_2PPRO
	DD			OFFSET fft8192p_3PPRO, OFFSET fft8192p_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			16, 1
	DD			16, 1, 1, 1, 0
	DD	257700,	12288,	0.00397,	46000
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft12Kp_1PPRO, OFFSET fft12Kp_2PPRO
	DD			OFFSET fft12Kp_3PPRO, OFFSET fft12Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			24, 1
	DD			24, 1, 1, 1, 0
	DD	344700,	16384,	0.00588,	58752
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft16Kp_1PPRO, OFFSET fft16Kp_2PPRO
	DD			OFFSET fft16Kp_3PPRO, OFFSET fft16Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			32, 1
	DD			32, 1, 1, 1, 0
	DD	508600,	24576,	0.00826,	84688
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft24Kp_1PPRO, OFFSET fft24Kp_2PPRO
	DD			OFFSET fft24Kp_3PPRO, OFFSET fft24Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			48, 1
	DD			48, 1, 1, 1, 0
	DD	679400,	32768,	0.01299,	110336
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft32Kp_1PPRO, OFFSET fft32Kp_2PPRO
	DD			OFFSET fft32Kp_3PPRO, OFFSET fft32Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			64, 1
	DD			64, 1, 1, 1, 0
	DD	1006000, 49152,	0.0184,		162352
	DD			RPF+4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48Kp_1P3, OFFSET fft48Kp_2P3
	DD			OFFSET fft48Kp_3P3, OFFSET fft48Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			96, 1
	DD			96, 1, 1, 1, 0
	DD	1006000, 49152,	0.0184,		162352
	DD			4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft48Kp_1PPRO, OFFSET fft48Kp_2PPRO
	DD			OFFSET fft48Kp_3PPRO, OFFSET fft48Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			96, 1
	DD			96, 1, 1, 1, 0
	DD	1345000, 65536,	0.03283,	213504
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64Kp_1P3, OFFSET fft64Kp_2P3
	DD			OFFSET fft64Kp_3P3, OFFSET fft64Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			128, 1
	DD			128, 1, 1, 1, 0
	DD	1345000, 65536,	0.03283,	213504
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft64Kp_1PPRO, OFFSET fft64Kp_2PPRO
	DD			OFFSET fft64Kp_3PPRO, OFFSET fft64Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			128, 1
	DD			128, 1, 1, 1, 0
	DD	1983000, 98304, 0.0445,		290512
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96Kp_1P3, OFFSET fft96Kp_2P3
	DD			OFFSET fft96Kp_3P3, OFFSET fft96Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			48, 1
	DD			48, 1, 1, 1, 0
	DD	1983000, 98304, 0.0445,		290512
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft96Kp_1PPRO, OFFSET fft96Kp_2PPRO
	DD			OFFSET fft96Kp_3PPRO, OFFSET fft96Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			48, 1
	DD			48, 1, 1, 1, 0
	DD	2652000, 131072, 0.0719,	377600
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128Kp_1P3, OFFSET fft128Kp_2P3
	DD			OFFSET fft128Kp_3P3, OFFSET fft128Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			64, 1
	DD			64, 1, 1, 1, 0
	DD	2652000, 131072, 0.0719,	377600
	DD			2		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DD			OFFSET fft128Kp_1PPRO, OFFSET fft128Kp_2PPRO
	DD			OFFSET fft128Kp_3PPRO, OFFSET fft128Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			64, 1
	DD			64, 1, 1, 1, 0
	DD	3924000, 196608, 0.0982,	552496
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			OFFSET fft192Kp_1P3, OFFSET fft192Kp_2P3
	DD			OFFSET fft192Kp_3P3, OFFSET fft192Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			96, 1
	DD			96, 1, 1, 1, 0
	DD	3924000, 196608, 0.0982,	552496
	DD			2		;; Flags, min_l2_cache, clm
	DD			6144		;; scratch area size
	DD			OFFSET fft192Kp_1PPRO, OFFSET fft192Kp_2PPRO
	DD			OFFSET fft192Kp_3PPRO, OFFSET fft192Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			96, 1
	DD			96, 1, 1, 1, 0
	DD	5242000, 262144, 0.155,		726528
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			OFFSET fft256Kp_1P3, OFFSET fft256Kp_2P3
	DD			OFFSET fft256Kp_3P3, OFFSET fft256Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			128, 1
	DD			128, 1, 1, 1, 0
	DD	5242000, 262144, 0.155,		726528
	DD			2		;; Flags, min_l2_cache, clm
	DD			8192+3*64	;; scratch area size
	DD			OFFSET fft256Kp_1PPRO, OFFSET fft256Kp_2PPRO
	DD			OFFSET fft256Kp_3PPRO, OFFSET fft256Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			128, 1
	DD			128, 1, 1, 1, 0
	DD	7733000, 393216, 0.2107,	1076176
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft384Kp_1P3, OFFSET fft384Kp_2P3
	DD			OFFSET fft384Kp_3P3, OFFSET fft384Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			192, 1
	DD			192, 1, 1, 1, 0
	DD	7733000, 393216, 0.2107,	1076176
	DD			2		;; Flags, min_l2_cache, clm
	DD			12288+5*64	;; scratch area size
	DD			OFFSET fft384Kp_1PPRO, OFFSET fft384Kp_2PPRO
	DD			OFFSET fft384Kp_3PPRO, OFFSET fft384Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			192, 1
	DD			192, 1, 1, 1, 0
	DD	10320000, 524288, 0.322,	1424384
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft512Kp_1P3, OFFSET fft512Kp_2P3
	DD			OFFSET fft512Kp_3P3, OFFSET fft512Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			256, 1
	DD			256, 1, 1, 1, 0
	DD	10320000, 524288, 0.322,	1424384
	DD			2		;; Flags, min_l2_cache, clm
	DD			16384+7*64	;; scratch area size
	DD			OFFSET fft512Kp_1PPRO, OFFSET fft512Kp_2PPRO
	DD			OFFSET fft512Kp_3PPRO, OFFSET fft512Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			256, 1
	DD			256, 1, 1, 1, 0
	DD	15260000, 786432, 0.453,	2123824
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft768Kp_1P3, OFFSET fft768Kp_2P3
	DD			OFFSET fft768Kp_3P3, OFFSET fft768Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			384, 1
	DD			384, 1, 1, 1, 0
	DD	15260000, 786432, 0.453,	2123824
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft768Kp_1PPRO, OFFSET fft768Kp_2PPRO
	DD			OFFSET fft768Kp_3PPRO, OFFSET fft768Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			384, 1
	DD			384, 1, 1, 1, 0
	DD	20360000, 1048576, 0.681,	2820096
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft1024Kp_1P3, OFFSET fft1024Kp_2P3
	DD			OFFSET fft1024Kp_3P3, OFFSET fft1024Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			512, 1
	DD			512, 1, 1, 1, 0
	DD	20360000, 1048576, 0.681,	2820096
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft1024Kp_1PPRO, OFFSET fft1024Kp_2PPRO
	DD			OFFSET fft1024Kp_3PPRO, OFFSET fft1024Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			512, 1
	DD			512, 1, 1, 1, 0
	DD	30070000, 1572864, 0.934,	4111312
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+5*64	;; scratch area size
	DD			OFFSET fft1536Kp_1P3, OFFSET fft1536Kp_2P3
	DD			OFFSET fft1536Kp_3P3, OFFSET fft1536Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			192, 1
	DD			192, 1, 1, 1, 0
	DD	30070000, 1572864, 0.934,	4111312
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+5*64	;; scratch area size
	DD			OFFSET fft1536Kp_1PPRO, OFFSET fft1536Kp_2PPRO
	DD			OFFSET fft1536Kp_3PPRO, OFFSET fft1536Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			192, 1
	DD			192, 1, 1, 1, 0
	DD	40110000, 2097152, 1.380,	5442560
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+7*64	;; scratch area size
	DD			OFFSET fft2048Kp_1P3, OFFSET fft2048Kp_2P3
	DD			OFFSET fft2048Kp_3P3, OFFSET fft2048Kp_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, 1
	DD			256, 1, 1, 1, 0
	DD	40110000, 2097152, 1.380,	5442560
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+7*64	;; scratch area size
	DD			OFFSET fft2048Kp_1PPRO, OFFSET fft2048Kp_2PPRO
	DD			OFFSET fft2048Kp_3PPRO, OFFSET fft2048Kp_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, 1
	DD			256, 1, 1, 1, 0
	DD	59360000, 3145728, 1.990,	8108080
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft3072Kp2_1P3, OFFSET fft3072Kp2_2P3
	DD			OFFSET fft3072Kp2_3P3, OFFSET fft3072Kp2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, 1
	DD			384, 1, 1, 1, 0
	DD	59360000, 3145728, 1.990,	8108080
	DD			2		;; Flags, min_l2_cache, clm
	DD			24576+11*64	;; scratch area size
	DD			OFFSET fft3072Kp2_1PPRO, OFFSET fft3072Kp2_2PPRO
	DD			OFFSET fft3072Kp2_3PPRO, OFFSET fft3072Kp2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, 1
	DD			384, 1, 1, 1, 0
	DD	79100000, 4194304, 2.919,	10770432
	DD			RPF+2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft4096Kp2_1P3, OFFSET fft4096Kp2_2P3
	DD			OFFSET fft4096Kp2_3P3, OFFSET fft4096Kp2_4P3
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, 1
	DD			512, 1, 1, 1, 0
	DD	79100000, 4194304, 2.919,	10770432
	DD			2		;; Flags, min_l2_cache, clm
	DD			32768+15*64	;; scratch area size
	DD			OFFSET fft4096Kp2_1PPRO, OFFSET fft4096Kp2_2PPRO
	DD			OFFSET fft4096Kp2_3PPRO, OFFSET fft4096Kp2_4PPRO
	DD			OFFSET prctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, 1
	DD			512, 1, 1, 1, 0
	DD	0
ENDIF

;; Jump tables for the Pentium 4 SSE2 optimized code

CELE_D	EQU	257		;; 256K L2 cache of a Celeron D
WILLI	EQU	256		;; 256K L2 cache of a Willamette P4

xjmptable DD	743,	32,	0.00000111,	896
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft32_1, OFFSET xfft32_2
	DP			OFFSET xfft32_3, OFFSET xfft32_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			4, 4
	DD			1, 1, 1, 1, 1, 0
	DD	1099,	48,	0.00000144,	1408
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft48_1, OFFSET xfft48_2
	DP			OFFSET xfft48_3, OFFSET xfft48_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			6, 6
	DD			2, 1, 1, 1, 1, 0
	DD	1469,	64,	0.00000178,	1920
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft64_1, OFFSET xfft64_2
	DP			OFFSET xfft64_3, OFFSET xfft64_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			8, 8
	DD			3, 1, 1, 1, 1, 0
	DD	1827,	80,	0.00000222,	2176
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft80_1, OFFSET xfft80_2
	DP			OFFSET xfft80_3, OFFSET xfft80_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			10, 8*2048+2
	DD			4, 2, 1, 1, 1, 0
	DD	2179,	96,	0.00000259,	2432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft96_1, OFFSET xfft96_2
	DP			OFFSET xfft96_3, OFFSET xfft96_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			12, 12
	DD			5, 2, 1, 1, 1, 0
	DD	2539,	112,	0.00000311,	2944
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft112_1, OFFSET xfft112_2
	DP			OFFSET xfft112_3, OFFSET xfft112_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			14, (8*2048+4)*2048+2
	DD			6, 3, 1, 1, 1, 0
	DD	2905,	128,	0.00000319,	3328
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft128_1, OFFSET xfft128_2
	DP			OFFSET xfft128_3, OFFSET xfft128_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			16, 16
	DD			7, 3, 1, 1, 1, 0
	DD	3613,	160,	0.00000450,	4736
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft160_1, OFFSET xfft160_2
	DP			OFFSET xfft160_3, OFFSET xfft160_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			20, 16*2048+4
	DD			9, 9, 2, 1, 4, 0
	DD	4311,	192,	0.00000542,	5632
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft192_1, OFFSET xfft192_2
	DP			OFFSET xfft192_3, OFFSET xfft192_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			24, 24
	DD			11, 11, 2, 1, 5, 0
	DD	5029,	224,	0.00000663,	6656
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft224_1, OFFSET xfft224_2
	DP			OFFSET xfft224_3, OFFSET xfft224_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			28, (16*2048+8)*2048+4
	DD			13, 13, 3, 1, 6, 0
	DD	5755,	256,	0.00000691,	7296
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft256_1, OFFSET xfft256_2
	DP			OFFSET xfft256_3, OFFSET xfft256_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			32, 32
	DD			15, 15, 3, 1, 7, 0
	DD	7149,	320,	0.00000928,	8448
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft320_1, OFFSET xfft320_2
	DP			OFFSET xfft320_3, OFFSET xfft320_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			40, 32*2048+8
	DD			19, 9, 2, 4, 1, 0
	DD	8527,	384,	0.0000111,	9984
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft384_1, OFFSET xfft384_2
	DP			OFFSET xfft384_3, OFFSET xfft384_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			48, 48
	DD			23, 11, 2, 5, 1, 0
	DD	9933,	448,	0.0000133,	11648
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft448_1, OFFSET xfft448_2
	DP			OFFSET xfft448_3, OFFSET xfft448_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			56, (32*2048+16)*2048+8
	DD			27, 13, 3, 6, 1, 0
	DD	11359,	512,	0.0000143,	13056
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft512_1, OFFSET xfft512_2
	DP			OFFSET xfft512_3, OFFSET xfft512_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			64, 64
	DD			31, 15, 3, 7, 1, 0
	DD	14119,	640,	0.0000215,	17408
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft640_1, OFFSET xfft640_2
	DP			OFFSET xfft640_3, OFFSET xfft640_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			80, 64*2048+16
	DD			39, 19, 9, 9, 4*256+2, 0
	DD	16839,	768,	0.0000260,	20736
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft768_1, OFFSET xfft768_2
	DP			OFFSET xfft768_3, OFFSET xfft768_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			96, 96
	DD			47, 23, 11, 11, 5*256+2, 0
	DD	19639,	896,	0.0000321,	24448
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft896_1, OFFSET xfft896_2
	DP			OFFSET xfft896_3, OFFSET xfft896_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			112, (64*2048+32)*2048+16
	DD			55, 27, 13, 13, 6*256+3, 0
	DD	22477,	1024,	0.0000349,	26112
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft1024_1, OFFSET xfft1024_2
	DP			OFFSET xfft1024_3, OFFSET xfft1024_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			128, 128
	DD			63, 31, 15, 15, 7*256+3, 0
	DD	27899,	1280,	0.0000494,	33664
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft1280_1, OFFSET xfft1280_2
	DP			OFFSET xfft1280_3, OFFSET xfft1280_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			160, 128*2048+32
	DD			79, 39, 9, 4*256+19, 2, 0
	DD	33289,	1536,	0.0000601,	40320
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft1536_1, OFFSET xfft1536_2
	DP			OFFSET xfft1536_3, OFFSET xfft1536_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			192, 192
	DD			95, 47, 11, 5*256+23, 2, 0
	DD	38799,	1792,	0.0000719,	47232
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft1792_1, OFFSET xfft1792_2
	DP			OFFSET xfft1792_3, OFFSET xfft1792_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			224, (128*2048+64)*2048+32
	DD			111, 55, 13, 6*256+27, 3, 0
	DD	44339,	2048,	0.0000773,	52224
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft2048_1, OFFSET xfft2048_2
	DP			OFFSET xfft2048_3, OFFSET xfft2048_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			256, 256
	DD			127, 63, 15, 7*256+31, 3, 0
	DD	55099,	2560,	0.000111,	68096
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft2560_1, OFFSET xfft2560_2
	DP			OFFSET xfft2560_3, OFFSET xfft2560_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			320, 256*2048+64
	DD			159, 79, 9*256+19, 9*256+39, 4*256+2, 0
	DD	65729,	3072,	0.000131,	81792
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft3072_1, OFFSET xfft3072_2
	DP			OFFSET xfft3072_3, OFFSET xfft3072_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			384, 384
	DD			191, 95, 11*256+23, 11*256+47, 5*256+2, 0
	DD	76559,	3584,	0.000165,	95488
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft3584_1, OFFSET xfft3584_2
	DP			OFFSET xfft3584_3, OFFSET xfft3584_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			448, (256*2048+128)*2048+64
	DD			223, 111, 13*256+27, 13*256+55, 6*256+3, 0
	DD	87549,	4096,	0.000163,	104448
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft4096_1, OFFSET xfft4096_2
	DP			OFFSET xfft4096_3, OFFSET xfft4096_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			512, 512
	DD			255, 127, 15*256+31, 15*256+63, 7*256+3, 0
	DD	108800,	5120,	0.000215,	135296
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft5120_1, OFFSET xfft5120_2
	DP			OFFSET xfft5120_3, OFFSET xfft5120_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			640, 512*2048+128
	DD			319, 159, 9*256+39, 19*256+79, 4*256+2, 0
	DD	129900,	6144,	0.000276,	162432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft6144_1, OFFSET xfft6144_2
	DP			OFFSET xfft6144_3, OFFSET xfft6144_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			768, 768
	DD			383, 191, 11*256+47, 23*256+95, 5*256+2, 0
	DD	151300,	7168,	0.000374,	189568
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft7168_1, OFFSET xfft7168_2
	DP			OFFSET xfft7168_3, OFFSET xfft7168_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			896, (512*2048+256)*2048+128
	DD			447, 223, 13*256+55, 27*256+111, 6*256+3, 0
	DD	172700,	8192,	0.000398,	208896
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft8192_1, OFFSET xfft8192_2
	DP			OFFSET xfft8192_3, OFFSET xfft8192_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1024, 1024
	DD			511, 255, 15*256+63, 31*256+127, 7*256+3, 0
	DD	214400,	10240,	0.000470,	59904
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft10K_1AMD, OFFSET xfft10K_2AMD
	DP			OFFSET xfft10K_3AMD, OFFSET xfft10K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			10, GAP2_8_4+1
	DD			9, 1, 8*2048+2, 1, 1, 0
	DD	214400,	10240,	0.000470,	59904
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft10K_1, OFFSET xfft10K_2
	DP			OFFSET xfft10K_3, OFFSET xfft10K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			10, GAP2_8_4+1
	DD			9, 1, 8*2048+2, 1, 1, 0
	DD	255300,	12288,	0.000590,	69632
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft12K_1AMD, OFFSET xfft12K_2AMD
	DP			OFFSET xfft12K_3AMD, OFFSET xfft12K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			12, GAP2_8_4+1
	DD			11, 1, 12, 1, 1, 0
	DD	255300,	12288,	0.000590,	69632
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft12K_1, OFFSET xfft12K_2
	DP			OFFSET xfft12K_3, OFFSET xfft12K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			12, GAP2_8_4+1
	DD			11, 1, 12, 1, 1, 0
	DD	297300,	14336,	0.000716,	79488
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft14K_1AMD, OFFSET xfft14K_2AMD
	DP			OFFSET xfft14K_3AMD, OFFSET xfft14K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			14, GAP2_8_4+1
	DD			13, 1, (8*2048+4)*2048+2, 1, 1, 0
	DD	297300,	14336,	0.000716,	79488
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft14K_1, OFFSET xfft14K_2
	DP			OFFSET xfft14K_3, OFFSET xfft14K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			14, GAP2_8_4+1
	DD			13, 1, (8*2048+4)*2048+2, 1, 1, 0
	DD	340400,	16384,	0.000787,	89088
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft16K_1AMD, OFFSET xfft16K_2AMD
	DP			OFFSET xfft16K_3AMD, OFFSET xfft16K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			16, GAP2_8_4+1
	DD			15, 1, 16, 1, 1, 0
	DD	340400,	16384,	0.000787,	89088
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft16K_1, OFFSET xfft16K_2
	DP			OFFSET xfft16K_3, OFFSET xfft16K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			16, GAP2_8_4+1
	DD			15, 1, 16, 1, 1, 0
	DD	423300,	20480,	0.00103,	107904
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft20K_1AMD, OFFSET xfft20K_2AMD
	DP			OFFSET xfft20K_3AMD, OFFSET xfft20K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			20, GAP2_8_4+1
	DD			19, 1, 16*2048+4, 1, 1, 0
	DD	423300,	20480,	0.00103,	107904
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft20K_1, OFFSET xfft20K_2
	DP			OFFSET xfft20K_3, OFFSET xfft20K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			20, GAP2_8_4+1
	DD			19, 1, 16*2048+4, 1, 1, 0
	DD	504600,	24576,	0.00132,	127232
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft24K_1AMD, OFFSET xfft24K_2AMD
	DP			OFFSET xfft24K_3AMD, OFFSET xfft24K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			24, GAP2_8_4+1
	DD			23, 1, 24, 1, 1, 0
	DD	504600,	24576,	0.00132,	127232
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft24K_1, OFFSET xfft24K_2
	DP			OFFSET xfft24K_3, OFFSET xfft24K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			24, GAP2_8_4+1
	DD			23, 1, 24, 1, 1, 0
	DD	587500,	28672,	0.00156,	146688
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft28K_1AMD, OFFSET xfft28K_2AMD
	DP			OFFSET xfft28K_3AMD, OFFSET xfft28K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			28, GAP2_8_4+1
	DD			27, 1, (16*2048+8)*2048+4, 1, 1, 0
	DD	587500,	28672,	0.00156,	146688
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft28K_1, OFFSET xfft28K_2
	DP			OFFSET xfft28K_3, OFFSET xfft28K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			28, GAP2_8_4+1
	DD			27, 1, (16*2048+8)*2048+4, 1, 1, 0
	DD	671400,	32768,	0.00175,	165888
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft32K_1AMD, OFFSET xfft32K_2AMD
	DP			OFFSET xfft32K_3AMD, OFFSET xfft32K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			32, GAP2_8_4+1
	DD			31, 1, 32, 1, 1, 0
	DD	671400,	32768,	0.00175,	165888
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft32K_1, OFFSET xfft32K_2
	DP			OFFSET xfft32K_3, OFFSET xfft32K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			32, GAP2_8_4+1
	DD			31, 1, 32, 1, 1, 0
	DD	835200,	40960,	0.00225,	189696
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft40K_1AMD, OFFSET xfft40K_2AMD
	DP			OFFSET xfft40K_3AMD, OFFSET xfft40K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			5, GAP2_11_4+1
	DD			4, 1, 4*2048+1, 1, 1, 0
	DD	835200,	40960,	0.00225,	189696
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft40K_1, OFFSET xfft40K_2
	DP			OFFSET xfft40K_3, OFFSET xfft40K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			5, GAP2_11_4+1
	DD			4, 1, 4*2048+1, 1, 1, 0
allfft	DD	835200,	40960,	0.00225,	156672
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft40K_1AMD, OFFSET xfft40K_2AMD
allfft	DP			OFFSET xfft40K_3AMD, OFFSET xfft40K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			10, GAP2_10_4+1
allfft	DD			9, 1, 8*2048+2, 1, 1, 0
allfft	DD	835200,	40960,	0.00225,	156672
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft40K_1, OFFSET xfft40K_2
allfft	DP			OFFSET xfft40K_3, OFFSET xfft40K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			10, GAP2_10_4+1
allfft	DD			9, 1, 8*2048+2, 1, 1, 0
	DD	995500, 49152,	0.00279,	206208
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft48K_1AMD, OFFSET xfft48K_2AMD
	DP			OFFSET xfft48K_3AMD, OFFSET xfft48K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			6, GAP2_11_4+1
	DD			5, 1, 6, 1, 1, 0
	DD	995500, 49152,	0.00279,	206208
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft48K_1, OFFSET xfft48K_2
	DP			OFFSET xfft48K_3, OFFSET xfft48K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			6, GAP2_11_4+1
	DD			5, 1, 6, 1, 1, 0
allfft	DD	995500, 49152,	0.00279,	177664
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft48K_1AMD, OFFSET xfft48K_2AMD
allfft	DP			OFFSET xfft48K_3AMD, OFFSET xfft48K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			12, GAP2_10_4+1
allfft	DD			11, 1, 12, 1, 1, 0
allfft	DD	995500, 49152,	0.00279,	177664
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft48K_1, OFFSET xfft48K_2
allfft	DP			OFFSET xfft48K_3, OFFSET xfft48K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			12, GAP2_10_4+1
allfft	DD			11, 1, 12, 1, 1, 0
	DD	1158000, 57344,	0.00327,	222976
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft56K_1AMD, OFFSET xfft56K_2AMD
	DP			OFFSET xfft56K_3AMD, OFFSET xfft56K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			7, GAP2_11_4+1
	DD			6, 1, (4*2048+2)*2048+1, 1, 1, 0
	DD	1158000, 57344,	0.00327,	222976
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft56K_1, OFFSET xfft56K_2
	DP			OFFSET xfft56K_3, OFFSET xfft56K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			7, GAP2_11_4+1
	DD			6, 1, (4*2048+2)*2048+1, 1, 1, 0
allfft	DD	1158000, 57344,	0.00327,	198784
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft56K_1AMD, OFFSET xfft56K_2AMD
allfft	DP			OFFSET xfft56K_3AMD, OFFSET xfft56K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			14, GAP2_10_4+1
allfft	DD			13, 1, (8*2048+4)*2048+2, 1, 1, 0
allfft	DD	1158000, 57344,	0.00327,	198784
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft56K_1, OFFSET xfft56K_2
allfft	DP			OFFSET xfft56K_3, OFFSET xfft56K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			14, GAP2_10_4+1
allfft	DD			13, 1, (8*2048+4)*2048+2, 1, 1, 0
	DD	1325000, 65536,	0.00367,	239488
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft64K_1AMD, OFFSET xfft64K_2AMD
	DP			OFFSET xfft64K_3AMD, OFFSET xfft64K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			8, GAP2_11_4+1
	DD			7, 1, 8, 1, 1, 0
	DD	1325000, 65536,	0.00367,	239488
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft64K_1, OFFSET xfft64K_2
	DP			OFFSET xfft64K_3, OFFSET xfft64K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			8, GAP2_11_4+1
	DD			7, 1, 8, 1, 1, 0
allfft	DD	1325000, 65536,	0.00367,	219648
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft64K_1AMD, OFFSET xfft64K_2AMD
allfft	DP			OFFSET xfft64K_3AMD, OFFSET xfft64K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			16, GAP2_10_4+1
allfft	DD			15, 1, 16, 1, 1, 0
allfft	DD	1325000, 65536,	0.00367,	219648
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft64K_1, OFFSET xfft64K_2
allfft	DP			OFFSET xfft64K_3, OFFSET xfft64K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			16, GAP2_10_4+1
allfft	DD			15, 1, 16, 1, 1, 0
	DD	1648000, 81920, 0.00474,	273408
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft80K_1AMD, OFFSET xfft80K_2AMD
	DP			OFFSET xfft80K_3AMD, OFFSET xfft80K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			10, GAP2_11_4+1
	DD			9, 1, 8*2048+2, 1, 1, 0
	DD	1648000, 81920, 0.00474,	273408
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft80K_1, OFFSET xfft80K_2
	DP			OFFSET xfft80K_3, OFFSET xfft80K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			10, GAP2_11_4+1
	DD			9, 1, 8*2048+2, 1, 1, 0
allfft	DD	1648000, 81920, 0.00474,	260992
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft80K_1AMD, OFFSET xfft80K_2AMD
allfft	DP			OFFSET xfft80K_3AMD, OFFSET xfft80K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			20, GAP2_10_4+1
allfft	DD			19, 1, 10*2048+4, 1, 1, 0
allfft	DD	1648000, 81920, 0.00474,	260992
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft80K_1, OFFSET xfft80K_2
allfft	DP			OFFSET xfft80K_3, OFFSET xfft80K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			20, GAP2_10_4+1
allfft	DD			19, 1, 16*2048+4, 1, 1, 0
	DD	1966000, 98304, 0.00584,	306688
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft96K_1AMD, OFFSET xfft96K_2AMD
	DP			OFFSET xfft96K_3AMD, OFFSET xfft96K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			12, GAP2_11_4+1
	DD			11, 1, 12, 1, 1, 0
	DD	1966000, 98304, 0.00584,	306688
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft96K_1, OFFSET xfft96K_2
	DP			OFFSET xfft96K_3, OFFSET xfft96K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			12, GAP2_11_4+1
	DD			11, 1, 12, 1, 1, 0
allfft	DD	1966000, 98304, 0.00584,	302848
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft96K_1AMD, OFFSET xfft96K_2AMD
allfft	DP			OFFSET xfft96K_3AMD, OFFSET xfft96K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			24, GAP2_10_4+1
allfft	DD			23, 1, 24, 1, 1, 0
allfft	DD	1966000, 98304, 0.00584,	302848
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft96K_1, OFFSET xfft96K_2
allfft	DP			OFFSET xfft96K_3, OFFSET xfft96K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			24, GAP2_10_4+1
allfft	DD			23, 1, 24, 1, 1, 0
	DD	2287000, 114688, 0.00693,	340096
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft112K_1AMD, OFFSET xfft112K_2AMD
	DP			OFFSET xfft112K_3AMD, OFFSET xfft112K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			14, GAP2_11_4+1
	DD			13, 1, (8*2048+4)*2048+2, 1, 1, 0
	DD	2287000, 114688, 0.00693,	340096
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft112K_1, OFFSET xfft112K_2
	DP			OFFSET xfft112K_3, OFFSET xfft112K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			14, GAP2_11_4+1
	DD			13, 1, (8*2048+4)*2048+2, 1, 1, 0
allfft	DD	2287000, 114688, 0.00693,	344832
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft112K_1AMD, OFFSET xfft112K_2AMD
allfft	DP			OFFSET xfft112K_3AMD, OFFSET xfft112K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			28, GAP2_10_4+1
allfft	DD			27, 1, (16*2048+8)*2048+4, 1, 1, 0
allfft	DD	2287000, 114688, 0.00693,	344832
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft112K_1, OFFSET xfft112K_2
allfft	DP			OFFSET xfft112K_3, OFFSET xfft112K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			28, GAP2_10_4+1
allfft	DD			27, 1, (16*2048+8)*2048+4, 1, 1, 0
	DD	2614000, 131072, 0.00779,	373248
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft128K_1AMD, OFFSET xfft128K_2AMD
	DP			OFFSET xfft128K_3AMD, OFFSET xfft128K_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			16, GAP2_11_4+1
	DD			15, 1, 16, 1, 1, 0
	DD	2614000, 131072, 0.00779,	373248
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft128K_1, OFFSET xfft128K_2
	DP			OFFSET xfft128K_3, OFFSET xfft128K_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			16, GAP2_11_4+1
	DD			15, 1, 16, 1, 1, 0
allfft	DD	2614000, 131072, 0.00779,	386560
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft128K_1AMD, OFFSET xfft128K_2AMD
allfft	DP			OFFSET xfft128K_3AMD, OFFSET xfft128K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			32, GAP2_10_4+1
allfft	DD			31, 1, 32, 1, 1, 0
allfft	DD	2614000, 131072, 0.00779,	386560
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft128K_1, OFFSET xfft128K_2
allfft	DP			OFFSET xfft128K_3, OFFSET xfft128K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			32, GAP2_10_4+1
allfft	DD			31, 1, 32, 1, 1, 0
allfft	DD	3251000, 163840, 0.00914,	439168
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft160K_1AMD, OFFSET xfft160K_2AMD
allfft	DP			OFFSET xfft160K_3AMD, OFFSET xfft160K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			20, GAP2_11_4+1
allfft	DD			19, 1, 16*2048+4, 1, 1, 0
allfft	DD	3251000, 163840, 0.00914,	471424
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET xfft160K410_1AMD
allfft	DP			OFFSET xfft160K410_2AMD
allfft	DP			OFFSET xfft160K410_3AMD
allfft	DP			OFFSET xfft160K410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			40, GAP2_10_4+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
	DD	3251000, 163840, 0.00914,	471424
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			10240+4*128	;; scratch area size
	DP			OFFSET xfft160K210_1AMD
	DP			OFFSET xfft160K210_2AMD
	DP			OFFSET xfft160K210_3AMD
	DP			OFFSET xfft160K210_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			40, GAP2_10_2+1
	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	3251000, 163840, 0.00914,	471424
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			5120+4*128	;; scratch area size
allfft	DP			OFFSET xfft160K110_1AMD
allfft	DP			OFFSET xfft160K110_2AMD
allfft	DP			OFFSET xfft160K110_3AMD
allfft	DP			OFFSET xfft160K110_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			40, GAP2_10_1+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	3251000, 163840, 0.00914,	439168
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft160K_1, OFFSET xfft160K_2
allfft	DP			OFFSET xfft160K_3, OFFSET xfft160K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			20, GAP2_11_4+1
allfft	DD			19, 1, 16*2048+4, 1, 1, 0
allfft	DD	3251000, 163840, 0.00914,	471424
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET xfft160K410_1, OFFSET xfft160K410_2
allfft	DP			OFFSET xfft160K410_3, OFFSET xfft160K410_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			40, GAP2_10_4+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
	DD	3251000, 163840, 0.00914,	471424
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			10240+4*128	;; scratch area size
	DP			OFFSET xfft160K210_1, OFFSET xfft160K210_2
	DP			OFFSET xfft160K210_3, OFFSET xfft160K210_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			40, GAP2_10_2+1
	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	3251000, 163840, 0.00914,	471424
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			5120+4*128	;; scratch area size
allfft	DP			OFFSET xfft160K110_1, OFFSET xfft160K110_2
allfft	DP			OFFSET xfft160K110_3, OFFSET xfft160K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			40, GAP2_10_1+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	3875000, 196608, 0.0114,	505600
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft192K_1AMD, OFFSET xfft192K_2AMD
allfft	DP			OFFSET xfft192K_3AMD, OFFSET xfft192K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			24, GAP2_11_4+1
allfft	DD			23, 1, 24, 1, 1, 0
allfft	DD	3875000, 196608, 0.0114,	555392
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft192K410_1AMD
allfft	DP			OFFSET xfft192K410_2AMD
allfft	DP			OFFSET xfft192K410_3AMD
allfft	DP			OFFSET xfft192K410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			48, GAP2_10_4+1
allfft	DD			47, 1, 48, 1, 1, 0
	DD	3875000, 196608, 0.0114,	555392
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET xfft192K210_1AMD
	DP			OFFSET xfft192K210_2AMD
	DP			OFFSET xfft192K210_3AMD
	DP			OFFSET xfft192K210_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			48, GAP2_10_2+1
	DD			47, 1, 48, 1, 1, 0
allfft	DD	3875000, 196608, 0.0114,	555392
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			6144+5*128	;; scratch area size
allfft	DP			OFFSET xfft192K110_1AMD
allfft	DP			OFFSET xfft192K110_2AMD
allfft	DP			OFFSET xfft192K110_3AMD
allfft	DP			OFFSET xfft192K110_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			48, GAP2_10_1+1
allfft	DD			47, 1, 48, 1, 1, 0
allfft	DD	3875000, 196608, 0.0114,	505600
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft192K_1, OFFSET xfft192K_2
allfft	DP			OFFSET xfft192K_3, OFFSET xfft192K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			24, GAP2_11_4+1
allfft	DD			23, 1, 24, 1, 1, 0
	DD	3875000, 196608, 0.0114,	555392
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			24576+5*128	;; scratch area size
	DP			OFFSET xfft192K410_1, OFFSET xfft192K410_2
	DP			OFFSET xfft192K410_3, OFFSET xfft192K410_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			48, GAP2_10_4+1
	DD			47, 1, 48, 1, 1, 0
allfft	DD	3875000, 196608, 0.0114,	555392
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			12288+5*128	;; scratch area size
allfft	DP			OFFSET xfft192K210_1, OFFSET xfft192K210_2
allfft	DP			OFFSET xfft192K210_3, OFFSET xfft192K210_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			48, GAP2_10_2+1
allfft	DD			47, 1, 48, 1, 1, 0
	DD	3875000, 196608, 0.0114,	555392
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			6144+5*128	;; scratch area size
	DP			OFFSET xfft192K110_1, OFFSET xfft192K110_2
	DP			OFFSET xfft192K110_3, OFFSET xfft192K110_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			48, GAP2_10_1+1
	DD			47, 1, 48, 1, 1, 0
allfft	DD	4512000, 229376, 0.0134,	572160
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft224K_1AMD, OFFSET xfft224K_2AMD
allfft	DP			OFFSET xfft224K_3AMD, OFFSET xfft224K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			28, GAP2_11_4+1
allfft	DD			27, 1, (16*2048+8)*2048+4, 1, 1, 0
allfft	DD	4512000, 229376, 0.0134,	639616
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET xfft224K410_1AMD
allfft	DP			OFFSET xfft224K410_2AMD
allfft	DP			OFFSET xfft224K410_3AMD
allfft	DP			OFFSET xfft224K410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			56, GAP2_10_4+1
allfft	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
	DD	4512000, 229376, 0.0134,	639616
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			14336+6*128	;; scratch area size
	DP			OFFSET xfft224K210_1AMD
	DP			OFFSET xfft224K210_2AMD
	DP			OFFSET xfft224K210_3AMD
	DP			OFFSET xfft224K210_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			56, GAP2_10_2+1
	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
allfft	DD	4512000, 229376, 0.0134,	639616
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			7168+6*128	;; scratch area size
allfft	DP			OFFSET xfft224K110_1AMD
allfft	DP			OFFSET xfft224K110_2AMD
allfft	DP			OFFSET xfft224K110_3AMD
allfft	DP			OFFSET xfft224K110_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			56, GAP2_10_1+1
allfft	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
allfft	DD	4512000, 229376, 0.0134,	572160
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft224K_1, OFFSET xfft224K_2
allfft	DP			OFFSET xfft224K_3, OFFSET xfft224K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			28, GAP2_11_4+1
allfft	DD			27, 1, (16*2048+8)*2048+4, 1, 1, 0
	DD	4512000, 229376, 0.0134,	639616
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			28672+6*128	;; scratch area size
	DP			OFFSET xfft224K410_1, OFFSET xfft224K410_2
	DP			OFFSET xfft224K410_3, OFFSET xfft224K410_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			56, GAP2_10_4+1
	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
	DD	4512000, 229376, 0.0134,	639616
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			14336+6*128	;; scratch area size
	DP			OFFSET xfft224K210_1, OFFSET xfft224K210_2
	DP			OFFSET xfft224K210_3, OFFSET xfft224K210_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			56, GAP2_10_2+1
	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
allfft	DD	4512000, 229376, 0.0134,	639616
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			7168+6*128	;; scratch area size
allfft	DP			OFFSET xfft224K110_1, OFFSET xfft224K110_2
allfft	DP			OFFSET xfft224K110_3, OFFSET xfft224K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			56, GAP2_10_1+1
allfft	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
allfft	DD	5158000, 262144, 0.0150,	638464
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft256K_1AMD, OFFSET xfft256K_2AMD
allfft	DP			OFFSET xfft256K_3AMD, OFFSET xfft256K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			32, GAP2_11_4+1
allfft	DD			31, 1, 32, 1, 1, 0
allfft	DD	5158000, 262144, 0.0150,	721920
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET xfft256K410_1AMD
allfft	DP			OFFSET xfft256K410_2AMD
allfft	DP			OFFSET xfft256K410_3AMD
allfft	DP			OFFSET xfft256K410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			64, GAP2_10_4+1
allfft	DD			63, 1, 64, 1, 1, 0
	DD	5158000, 262144, 0.0150,	721920
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			16384+7*128	;; scratch area size
	DP			OFFSET xfft256K210_1AMD
	DP			OFFSET xfft256K210_2AMD
	DP			OFFSET xfft256K210_3AMD
	DP			OFFSET xfft256K210_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			64, GAP2_10_2+1
	DD			63, 1, 64, 1, 1, 0
allfft	DD	5158000, 262144, 0.0150,	721920
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			8192+7*128	;; scratch area size
allfft	DP			OFFSET xfft256K110_1AMD
allfft	DP			OFFSET xfft256K110_2AMD
allfft	DP			OFFSET xfft256K110_3AMD
allfft	DP			OFFSET xfft256K110_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			64, GAP2_10_1+1
allfft	DD			63, 1, 64, 1, 1, 0
allfft	DD	5158000, 262144, 0.0150,	638464
allfft	DD			4096*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft256K_1, OFFSET xfft256K_2
allfft	DP			OFFSET xfft256K_3, OFFSET xfft256K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			32, GAP2_11_4+1
allfft	DD			31, 1, 32, 1, 1, 0
	DD	5158000, 262144, 0.0150,	721920
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			32768+7*128	;; scratch area size
	DP			OFFSET xfft256K410_1, OFFSET xfft256K410_2
	DP			OFFSET xfft256K410_3, OFFSET xfft256K410_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			64, GAP2_10_4+1
	DD			63, 1, 64, 1, 1, 0
	DD	5158000, 262144, 0.0150,	721920
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			16384+7*128	;; scratch area size
	DP			OFFSET xfft256K210_1, OFFSET xfft256K210_2
	DP			OFFSET xfft256K210_3, OFFSET xfft256K210_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			64, GAP2_10_2+1
	DD			63, 1, 64, 1, 1, 0
allfft	DD	5158000, 262144, 0.0150,	721920
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			8192+7*128	;; scratch area size
allfft	DP			OFFSET xfft256K110_1, OFFSET xfft256K110_2
allfft	DP			OFFSET xfft256K110_3, OFFSET xfft256K110_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			64, GAP2_10_1+1
allfft	DD			63, 1, 64, 1, 1, 0
allfft	DD	6421000, 327680, 0.0192,	772480
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET xfft320K411_1AMD
allfft	DP			OFFSET xfft320K411_2AMD
allfft	DP			OFFSET xfft320K411_3AMD
allfft	DP			OFFSET xfft320K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			40, GAP2_11_4+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	6421000, 327680, 0.0192,	772480
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft320K_1AMD, OFFSET xfft320K_2AMD
allfft	DP			OFFSET xfft320K_3AMD, OFFSET xfft320K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			40, GAP2_11_4+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	6421000, 327680, 0.0192,	890624
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			40960+9*128	;; scratch area size
allfft	DP			OFFSET xfft320K410_1AMD
allfft	DP			OFFSET xfft320K410_2AMD
allfft	DP			OFFSET xfft320K410_3AMD
allfft	DP			OFFSET xfft320K410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			80, GAP2_10_4+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	6421000, 327680, 0.0192,	890624
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET xfft320K210_1AMD
allfft	DP			OFFSET xfft320K210_2AMD
allfft	DP			OFFSET xfft320K210_3AMD
allfft	DP			OFFSET xfft320K210_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			80, GAP2_10_2+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
	DD	6421000, 327680, 0.0192,	890624
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			10240+9*128	;; scratch area size
	DP			OFFSET xfft320K110_1AMD
	DP			OFFSET xfft320K110_2AMD
	DP			OFFSET xfft320K110_3AMD
	DP			OFFSET xfft320K110_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			80, GAP2_10_1+1
	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	6421000, 327680, 0.0192,	772480
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET xfft320K411_1, OFFSET xfft320K411_2
allfft	DP			OFFSET xfft320K411_3, OFFSET xfft320K411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			40, GAP2_11_4+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	6421000, 327680, 0.0192,	772480
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft320K_1, OFFSET xfft320K_2
allfft	DP			OFFSET xfft320K_3, OFFSET xfft320K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			40, GAP2_11_4+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
	DD	6421000, 327680, 0.0192,	890624
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			40960+9*128	;; scratch area size
	DP			OFFSET xfft320K410_1, OFFSET xfft320K410_2
	DP			OFFSET xfft320K410_3, OFFSET xfft320K410_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			80, GAP2_10_4+1
	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	6421000, 327680, 0.0192,	890624
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET xfft320K210_1, OFFSET xfft320K210_2
allfft	DP			OFFSET xfft320K210_3, OFFSET xfft320K210_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			80, GAP2_10_2+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
	DD	6421000, 327680, 0.0192,	890624
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			10240+9*128	;; scratch area size
	DP			OFFSET xfft320K110_1, OFFSET xfft320K110_2
	DP			OFFSET xfft320K110_3, OFFSET xfft320K110_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			80, GAP2_10_1+1
	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	7651000, 393216, 0.0238,	905600
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft384K411_1AMD
allfft	DP			OFFSET xfft384K411_2AMD
allfft	DP			OFFSET xfft384K411_3AMD
allfft	DP			OFFSET xfft384K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			48, GAP2_11_4+1
allfft	DD			47, 1, 48, 1, 1, 0
allfft	DD	7651000, 393216, 0.0238,	905600
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft384K_1AMD, OFFSET xfft384K_2AMD
allfft	DP			OFFSET xfft384K_3AMD, OFFSET xfft384K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			48, GAP2_11_4+1
allfft	DD			47, 1, 48, 1, 1, 0
allfft	DD	7651000, 393216, 0.0238,	1058432
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET xfft384K410_1AMD
allfft	DP			OFFSET xfft384K410_2AMD
allfft	DP			OFFSET xfft384K410_3AMD
allfft	DP			OFFSET xfft384K410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			96, GAP2_10_4+1
allfft	DD			95, 1, 96, 1, 1, 0
allfft	DD	7651000, 393216, 0.0238,	1058432
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET xfft384K210_1AMD
allfft	DP			OFFSET xfft384K210_2AMD
allfft	DP			OFFSET xfft384K210_3AMD
allfft	DP			OFFSET xfft384K210_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			96, GAP2_10_2+1
allfft	DD			95, 1, 96, 1, 1, 0
	DD	7651000, 393216, 0.0238,	1058432
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET xfft384K110_1AMD
	DP			OFFSET xfft384K110_2AMD
	DP			OFFSET xfft384K110_3AMD
	DP			OFFSET xfft384K110_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			96, GAP2_10_1+1
	DD			95, 1, 96, 1, 1, 0
allfft	DD	7651000, 393216, 0.0238,	905600
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft384K411_1, OFFSET xfft384K411_2
allfft	DP			OFFSET xfft384K411_3, OFFSET xfft384K411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			48, GAP2_11_4+1
allfft	DD			47, 1, 48, 1, 1, 0
allfft	DD	7651000, 393216, 0.0238,	905600
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft384K_1, OFFSET xfft384K_2
allfft	DP			OFFSET xfft384K_3, OFFSET xfft384K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			48, GAP2_11_4+1
allfft	DD			47, 1, 48, 1, 1, 0
	DD	7651000, 393216, 0.0238,	1058432
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			49152+11*128	;; scratch area size
	DP			OFFSET xfft384K410_1, OFFSET xfft384K410_2
	DP			OFFSET xfft384K410_3, OFFSET xfft384K410_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			96, GAP2_10_4+1
	DD			95, 1, 96, 1, 1, 0
	DD	7651000, 393216, 0.0238,	1058432
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			24576+11*128	;; scratch area size
	DP			OFFSET xfft384K210_1, OFFSET xfft384K210_2
	DP			OFFSET xfft384K210_3, OFFSET xfft384K210_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			96, GAP2_10_2+1
	DD			95, 1, 96, 1, 1, 0
allfft	DD	7651000, 393216, 0.0238,	1058432
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET xfft384K110_1, OFFSET xfft384K110_2
allfft	DP			OFFSET xfft384K110_3, OFFSET xfft384K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			96, GAP2_10_1+1
allfft	DD			95, 1, 96, 1, 1, 0
allfft	DD	8908000, 458752, 0.0283,	1038976
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET xfft448K411_1AMD
allfft	DP			OFFSET xfft448K411_2AMD
allfft	DP			OFFSET xfft448K411_3AMD
allfft	DP			OFFSET xfft448K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			56, GAP2_11_4+1
allfft	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
allfft	DD	8908000, 458752, 0.0283,	1038976
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft448K_1AMD, OFFSET xfft448K_2AMD
allfft	DP			OFFSET xfft448K_3AMD, OFFSET xfft448K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			56, GAP2_11_4+1
allfft	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
allfft	DD	8908000, 458752, 0.0283,	1226496
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			57344+13*128	;; scratch area size
allfft	DP			OFFSET xfft448K410_1AMD
allfft	DP			OFFSET xfft448K410_2AMD
allfft	DP			OFFSET xfft448K410_3AMD
allfft	DP			OFFSET xfft448K410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			112, GAP2_10_4+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	8908000, 458752, 0.0283,	1226496
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET xfft448K210_1AMD
allfft	DP			OFFSET xfft448K210_2AMD
allfft	DP			OFFSET xfft448K210_3AMD
allfft	DP			OFFSET xfft448K210_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			112, GAP2_10_2+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
	DD	8908000, 458752, 0.0283,	1226496
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			14336+13*128	;; scratch area size
	DP			OFFSET xfft448K110_1AMD
	DP			OFFSET xfft448K110_2AMD
	DP			OFFSET xfft448K110_3AMD
	DP			OFFSET xfft448K110_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			112, GAP2_10_1+1
	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	8908000, 458752, 0.0283,	1038976
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET xfft448K411_1, OFFSET xfft448K411_2
allfft	DP			OFFSET xfft448K411_3, OFFSET xfft448K411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			56, GAP2_11_4+1
allfft	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
allfft	DD	8908000, 458752, 0.0283,	1038976
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft448K_1, OFFSET xfft448K_2
allfft	DP			OFFSET xfft448K_3, OFFSET xfft448K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			56, GAP2_11_4+1
allfft	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
	DD	8908000, 458752, 0.0283,	1226496
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			57344+13*128	;; scratch area size
	DP			OFFSET xfft448K410_1, OFFSET xfft448K410_2
	DP			OFFSET xfft448K410_3, OFFSET xfft448K410_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			112, GAP2_10_4+1
	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
	DD	8908000, 458752, 0.0283,	1226496
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			28672+13*128	;; scratch area size
	DP			OFFSET xfft448K210_1, OFFSET xfft448K210_2
	DP			OFFSET xfft448K210_3, OFFSET xfft448K210_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			112, GAP2_10_2+1
	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	8908000, 458752, 0.0283,	1226496
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET xfft448K110_1, OFFSET xfft448K110_2
allfft	DP			OFFSET xfft448K110_3, OFFSET xfft448K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			112, GAP2_10_1+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET xfft512K411_1AMD
allfft	DP			OFFSET xfft512K411_2AMD
allfft	DP			OFFSET xfft512K411_3AMD
allfft	DP			OFFSET xfft512K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			64, GAP2_11_4+1
allfft	DD			63, 1, 64, 1, 1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft512K_1AMD, OFFSET xfft512K_2AMD
allfft	DP			OFFSET xfft512K_3AMD, OFFSET xfft512K_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			64, GAP2_11_4+1
allfft	DD			63, 1, 64, 1, 1, 0
allfft	DD	10180000, 524288, 0.0319,	1392640
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET xfft512K410_1AMD
allfft	DP			OFFSET xfft512K410_2AMD
allfft	DP			OFFSET xfft512K410_3AMD
allfft	DP			OFFSET xfft512K410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			128, GAP2_10_4+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	10180000, 524288, 0.0319,	1392640
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET xfft512K210_1AMD
allfft	DP			OFFSET xfft512K210_2AMD
allfft	DP			OFFSET xfft512K210_3AMD
allfft	DP			OFFSET xfft512K210_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			128, GAP2_10_2+1
allfft	DD			127, 1, 128, 1, 1, 0
	DD	10180000, 524288, 0.0319,	1392640
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET xfft512K110_1AMD
	DP			OFFSET xfft512K110_2AMD
	DP			OFFSET xfft512K110_3AMD
	DP			OFFSET xfft512K110_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			128, GAP2_10_1+1
	DD			127, 1, 128, 1, 1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET xfft512K411_1, OFFSET xfft512K411_2
allfft	DP			OFFSET xfft512K411_3, OFFSET xfft512K411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			64, GAP2_11_4+1
allfft	DD			63, 1, 64, 1, 1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			16384+7*128	;; scratch area size
allfft	DP			OFFSET xfft512K211_1, OFFSET xfft512K211_2
allfft	DP			OFFSET xfft512K211_3, OFFSET xfft512K211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			64, GAP2_11_2+1
allfft	DD			63, 1, 64, 1, 1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			8192+7*128	;; scratch area size
allfft	DP			OFFSET xfft512K111_1, OFFSET xfft512K111_2
allfft	DP			OFFSET xfft512K111_3, OFFSET xfft512K111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			64, GAP2_11_1+1
allfft	DD			63, 1, 64, 1, 1, 0
allfft	DD	10180000, 524288, 0.0319,	1170432
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft512K_1, OFFSET xfft512K_2
allfft	DP			OFFSET xfft512K_3, OFFSET xfft512K_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			64, GAP2_11_4+1
allfft	DD			63, 1, 64, 1, 1, 0
	DD	10180000, 524288, 0.0319,	1392640
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET xfft512K410_1, OFFSET xfft512K410_2
	DP			OFFSET xfft512K410_3, OFFSET xfft512K410_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			128, GAP2_10_4+1
	DD			127, 1, 128, 1, 1, 0
	DD	10180000, 524288, 0.0319,	1392640
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			32768+15*128	;; scratch area size
	DP			OFFSET xfft512K210_1, OFFSET xfft512K210_2
	DP			OFFSET xfft512K210_3, OFFSET xfft512K210_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			128, GAP2_10_2+1
	DD			127, 1, 128, 1, 1, 0
allfft	DD	10180000, 524288, 0.0319,	1392640
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET xfft512K110_1, OFFSET xfft512K110_2
allfft	DP			OFFSET xfft512K110_3, OFFSET xfft512K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			128, GAP2_10_1+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	10180000, 524288, 0.0319,	1392640
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET xfft512K010_1, OFFSET xfft512K010_2
allfft	DP			OFFSET xfft512K010_3, OFFSET xfft512K010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			128, GAP2_10_0+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1214848
allfft	DD			RPFW+512*65536+2*8;; Flags, min_l2_cache, clm
allfft	DD			40960+4*128	;; scratch area size
allfft	DP			OFFSET xfft640K812_1AMD
allfft	DP			OFFSET xfft640K812_2AMD
allfft	DP			OFFSET xfft640K812_3AMD
allfft	DP			OFFSET xfft640K812_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			40, GAP2_12_8+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1214848
allfft	DD			RPFW+512*65536+2*4;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET xfft640K412_1AMD
allfft	DP			OFFSET xfft640K412_2AMD
allfft	DP			OFFSET xfft640K412_3AMD
allfft	DP			OFFSET xfft640K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			40, GAP2_12_4+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
	DD	12650000, 655360, 0.0410,	1214848
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			10240+4*128	;; scratch area size
	DP			OFFSET xfft640K212_1AMD
	DP			OFFSET xfft640K212_2AMD
	DP			OFFSET xfft640K212_3AMD
	DP			OFFSET xfft640K212_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			40, GAP2_12_2+1
	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1214848
allfft	DD			RPFW+512*65536+2*1;; Flags, min_l2_cache, clm
allfft	DD			5120+4*128	;; scratch area size
allfft	DP			OFFSET xfft640K112_1AMD
allfft	DP			OFFSET xfft640K112_2AMD
allfft	DP			OFFSET xfft640K112_3AMD
allfft	DP			OFFSET xfft640K112_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			40, GAP2_12_1+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1437440
allfft	DD			RPFW+2*8	;; Flags, min_l2_cache, clm
allfft	DD			81920+9*128	;; scratch area size
allfft	DP			OFFSET xfft640K811_1AMD
allfft	DP			OFFSET xfft640K811_2AMD
allfft	DP			OFFSET xfft640K811_3AMD
allfft	DP			OFFSET xfft640K811_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			80, GAP2_11_8+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1437440
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			40960+9*128	;; scratch area size
allfft	DP			OFFSET xfft640K411_1AMD
allfft	DP			OFFSET xfft640K411_2AMD
allfft	DP			OFFSET xfft640K411_3AMD
allfft	DP			OFFSET xfft640K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			80, GAP2_11_4+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1214848
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET xfft640K412_1, OFFSET xfft640K412_2
allfft	DP			OFFSET xfft640K412_3, OFFSET xfft640K412_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			40, GAP2_12_4+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
	DD	12650000, 655360, 0.0410,	1214848
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			10240+4*128	;; scratch area size
	DP			OFFSET xfft640K212_1, OFFSET xfft640K212_2
	DP			OFFSET xfft640K212_3, OFFSET xfft640K212_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			40, GAP2_12_2+1
	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1437440
allfft	DD			512*65536+2*8	;; Flags, min_l2_cache, clm
allfft	DD			81920+9*128	;; scratch area size
allfft	DP			OFFSET xfft640K811_1, OFFSET xfft640K811_2
allfft	DP			OFFSET xfft640K811_3, OFFSET xfft640K811_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			80, GAP2_11_8+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
	DD	12650000, 655360, 0.0410,	1437440
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			40960+9*128	;; scratch area size
	DP			OFFSET xfft640K411_1, OFFSET xfft640K411_2
	DP			OFFSET xfft640K411_3, OFFSET xfft640K411_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			80, GAP2_11_4+1
	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1437440
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET xfft640K211_1, OFFSET xfft640K211_2
allfft	DP			OFFSET xfft640K211_3, OFFSET xfft640K211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			80, GAP2_11_2+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
	DD	12650000, 655360, 0.0410,	1437440
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			10240+9*128	;; scratch area size
	DP			OFFSET xfft640K111_1, OFFSET xfft640K111_2
	DP			OFFSET xfft640K111_3, OFFSET xfft640K111_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			80, GAP2_11_1+1
	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1437440
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			10240+9*128	;; scratch area size
allfft	DP			OFFSET xfft640K9111_1, OFFSET xfft640K9111_2
allfft	DP			OFFSET xfft640K9111_3, OFFSET xfft640K9111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			80, GAP2_11_1+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1730816
allfft	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET xfft640K410_1, OFFSET xfft640K410_2
allfft	DP			OFFSET xfft640K410_3, OFFSET xfft640K410_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			160, GAP2_10_4+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1730816
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET xfft640K210_1, OFFSET xfft640K210_2
allfft	DP			OFFSET xfft640K210_3, OFFSET xfft640K210_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			160, GAP2_10_2+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	1730816
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET xfft640K110_1, OFFSET xfft640K110_2
allfft	DP			OFFSET xfft640K110_3, OFFSET xfft640K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			160, GAP2_10_1+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	12650000, 655360, 0.0410,	3123840
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft640K18_1, OFFSET xfft640K18_2
allfft	DP			OFFSET xfft640K18_3, OFFSET xfft640K18_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			8		;; FFT levels done in pass2
allfft	DD			640, GAP2_8_1+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	1413504
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft768K412_1AMD
allfft	DP			OFFSET xfft768K412_2AMD
allfft	DP			OFFSET xfft768K412_3AMD
allfft	DP			OFFSET xfft768K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			48, GAP2_12_4+1
allfft	DD			47, 1, 48, 1, 1, 0
	DD	15070000, 786432, 0.0507,	1413504
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET xfft768K212_1AMD
	DP			OFFSET xfft768K212_2AMD
	DP			OFFSET xfft768K212_3AMD
	DP			OFFSET xfft768K212_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			48, GAP2_12_2+1
	DD			47, 1, 48, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	1413504
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			6144+5*128	;; scratch area size
allfft	DP			OFFSET xfft768K112_1AMD
allfft	DP			OFFSET xfft768K112_2AMD
allfft	DP			OFFSET xfft768K112_3AMD
allfft	DP			OFFSET xfft768K112_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			48, GAP2_12_1+1
allfft	DD			47, 1, 48, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET xfft768K411_1AMD
allfft	DP			OFFSET xfft768K411_2AMD
allfft	DP			OFFSET xfft768K411_3AMD
allfft	DP			OFFSET xfft768K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			96, GAP2_11_4+1
allfft	DD			95, 1, 96, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET xfft768K211_1AMD
allfft	DP			OFFSET xfft768K211_2AMD
allfft	DP			OFFSET xfft768K211_3AMD
allfft	DP			OFFSET xfft768K211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			96, GAP2_11_2+1
allfft	DD			95, 1, 96, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	1413504
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft768K412_1, OFFSET xfft768K412_2
allfft	DP			OFFSET xfft768K412_3, OFFSET xfft768K412_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			48, GAP2_12_4+1
allfft	DD			47, 1, 48, 1, 1, 0
	DD	15070000, 786432, 0.0507,	1413504
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET xfft768K212_1, OFFSET xfft768K212_2
	DP			OFFSET xfft768K212_3, OFFSET xfft768K212_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			48, GAP2_12_2+1
	DD			47, 1, 48, 1, 1, 0
	DD	15070000, 786432, 0.0507,	1703552
	DD			CELE_D*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			49152+11*128	;; scratch area size
	DP			OFFSET xfft768K411_1, OFFSET xfft768K411_2
	DP			OFFSET xfft768K411_3, OFFSET xfft768K411_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			96, GAP2_11_4+1
	DD			95, 1, 96, 1, 1, 0
	DD	15070000, 786432, 0.0507,	1703552
	DD			WILLI*65536+2*2	;; Flags, min_l2_cache, clm
	DD			24576+11*128	;; scratch area size
	DP			OFFSET xfft768K211_1, OFFSET xfft768K211_2
	DP			OFFSET xfft768K211_3, OFFSET xfft768K211_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			96, GAP2_11_2+1
	DD			95, 1, 96, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET xfft768K111_1, OFFSET xfft768K111_2
allfft	DP			OFFSET xfft768K111_3, OFFSET xfft768K111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			96, GAP2_11_1+1
allfft	DD			95, 1, 96, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET xfft768K9111_1, OFFSET xfft768K9111_2
allfft	DP			OFFSET xfft768K9111_3, OFFSET xfft768K9111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			96, GAP2_11_1+1
allfft	DD			95, 1, 96, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET xfft768K011_1, OFFSET xfft768K011_2
allfft	DP			OFFSET xfft768K011_3, OFFSET xfft768K011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			96, GAP2_11_0+1
allfft	DD			95, 1, 96, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	1703552
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET xfft768K9011_1, OFFSET xfft768K9011_2
allfft	DP			OFFSET xfft768K9011_3, OFFSET xfft768K9011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			96, GAP2_11_0+1
allfft	DD			95, 1, 96, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	2066816
allfft	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET xfft768K410_1, OFFSET xfft768K410_2
allfft	DP			OFFSET xfft768K410_3, OFFSET xfft768K410_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			192, GAP2_10_4+1
allfft	DD			191, 1, 192, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	2066816
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET xfft768K210_1, OFFSET xfft768K210_2
allfft	DP			OFFSET xfft768K210_3, OFFSET xfft768K210_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			192, GAP2_10_2+1
allfft	DD			191, 1, 192, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	2066816
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET xfft768K110_1, OFFSET xfft768K110_2
allfft	DP			OFFSET xfft768K110_3, OFFSET xfft768K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			192, GAP2_10_1+1
allfft	DD			191, 1, 192, 1, 1, 0
	DD	15070000, 786432, 0.0507,	2066816
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET xfft768K010_1, OFFSET xfft768K010_2
	DP			OFFSET xfft768K010_3, OFFSET xfft768K010_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			192, GAP2_10_0+1
	DD			191, 1, 192, 1, 1, 0
allfft	DD	15070000, 786432, 0.0507,	3746560
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft768K18_1, OFFSET xfft768K18_2
allfft	DP			OFFSET xfft768K18_3, OFFSET xfft768K18_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			8		;; FFT levels done in pass2
allfft	DD			768, GAP2_8_1+1
allfft	DD			767, 1, 768, 1, 1, 0
allfft	DD	17550000, 917504, 0.0607,	1612416
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET xfft896K412_1AMD
allfft	DP			OFFSET xfft896K412_2AMD
allfft	DP			OFFSET xfft896K412_3AMD
allfft	DP			OFFSET xfft896K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			56, GAP2_12_4+1
allfft	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
	DD	17550000, 917504, 0.0607,	1612416
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			14336+6*128	;; scratch area size
	DP			OFFSET xfft896K212_1AMD
	DP			OFFSET xfft896K212_2AMD
	DP			OFFSET xfft896K212_3AMD
	DP			OFFSET xfft896K212_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			56, GAP2_12_2+1
	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
	DD	17550000, 917504, 0.0607,	1612416
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			7168+6*128	;; scratch area size
	DP			OFFSET xfft896K112_1AMD
	DP			OFFSET xfft896K112_2AMD
	DP			OFFSET xfft896K112_3AMD
	DP			OFFSET xfft896K112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			56, GAP2_12_1+1
	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
allfft	DD	17550000, 917504, 0.0607,	1969920
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			57344+13*128	;; scratch area size
allfft	DP			OFFSET xfft896K411_1AMD
allfft	DP			OFFSET xfft896K411_2AMD
allfft	DP			OFFSET xfft896K411_3AMD
allfft	DP			OFFSET xfft896K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			112, GAP2_11_4+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	17550000, 917504, 0.0607,	1969920
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET xfft896K211_1AMD
allfft	DP			OFFSET xfft896K211_2AMD
allfft	DP			OFFSET xfft896K211_3AMD
allfft	DP			OFFSET xfft896K211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			112, GAP2_11_2+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	17550000, 917504, 0.0607,	1612416
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET xfft896K412_1, OFFSET xfft896K412_2
allfft	DP			OFFSET xfft896K412_3, OFFSET xfft896K412_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			56, GAP2_12_4+1
allfft	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
	DD	17550000, 917504, 0.0607,	1612416
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			14336+6*128	;; scratch area size
	DP			OFFSET xfft896K212_1, OFFSET xfft896K212_2
	DP			OFFSET xfft896K212_3, OFFSET xfft896K212_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			56, GAP2_12_2+1
	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
	DD	17550000, 917504, 0.0607,	1969920
	DD			CELE_D*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			57344+13*128	;; scratch area size
	DP			OFFSET xfft896K411_1, OFFSET xfft896K411_2
	DP			OFFSET xfft896K411_3, OFFSET xfft896K411_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			112, GAP2_11_4+1
	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
	DD	17550000, 917504, 0.0607,	1969920
	DD			WILLI*65536+2*2	;; Flags, min_l2_cache, clm
	DD			28672+13*128	;; scratch area size
	DP			OFFSET xfft896K211_1, OFFSET xfft896K211_2
	DP			OFFSET xfft896K211_3, OFFSET xfft896K211_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			112, GAP2_11_2+1
	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	17550000, 917504, 0.0607,	1969920
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET xfft896K111_1, OFFSET xfft896K111_2
allfft	DP			OFFSET xfft896K111_3, OFFSET xfft896K111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			112, GAP2_11_1+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	17550000, 917504, 0.0607,	1969920
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET xfft896K9111_1, OFFSET xfft896K9111_2
allfft	DP			OFFSET xfft896K9111_3, OFFSET xfft896K9111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			112, GAP2_11_1+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
	DD	17550000, 917504, 0.0607,	1969920
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			14336+13*128	;; scratch area size
	DP			OFFSET xfft896K011_1, OFFSET xfft896K011_2
	DP			OFFSET xfft896K011_3, OFFSET xfft896K011_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			112, GAP2_11_0+1
	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	17550000, 917504, 0.0607,	1969920
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET xfft896K9011_1, OFFSET xfft896K9011_2
allfft	DP			OFFSET xfft896K9011_3, OFFSET xfft896K9011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			112, GAP2_11_0+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	17550000, 917504, 0.0607,	2402816
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET xfft896K210_1, OFFSET xfft896K210_2
allfft	DP			OFFSET xfft896K210_3, OFFSET xfft896K210_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			224, GAP2_10_2+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	17550000, 917504, 0.0607,	2402816
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET xfft896K110_1, OFFSET xfft896K110_2
allfft	DP			OFFSET xfft896K110_3, OFFSET xfft896K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			224, GAP2_10_1+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	17550000, 917504, 0.0607,	2402816
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET xfft896K010_1, OFFSET xfft896K010_2
allfft	DP			OFFSET xfft896K010_3, OFFSET xfft896K010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			224, GAP2_10_0+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	17550000, 917504, 0.0607,	4369280
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET xfft896K18_1, OFFSET xfft896K18_2
allfft	DP			OFFSET xfft896K18_3, OFFSET xfft896K18_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			8		;; FFT levels done in pass2
allfft	DD			896, GAP2_8_1+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	1809408
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET xfft1024K412_1AMD
allfft	DP			OFFSET xfft1024K412_2AMD
allfft	DP			OFFSET xfft1024K412_3AMD
allfft	DP			OFFSET xfft1024K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			64, GAP2_12_4+1
allfft	DD			63, 1, 64, 1, 1, 0
	DD	20050000, 1048576, 0.0676,	1809408
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			16384+7*128	;; scratch area size
	DP			OFFSET xfft1024K212_1AMD
	DP			OFFSET xfft1024K212_2AMD
	DP			OFFSET xfft1024K212_3AMD
	DP			OFFSET xfft1024K212_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			64, GAP2_12_2+1
	DD			63, 1, 64, 1, 1, 0
	DD	20050000, 1048576, 0.0676,	1809408
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			8192+7*128	;; scratch area size
	DP			OFFSET xfft1024K112_1AMD
	DP			OFFSET xfft1024K112_2AMD
	DP			OFFSET xfft1024K112_3AMD
	DP			OFFSET xfft1024K112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			64, GAP2_12_1+1
	DD			63, 1, 64, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	2234368
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET xfft1024K411_1AMD
allfft	DP			OFFSET xfft1024K411_2AMD
allfft	DP			OFFSET xfft1024K411_3AMD
allfft	DP			OFFSET xfft1024K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			128, GAP2_11_4+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	2234368
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET xfft1024K211_1AMD
allfft	DP			OFFSET xfft1024K211_2AMD
allfft	DP			OFFSET xfft1024K211_3AMD
allfft	DP			OFFSET xfft1024K211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			128, GAP2_11_2+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	1809408
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET xfft1024K412_1, OFFSET xfft1024K412_2
allfft	DP			OFFSET xfft1024K412_3, OFFSET xfft1024K412_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			64, GAP2_12_4+1
allfft	DD			63, 1, 64, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	1809408
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			16384+7*128	;; scratch area size
allfft	DP			OFFSET xfft1024K212_1, OFFSET xfft1024K212_2
allfft	DP			OFFSET xfft1024K212_3, OFFSET xfft1024K212_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			64, GAP2_12_2+1
allfft	DD			63, 1, 64, 1, 1, 0
	DD	20050000, 1048576, 0.0676,	1809408
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			8192+7*128	;; scratch area size
	DP			OFFSET xfft1024K112_1, OFFSET xfft1024K112_2
	DP			OFFSET xfft1024K112_3, OFFSET xfft1024K112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			64, GAP2_12_1+1
	DD			63, 1, 64, 1, 1, 0
	DD	20050000, 1048576, 0.0676,	2234368
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET xfft1024K411_1, OFFSET xfft1024K411_2
	DP			OFFSET xfft1024K411_3, OFFSET xfft1024K411_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			128, GAP2_11_4+1
	DD			127, 1, 128, 1, 1, 0
	DD	20050000, 1048576, 0.0676,	2234368
	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
	DD			32768+15*128	;; scratch area size
	DP			OFFSET xfft1024K211_1, OFFSET xfft1024K211_2
	DP			OFFSET xfft1024K211_3, OFFSET xfft1024K211_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			128, GAP2_11_2+1
	DD			127, 1, 128, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	2234368
allfft	DD			705*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET xfft1024K111_1, OFFSET xfft1024K111_2
allfft	DP			OFFSET xfft1024K111_3, OFFSET xfft1024K111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			128, GAP2_11_1+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	2234368
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET xfft1024K9111_1, OFFSET xfft1024K9111_2
allfft	DP			OFFSET xfft1024K9111_3, OFFSET xfft1024K9111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			128, GAP2_11_1+1
allfft	DD			127, 1, 128, 1, 1, 0
	DD	20050000, 1048576, 0.0676,	2234368
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET xfft1024K011_1, OFFSET xfft1024K011_2
	DP			OFFSET xfft1024K011_3, OFFSET xfft1024K011_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			128, GAP2_11_0+1
	DD			127, 1, 128, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	2234368
allfft	DD			600*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET xfft1024K9011_1, OFFSET xfft1024K9011_2
allfft	DP			OFFSET xfft1024K9011_3, OFFSET xfft1024K9011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			128, GAP2_11_0+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	2734080
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET xfft1024K210_1, OFFSET xfft1024K210_2
allfft	DP			OFFSET xfft1024K210_3, OFFSET xfft1024K210_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			256, GAP2_10_2+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	2734080
allfft	DD			200*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft1024K110_1, OFFSET xfft1024K110_2
allfft	DP			OFFSET xfft1024K110_3, OFFSET xfft1024K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			256, GAP2_10_1+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	2734080
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft1024K010_1, OFFSET xfft1024K010_2
allfft	DP			OFFSET xfft1024K010_3, OFFSET xfft1024K010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			256, GAP2_10_0+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	20050000, 1048576, 0.0676,	4975104
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET xfft1024K18_1, OFFSET xfft1024K18_2
allfft	DP			OFFSET xfft1024K18_3, OFFSET xfft1024K18_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			8		;; FFT levels done in pass2
allfft	DD			1024, GAP2_8_1+1
allfft	DD			1023, 1, 1024, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2099584
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET xfft1280K413_1AMD
allfft	DP			OFFSET xfft1280K413_2AMD
allfft	DP			OFFSET xfft1280K413_3AMD
allfft	DP			OFFSET xfft1280K413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			40, GAP2_13_4+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2207488
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			40960+9*128	;; scratch area size
allfft	DP			OFFSET xfft1280K412_1AMD
allfft	DP			OFFSET xfft1280K412_2AMD
allfft	DP			OFFSET xfft1280K412_3AMD
allfft	DP			OFFSET xfft1280K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			80, GAP2_12_4+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2207488
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET xfft1280K212_1AMD
allfft	DP			OFFSET xfft1280K212_2AMD
allfft	DP			OFFSET xfft1280K212_3AMD
allfft	DP			OFFSET xfft1280K212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			80, GAP2_12_2+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
	DD	24930000, 1310720, 0.0892,	2207488
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			10240+9*128	;; scratch area size
	DP			OFFSET xfft1280K112_1AMD
	DP			OFFSET xfft1280K112_2AMD
	DP			OFFSET xfft1280K112_3AMD
	DP			OFFSET xfft1280K112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			80, GAP2_12_1+1
	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET xfft1280K411_1AMD
allfft	DP			OFFSET xfft1280K411_2AMD
allfft	DP			OFFSET xfft1280K411_3AMD
allfft	DP			OFFSET xfft1280K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			160, GAP2_11_4+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			RPFW+256*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET xfft1280K211_1AMD
allfft	DP			OFFSET xfft1280K211_2AMD
allfft	DP			OFFSET xfft1280K211_3AMD
allfft	DP			OFFSET xfft1280K211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			160, GAP2_11_2+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET xfft1280K111_1AMD
allfft	DP			OFFSET xfft1280K111_2AMD
allfft	DP			OFFSET xfft1280K111_3AMD
allfft	DP			OFFSET xfft1280K111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			160, GAP2_11_1+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2099584
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			20480+4*128	;; scratch area size
allfft	DP			OFFSET xfft1280K413_1, OFFSET xfft1280K413_2
allfft	DP			OFFSET xfft1280K413_3, OFFSET xfft1280K413_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			40, GAP2_13_4+1
allfft	DD			39, 1, 32*2048+8, 1, 1, 0
	DD	24930000, 1310720, 0.0892,	2207488
	DD			CELE_D*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			40960+9*128	;; scratch area size
	DP			OFFSET xfft1280K412_1, OFFSET xfft1280K412_2
	DP			OFFSET xfft1280K412_3, OFFSET xfft1280K412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			80, GAP2_12_4+1
	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2207488
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET xfft1280K212_1, OFFSET xfft1280K212_2
allfft	DP			OFFSET xfft1280K212_3, OFFSET xfft1280K212_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			80, GAP2_12_2+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2207488
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			10240+9*128	;; scratch area size
allfft	DP			OFFSET xfft1280K112_1, OFFSET xfft1280K112_2
allfft	DP			OFFSET xfft1280K112_3, OFFSET xfft1280K112_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			80, GAP2_12_1+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET xfft1280K411_1, OFFSET xfft1280K411_2
allfft	DP			OFFSET xfft1280K411_3, OFFSET xfft1280K411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			160, GAP2_11_4+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
	DD	24930000, 1310720, 0.0892,	2769152
	DD			WILLI*65536+2*2	;; Flags, min_l2_cache, clm
	DD			40960+19*128	;; scratch area size
	DP			OFFSET xfft1280K211_1, OFFSET xfft1280K211_2
	DP			OFFSET xfft1280K211_3, OFFSET xfft1280K211_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			160, GAP2_11_2+1
	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET xfft1280K111_1, OFFSET xfft1280K111_2
allfft	DP			OFFSET xfft1280K111_3, OFFSET xfft1280K111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			160, GAP2_11_1+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET xfft1280K9111_1, OFFSET xfft1280K9111_2
allfft	DP			OFFSET xfft1280K9111_3, OFFSET xfft1280K9111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			160, GAP2_11_1+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
	DD	24930000, 1310720, 0.0892,	2769152
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			20480+19*128	;; scratch area size
	DP			OFFSET xfft1280K011_1, OFFSET xfft1280K011_2
	DP			OFFSET xfft1280K011_3, OFFSET xfft1280K011_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			160, GAP2_11_0+1
	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	2769152
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET xfft1280K9011_1, OFFSET xfft1280K9011_2
allfft	DP			OFFSET xfft1280K9011_3, OFFSET xfft1280K9011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			160, GAP2_11_0+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	3409536
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET xfft1280K110_1, OFFSET xfft1280K110_2
allfft	DP			OFFSET xfft1280K110_3, OFFSET xfft1280K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			320, GAP2_10_1+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	24930000, 1310720, 0.0892,	3409536
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET xfft1280K010_1, OFFSET xfft1280K010_2
allfft	DP			OFFSET xfft1280K010_3, OFFSET xfft1280K010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			320, GAP2_10_0+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	2429312
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft1536K413_1AMD
allfft	DP			OFFSET xfft1536K413_2AMD
allfft	DP			OFFSET xfft1536K413_3AMD
allfft	DP			OFFSET xfft1536K413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			48, GAP2_13_4+1
allfft	DD			47, 1, 48, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	2604672
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET xfft1536K412_1AMD
allfft	DP			OFFSET xfft1536K412_2AMD
allfft	DP			OFFSET xfft1536K412_3AMD
allfft	DP			OFFSET xfft1536K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			96, GAP2_12_4+1
allfft	DD			95, 1, 96, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	2604672
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET xfft1536K212_1AMD
allfft	DP			OFFSET xfft1536K212_2AMD
allfft	DP			OFFSET xfft1536K212_3AMD
allfft	DP			OFFSET xfft1536K212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			96, GAP2_12_2+1
allfft	DD			95, 1, 96, 1, 1, 0
	DD	29690000, 1572864, 0.113,	2604672
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET xfft1536K112_1AMD
	DP			OFFSET xfft1536K112_2AMD
	DP			OFFSET xfft1536K112_3AMD
	DP			OFFSET xfft1536K112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			96, GAP2_12_1+1
	DD			95, 1, 96, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536K411_1AMD
allfft	DP			OFFSET xfft1536K411_2AMD
allfft	DP			OFFSET xfft1536K411_3AMD
allfft	DP			OFFSET xfft1536K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_4+1
allfft	DD			191, 1, 192, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536K211_1AMD
allfft	DP			OFFSET xfft1536K211_2AMD
allfft	DP			OFFSET xfft1536K211_3AMD
allfft	DP			OFFSET xfft1536K211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_2+1
allfft	DD			191, 1, 192, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536K111_1AMD
allfft	DP			OFFSET xfft1536K111_2AMD
allfft	DP			OFFSET xfft1536K111_3AMD
allfft	DP			OFFSET xfft1536K111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_1+1
allfft	DD			191, 1, 192, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	2429312
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft1536K413_1, OFFSET xfft1536K413_2
allfft	DP			OFFSET xfft1536K413_3, OFFSET xfft1536K413_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			48, GAP2_13_4+1
allfft	DD			47, 1, 48, 1, 1, 0
	DD	29690000, 1572864, 0.113,	2604672
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			49152+11*128	;; scratch area size
	DP			OFFSET xfft1536K412_1, OFFSET xfft1536K412_2
	DP			OFFSET xfft1536K412_3, OFFSET xfft1536K412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			96, GAP2_12_4+1
	DD			95, 1, 96, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	2604672
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET xfft1536K212_1, OFFSET xfft1536K212_2
allfft	DP			OFFSET xfft1536K212_3, OFFSET xfft1536K212_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			96, GAP2_12_2+1
allfft	DD			95, 1, 96, 1, 1, 0
	DD	29690000, 1572864, 0.113,	2604672
	DD			CELE_D*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET xfft1536K112_1, OFFSET xfft1536K112_2
	DP			OFFSET xfft1536K112_3, OFFSET xfft1536K112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			96, GAP2_12_1+1
	DD			95, 1, 96, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536K411_1, OFFSET xfft1536K411_2
allfft	DP			OFFSET xfft1536K411_3, OFFSET xfft1536K411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_4+1
allfft	DD			191, 1, 192, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536K211_1, OFFSET xfft1536K211_2
allfft	DP			OFFSET xfft1536K211_3, OFFSET xfft1536K211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_2+1
allfft	DD			191, 1, 192, 1, 1, 0
	DD	29690000, 1572864, 0.113,	3301760
	DD			WILLI*65536+2*1	;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET xfft1536K111_1, OFFSET xfft1536K111_2
	DP			OFFSET xfft1536K111_3, OFFSET xfft1536K111_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			192, GAP2_11_1+1
	DD			191, 1, 192, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536K9111_1, OFFSET xfft1536K9111_2
allfft	DP			OFFSET xfft1536K9111_3, OFFSET xfft1536K9111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_1+1
allfft	DD			191, 1, 192, 1, 1, 0
	DD	29690000, 1572864, 0.113,	3301760
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET xfft1536K011_1, OFFSET xfft1536K011_2
	DP			OFFSET xfft1536K011_3, OFFSET xfft1536K011_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			192, GAP2_11_0+1
	DD			191, 1, 192, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	3301760
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536K9011_1, OFFSET xfft1536K9011_2
allfft	DP			OFFSET xfft1536K9011_3, OFFSET xfft1536K9011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_0+1
allfft	DD			191, 1, 192, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	4081280
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET xfft1536K110_1, OFFSET xfft1536K110_2
allfft	DP			OFFSET xfft1536K110_3, OFFSET xfft1536K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			384, GAP2_10_1+1
allfft	DD			383, 1, 384, 1, 1, 0
allfft	DD	29690000, 1572864, 0.113,	4081280
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET xfft1536K010_1, OFFSET xfft1536K010_2
allfft	DP			OFFSET xfft1536K010_3, OFFSET xfft1536K010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			384, GAP2_10_0+1
allfft	DD			383, 1, 384, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	2759296
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET xfft1792K413_1AMD
allfft	DP			OFFSET xfft1792K413_2AMD
allfft	DP			OFFSET xfft1792K413_3AMD
allfft	DP			OFFSET xfft1792K413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			56, GAP2_13_4+1
allfft	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3002112
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			57344+13*128	;; scratch area size
allfft	DP			OFFSET xfft1792K412_1AMD
allfft	DP			OFFSET xfft1792K412_2AMD
allfft	DP			OFFSET xfft1792K412_3AMD
allfft	DP			OFFSET xfft1792K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			112, GAP2_12_4+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3002112
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET xfft1792K212_1AMD
allfft	DP			OFFSET xfft1792K212_2AMD
allfft	DP			OFFSET xfft1792K212_3AMD
allfft	DP			OFFSET xfft1792K212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			112, GAP2_12_2+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
	DD	34560000, 1835008, 0.135,	3002112
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			14336+13*128	;; scratch area size
	DP			OFFSET xfft1792K112_1AMD
	DP			OFFSET xfft1792K112_2AMD
	DP			OFFSET xfft1792K112_3AMD
	DP			OFFSET xfft1792K112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			112, GAP2_12_1+1
	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			114688+27*128	;; scratch area size
allfft	DP			OFFSET xfft1792K411_1AMD
allfft	DP			OFFSET xfft1792K411_2AMD
allfft	DP			OFFSET xfft1792K411_3AMD
allfft	DP			OFFSET xfft1792K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			224, GAP2_11_4+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET xfft1792K211_1AMD
allfft	DP			OFFSET xfft1792K211_2AMD
allfft	DP			OFFSET xfft1792K211_3AMD
allfft	DP			OFFSET xfft1792K211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			224, GAP2_11_2+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET xfft1792K111_1AMD
allfft	DP			OFFSET xfft1792K111_2AMD
allfft	DP			OFFSET xfft1792K111_3AMD
allfft	DP			OFFSET xfft1792K111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			224, GAP2_11_1+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET xfft1792K011_1AMD
allfft	DP			OFFSET xfft1792K011_2AMD
allfft	DP			OFFSET xfft1792K011_3AMD
allfft	DP			OFFSET xfft1792K011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			224, GAP2_11_0+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	2759296
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			28672+6*128	;; scratch area size
allfft	DP			OFFSET xfft1792K413_1, OFFSET xfft1792K413_2
allfft	DP			OFFSET xfft1792K413_3, OFFSET xfft1792K413_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			56, GAP2_13_4+1
allfft	DD			55, 1, (32*2048+16)*2048+8, 1, 1, 0
	DD	34560000, 1835008, 0.135,	3002112
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			57344+13*128	;; scratch area size
	DP			OFFSET xfft1792K412_1, OFFSET xfft1792K412_2
	DP			OFFSET xfft1792K412_3, OFFSET xfft1792K412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			112, GAP2_12_4+1
	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3002112
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET xfft1792K212_1, OFFSET xfft1792K212_2
allfft	DP			OFFSET xfft1792K212_3, OFFSET xfft1792K212_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			112, GAP2_12_2+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
	DD	34560000, 1835008, 0.135,	3002112
	DD			CELE_D*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			14336+13*128	;; scratch area size
	DP			OFFSET xfft1792K112_1, OFFSET xfft1792K112_2
	DP			OFFSET xfft1792K112_3, OFFSET xfft1792K112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			112, GAP2_12_1+1
	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3002112
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET xfft1792K9112_1, OFFSET xfft1792K9112_2
allfft	DP			OFFSET xfft1792K9112_3, OFFSET xfft1792K9112_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			112, GAP2_12_1+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3002112
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			14336+13*128	;; scratch area size
allfft	DP			OFFSET xfft1792K9012_1, OFFSET xfft1792K9012_2
allfft	DP			OFFSET xfft1792K9012_3, OFFSET xfft1792K9012_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			112, GAP2_12_0+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			114688+27*128	;; scratch area size
allfft	DP			OFFSET xfft1792K411_1, OFFSET xfft1792K411_2
allfft	DP			OFFSET xfft1792K411_3, OFFSET xfft1792K411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			224, GAP2_11_4+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET xfft1792K211_1, OFFSET xfft1792K211_2
allfft	DP			OFFSET xfft1792K211_3, OFFSET xfft1792K211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			224, GAP2_11_2+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
	DD	34560000, 1835008, 0.135,	3834368
	DD			WILLI*65536+2*1	;; Flags, min_l2_cache, clm
	DD			28672+27*128	;; scratch area size
	DP			OFFSET xfft1792K111_1, OFFSET xfft1792K111_2
	DP			OFFSET xfft1792K111_3, OFFSET xfft1792K111_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			224, GAP2_11_1+1
	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET xfft1792K9111_1, OFFSET xfft1792K9111_2
allfft	DP			OFFSET xfft1792K9111_3, OFFSET xfft1792K9111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			224, GAP2_11_1+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
	DD	34560000, 1835008, 0.135,	3834368
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			28672+27*128	;; scratch area size
	DP			OFFSET xfft1792K011_1, OFFSET xfft1792K011_2
	DP			OFFSET xfft1792K011_3, OFFSET xfft1792K011_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			224, GAP2_11_0+1
	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	3834368
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET xfft1792K9011_1, OFFSET xfft1792K9011_2
allfft	DP			OFFSET xfft1792K9011_3, OFFSET xfft1792K9011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			224, GAP2_11_0+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	4753024
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET xfft1792K110_1, OFFSET xfft1792K110_2
allfft	DP			OFFSET xfft1792K110_3, OFFSET xfft1792K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			448, GAP2_10_1+1
allfft	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	34560000, 1835008, 0.135,	4753024
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET xfft1792K010_1, OFFSET xfft1792K010_2
allfft	DP			OFFSET xfft1792K010_3, OFFSET xfft1792K010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			448, GAP2_10_0+1
allfft	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	3087360
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET xfft2048K413_1AMD
allfft	DP			OFFSET xfft2048K413_2AMD
allfft	DP			OFFSET xfft2048K413_3AMD
allfft	DP			OFFSET xfft2048K413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			64, GAP2_13_4+1
allfft	DD			63, 1, 64, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	3397632
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET xfft2048K412_1AMD
allfft	DP			OFFSET xfft2048K412_2AMD
allfft	DP			OFFSET xfft2048K412_3AMD
allfft	DP			OFFSET xfft2048K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			128, GAP2_12_4+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	3397632
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET xfft2048K212_1AMD
allfft	DP			OFFSET xfft2048K212_2AMD
allfft	DP			OFFSET xfft2048K212_3AMD
allfft	DP			OFFSET xfft2048K212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			128, GAP2_12_2+1
allfft	DD			127, 1, 128, 1, 1, 0
	DD	39500000, 2097152, 0.155,	3397632
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET xfft2048K112_1AMD
	DP			OFFSET xfft2048K112_2AMD
	DP			OFFSET xfft2048K112_3AMD
	DP			OFFSET xfft2048K112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			128, GAP2_12_1+1
	DD			127, 1, 128, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048K411_1AMD
allfft	DP			OFFSET xfft2048K411_2AMD
allfft	DP			OFFSET xfft2048K411_3AMD
allfft	DP			OFFSET xfft2048K411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_4+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048K211_1AMD
allfft	DP			OFFSET xfft2048K211_2AMD
allfft	DP			OFFSET xfft2048K211_3AMD
allfft	DP			OFFSET xfft2048K211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_2+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048K111_1AMD
allfft	DP			OFFSET xfft2048K111_2AMD
allfft	DP			OFFSET xfft2048K111_3AMD
allfft	DP			OFFSET xfft2048K111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_1+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048K011_1AMD
allfft	DP			OFFSET xfft2048K011_2AMD
allfft	DP			OFFSET xfft2048K011_3AMD
allfft	DP			OFFSET xfft2048K011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_0+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	3087360
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET xfft2048K413_1, OFFSET xfft2048K413_2
allfft	DP			OFFSET xfft2048K413_3, OFFSET xfft2048K413_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			64, GAP2_13_4+1
allfft	DD			63, 1, 64, 1, 1, 0
	DD	39500000, 2097152, 0.155,	3397632
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET xfft2048K412_1, OFFSET xfft2048K412_2
	DP			OFFSET xfft2048K412_3, OFFSET xfft2048K412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			128, GAP2_12_4+1
	DD			127, 1, 128, 1, 1, 0
	DD	39500000, 2097152, 0.155,	3397632
	DD			CELE_D*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET xfft2048K112_1, OFFSET xfft2048K112_2
	DP			OFFSET xfft2048K112_3, OFFSET xfft2048K112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			128, GAP2_12_1+1
	DD			127, 1, 128, 1, 1, 0
	DD	39500000, 2097152, 0.155,	3397632
	DD			WILLI*65536+2*2	;; Flags, min_l2_cache, clm
	DD			32768+15*128	;; scratch area size
	DP			OFFSET xfft2048K212_1, OFFSET xfft2048K212_2
	DP			OFFSET xfft2048K212_3, OFFSET xfft2048K212_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			128, GAP2_12_2+1
	DD			127, 1, 128, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	3397632
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET xfft2048K9112_1, OFFSET xfft2048K9112_2
allfft	DP			OFFSET xfft2048K9112_3, OFFSET xfft2048K9112_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			128, GAP2_12_1+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	3397632
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET xfft2048K9012_1, OFFSET xfft2048K9012_2
allfft	DP			OFFSET xfft2048K9012_3, OFFSET xfft2048K9012_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			128, GAP2_12_0+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048K411_1, OFFSET xfft2048K411_2
allfft	DP			OFFSET xfft2048K411_3, OFFSET xfft2048K411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_4+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048K211_1, OFFSET xfft2048K211_2
allfft	DP			OFFSET xfft2048K211_3, OFFSET xfft2048K211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_2+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048K111_1, OFFSET xfft2048K111_2
allfft	DP			OFFSET xfft2048K111_3, OFFSET xfft2048K111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_1+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048K9111_1, OFFSET xfft2048K9111_2
allfft	DP			OFFSET xfft2048K9111_3, OFFSET xfft2048K9111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_1+1
allfft	DD			255, 1, 256, 1, 1, 0
	DD	39500000, 2097152, 0.155,	4362240
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET xfft2048K011_1, OFFSET xfft2048K011_2
	DP			OFFSET xfft2048K011_3, OFFSET xfft2048K011_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			256, GAP2_11_0+1
	DD			255, 1, 256, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	4362240
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048K9011_1, OFFSET xfft2048K9011_2
allfft	DP			OFFSET xfft2048K9011_3, OFFSET xfft2048K9011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_0+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	5416960
allfft	DD			128*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET xfft2048K110_1, OFFSET xfft2048K110_2
allfft	DP			OFFSET xfft2048K110_3, OFFSET xfft2048K110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			512, GAP2_10_1+1
allfft	DD			511, 1, 512, 1, 1, 0
allfft	DD	39500000, 2097152, 0.155,	5416960
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET xfft2048K010_1, OFFSET xfft2048K010_2
allfft	DP			OFFSET xfft2048K010_3, OFFSET xfft2048K010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			512, GAP2_10_0+1
allfft	DD			511, 1, 512, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	3747584
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			40960+9*128	;; scratch area size
allfft	DP			OFFSET xfft2560K413_1AMD
allfft	DP			OFFSET xfft2560K413_2AMD
allfft	DP			OFFSET xfft2560K413_3AMD
allfft	DP			OFFSET xfft2560K413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			80, GAP2_13_4+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	3747584
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			20480+9*128	;; scratch area size
allfft	DP			OFFSET xfft2560K213_1AMD
allfft	DP			OFFSET xfft2560K213_2AMD
allfft	DP			OFFSET xfft2560K213_3AMD
allfft	DP			OFFSET xfft2560K213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			80, GAP2_13_2+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
	DD	49100000, 2621440, 0.204,	3747584
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			10240+9*128	;; scratch area size
	DP			OFFSET xfft2560K113_1AMD
	DP			OFFSET xfft2560K113_2AMD
	DP			OFFSET xfft2560K113_3AMD
	DP			OFFSET xfft2560K113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			80, GAP2_13_1+1
	DD			79, 1, 64*2048+16, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	4194560
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET xfft2560K412_1AMD
allfft	DP			OFFSET xfft2560K412_2AMD
allfft	DP			OFFSET xfft2560K412_3AMD
allfft	DP			OFFSET xfft2560K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			160, GAP2_12_4+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	4194560
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET xfft2560K212_1AMD
allfft	DP			OFFSET xfft2560K212_2AMD
allfft	DP			OFFSET xfft2560K212_3AMD
allfft	DP			OFFSET xfft2560K212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			160, GAP2_12_2+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
	DD	49100000, 2621440, 0.204,	4194560
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			20480+19*128	;; scratch area size
	DP			OFFSET xfft2560K112_1AMD
	DP			OFFSET xfft2560K112_2AMD
	DP			OFFSET xfft2560K112_3AMD
	DP			OFFSET xfft2560K112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			160, GAP2_12_1+1
	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET xfft2560K211_1AMD
allfft	DP			OFFSET xfft2560K211_2AMD
allfft	DP			OFFSET xfft2560K211_3AMD
allfft	DP			OFFSET xfft2560K211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			320, GAP2_11_2+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET xfft2560K111_1AMD
allfft	DP			OFFSET xfft2560K111_2AMD
allfft	DP			OFFSET xfft2560K111_3AMD
allfft	DP			OFFSET xfft2560K111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			320, GAP2_11_1+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET xfft2560K011_1AMD
allfft	DP			OFFSET xfft2560K011_2AMD
allfft	DP			OFFSET xfft2560K011_3AMD
allfft	DP			OFFSET xfft2560K011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			320, GAP2_11_0+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	3747584
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			40960+9*128	;; scratch area size
allfft	DP			OFFSET xfft2560K413_1, OFFSET xfft2560K413_2
allfft	DP			OFFSET xfft2560K413_3, OFFSET xfft2560K413_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			80, GAP2_13_4+1
allfft	DD			79, 1, 64*2048+16, 1, 1, 0
	DD	49100000, 2621440, 0.204,	4194560
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			81920+19*128	;; scratch area size
	DP			OFFSET xfft2560K412_1, OFFSET xfft2560K412_2
	DP			OFFSET xfft2560K412_3, OFFSET xfft2560K412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			160, GAP2_12_4+1
	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	4194560
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET xfft2560K212_1, OFFSET xfft2560K212_2
allfft	DP			OFFSET xfft2560K212_3, OFFSET xfft2560K212_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			160, GAP2_12_2+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
	DD	49100000, 2621440, 0.204,	4194560
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			20480+19*128	;; scratch area size
	DP			OFFSET xfft2560K112_1, OFFSET xfft2560K112_2
	DP			OFFSET xfft2560K112_3, OFFSET xfft2560K112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			160, GAP2_12_1+1
	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET xfft2560K211_1, OFFSET xfft2560K211_2
allfft	DP			OFFSET xfft2560K211_3, OFFSET xfft2560K211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			320, GAP2_11_2+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET xfft2560K111_1, OFFSET xfft2560K111_2
allfft	DP			OFFSET xfft2560K111_3, OFFSET xfft2560K111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			320, GAP2_11_1+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	5430912
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET xfft2560K011_1, OFFSET xfft2560K011_2
allfft	DP			OFFSET xfft2560K011_3, OFFSET xfft2560K011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			320, GAP2_11_0+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	49100000, 2621440, 0.204,	6768768
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft2560K010_1, OFFSET xfft2560K010_2
allfft	DP			OFFSET xfft2560K010_3, OFFSET xfft2560K010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			640, GAP2_10_0+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	4406912
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET xfft3072K413_1AMD
allfft	DP			OFFSET xfft3072K413_2AMD
allfft	DP			OFFSET xfft3072K413_3AMD
allfft	DP			OFFSET xfft3072K413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			96, GAP2_13_4+1
allfft	DD			95, 1, 96, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	4406912
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET xfft3072K213_1AMD
allfft	DP			OFFSET xfft3072K213_2AMD
allfft	DP			OFFSET xfft3072K213_3AMD
allfft	DP			OFFSET xfft3072K213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			96, GAP2_13_2+1
allfft	DD			95, 1, 96, 1, 1, 0
	DD	58520000, 3145728, 0.259,	4406912
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET xfft3072K113_1AMD
	DP			OFFSET xfft3072K113_2AMD
	DP			OFFSET xfft3072K113_3AMD
	DP			OFFSET xfft3072K113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			96, GAP2_13_1+1
	DD			95, 1, 96, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	4989312
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET xfft3072K412_1AMD
allfft	DP			OFFSET xfft3072K412_2AMD
allfft	DP			OFFSET xfft3072K412_3AMD
allfft	DP			OFFSET xfft3072K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			192, GAP2_12_4+1
allfft	DD			191, 1, 192, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	4989312
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET xfft3072K212_1AMD
allfft	DP			OFFSET xfft3072K212_2AMD
allfft	DP			OFFSET xfft3072K212_3AMD
allfft	DP			OFFSET xfft3072K212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			192, GAP2_12_2+1
allfft	DD			191, 1, 192, 1, 1, 0
	DD	58520000, 3145728, 0.259,	4989312
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET xfft3072K112_1AMD
	DP			OFFSET xfft3072K112_2AMD
	DP			OFFSET xfft3072K112_3AMD
	DP			OFFSET xfft3072K112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			192, GAP2_12_1+1
	DD			191, 1, 192, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072K211_1AMD
allfft	DP			OFFSET xfft3072K211_2AMD
allfft	DP			OFFSET xfft3072K211_3AMD
allfft	DP			OFFSET xfft3072K211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_2+1
allfft	DD			383, 1, 384, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072K111_1AMD
allfft	DP			OFFSET xfft3072K111_2AMD
allfft	DP			OFFSET xfft3072K111_3AMD
allfft	DP			OFFSET xfft3072K111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_1+1
allfft	DD			383, 1, 384, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072K011_1AMD
allfft	DP			OFFSET xfft3072K011_2AMD
allfft	DP			OFFSET xfft3072K011_3AMD
allfft	DP			OFFSET xfft3072K011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_0+1
allfft	DD			383, 1, 384, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	4406912
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET xfft3072K413_1, OFFSET xfft3072K413_2
allfft	DP			OFFSET xfft3072K413_3, OFFSET xfft3072K413_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			96, GAP2_13_4+1
allfft	DD			95, 1, 96, 1, 1, 0
	DD	58520000, 3145728, 0.259,	4989312
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			98304+23*128	;; scratch area size
	DP			OFFSET xfft3072K412_1, OFFSET xfft3072K412_2
	DP			OFFSET xfft3072K412_3, OFFSET xfft3072K412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			192, GAP2_12_4+1
	DD			191, 1, 192, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	4989312
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET xfft3072K212_1, OFFSET xfft3072K212_2
allfft	DP			OFFSET xfft3072K212_3, OFFSET xfft3072K212_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			192, GAP2_12_2+1
allfft	DD			191, 1, 192, 1, 1, 0
	DD	58520000, 3145728, 0.259,	4989312
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET xfft3072K112_1, OFFSET xfft3072K112_2
	DP			OFFSET xfft3072K112_3, OFFSET xfft3072K112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			192, GAP2_12_1+1
	DD			191, 1, 192, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072K211_1, OFFSET xfft3072K211_2
allfft	DP			OFFSET xfft3072K211_3, OFFSET xfft3072K211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_2+1
allfft	DD			383, 1, 384, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072K111_1, OFFSET xfft3072K111_2
allfft	DP			OFFSET xfft3072K111_3, OFFSET xfft3072K111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_1+1
allfft	DD			383, 1, 384, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	6495872
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072K011_1, OFFSET xfft3072K011_2
allfft	DP			OFFSET xfft3072K011_3, OFFSET xfft3072K011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_0+1
allfft	DD			383, 1, 384, 1, 1, 0
allfft	DD	58520000, 3145728, 0.259,	8112384
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft3072K010_1, OFFSET xfft3072K010_2
allfft	DP			OFFSET xfft3072K010_3, OFFSET xfft3072K010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			768, GAP2_10_0+1
allfft	DD			767, 1, 768, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	5066496
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			57344+13*128	;; scratch area size
allfft	DP			OFFSET xfft3584K413_1AMD
allfft	DP			OFFSET xfft3584K413_2AMD
allfft	DP			OFFSET xfft3584K413_3AMD
allfft	DP			OFFSET xfft3584K413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			112, GAP2_13_4+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	5066496
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			28672+13*128	;; scratch area size
allfft	DP			OFFSET xfft3584K213_1AMD
allfft	DP			OFFSET xfft3584K213_2AMD
allfft	DP			OFFSET xfft3584K213_3AMD
allfft	DP			OFFSET xfft3584K213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			112, GAP2_13_2+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
	DD	68130000, 3670016, 0.323,	5066496
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			14336+13*128	;; scratch area size
	DP			OFFSET xfft3584K113_1AMD
	DP			OFFSET xfft3584K113_2AMD
	DP			OFFSET xfft3584K113_3AMD
	DP			OFFSET xfft3584K113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			112, GAP2_13_1+1
	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	5784064
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			114688+27*128	;; scratch area size
allfft	DP			OFFSET xfft3584K412_1AMD
allfft	DP			OFFSET xfft3584K412_2AMD
allfft	DP			OFFSET xfft3584K412_3AMD
allfft	DP			OFFSET xfft3584K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			224, GAP2_12_4+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	5784064
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET xfft3584K212_1AMD
allfft	DP			OFFSET xfft3584K212_2AMD
allfft	DP			OFFSET xfft3584K212_3AMD
allfft	DP			OFFSET xfft3584K212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			224, GAP2_12_2+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
	DD	68130000, 3670016, 0.323,	5784064
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			28672+27*128	;; scratch area size
	DP			OFFSET xfft3584K112_1AMD
	DP			OFFSET xfft3584K112_2AMD
	DP			OFFSET xfft3584K112_3AMD
	DP			OFFSET xfft3584K112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			224, GAP2_12_1+1
	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			114688+55*128	;; scratch area size
allfft	DP			OFFSET xfft3584K211_1AMD
allfft	DP			OFFSET xfft3584K211_2AMD
allfft	DP			OFFSET xfft3584K211_3AMD
allfft	DP			OFFSET xfft3584K211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			448, GAP2_11_2+1
allfft	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET xfft3584K111_1AMD
allfft	DP			OFFSET xfft3584K111_2AMD
allfft	DP			OFFSET xfft3584K111_3AMD
allfft	DP			OFFSET xfft3584K111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			448, GAP2_11_1+1
allfft	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET xfft3584K011_1AMD
allfft	DP			OFFSET xfft3584K011_2AMD
allfft	DP			OFFSET xfft3584K011_3AMD
allfft	DP			OFFSET xfft3584K011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			448, GAP2_11_0+1
allfft	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	5066496
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			57344+13*128	;; scratch area size
allfft	DP			OFFSET xfft3584K413_1, OFFSET xfft3584K413_2
allfft	DP			OFFSET xfft3584K413_3, OFFSET xfft3584K413_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			112, GAP2_13_4+1
allfft	DD			111, 1, (64*2048+32)*2048+16, 1, 1, 0
	DD	68130000, 3670016, 0.323,	5784064
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			114688+27*128	;; scratch area size
	DP			OFFSET xfft3584K412_1, OFFSET xfft3584K412_2
	DP			OFFSET xfft3584K412_3, OFFSET xfft3584K412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			224, GAP2_12_4+1
	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
	DD	68130000, 3670016, 0.323,	5784064
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			57344+27*128	;; scratch area size
	DP			OFFSET xfft3584K212_1, OFFSET xfft3584K212_2
	DP			OFFSET xfft3584K212_3, OFFSET xfft3584K212_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			224, GAP2_12_2+1
	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
	DD	68130000, 3670016, 0.323,	5784064
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			28672+27*128	;; scratch area size
	DP			OFFSET xfft3584K112_1, OFFSET xfft3584K112_2
	DP			OFFSET xfft3584K112_3, OFFSET xfft3584K112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			224, GAP2_12_1+1
	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			114688+55*128	;; scratch area size
allfft	DP			OFFSET xfft3584K211_1, OFFSET xfft3584K211_2
allfft	DP			OFFSET xfft3584K211_3, OFFSET xfft3584K211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			448, GAP2_11_2+1
allfft	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET xfft3584K111_1, OFFSET xfft3584K111_2
allfft	DP			OFFSET xfft3584K111_3, OFFSET xfft3584K111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			448, GAP2_11_1+1
allfft	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	7560832
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			57344+55*128	;; scratch area size
allfft	DP			OFFSET xfft3584K011_1, OFFSET xfft3584K011_2
allfft	DP			OFFSET xfft3584K011_3, OFFSET xfft3584K011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			448, GAP2_11_0+1
allfft	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	68130000, 3670016, 0.323,	9456000
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET xfft3584K010_1, OFFSET xfft3584K010_2
allfft	DP			OFFSET xfft3584K010_3, OFFSET xfft3584K010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			896, GAP2_10_0+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	5724160
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET xfft4096K413_1AMD
allfft	DP			OFFSET xfft4096K413_2AMD
allfft	DP			OFFSET xfft4096K413_3AMD
allfft	DP			OFFSET xfft4096K413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			128, GAP2_13_4+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	5724160
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET xfft4096K213_1AMD
allfft	DP			OFFSET xfft4096K213_2AMD
allfft	DP			OFFSET xfft4096K213_3AMD
allfft	DP			OFFSET xfft4096K213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			128, GAP2_13_2+1
allfft	DD			127, 1, 128, 1, 1, 0
	DD	77910000, 4194304, 0.382,	5724160
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET xfft4096K113_1AMD
	DP			OFFSET xfft4096K113_2AMD
	DP			OFFSET xfft4096K113_3AMD
	DP			OFFSET xfft4096K113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			128, GAP2_13_1+1
	DD			127, 1, 128, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	6574080
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET xfft4096K412_1AMD
allfft	DP			OFFSET xfft4096K412_2AMD
allfft	DP			OFFSET xfft4096K412_3AMD
allfft	DP			OFFSET xfft4096K412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			256, GAP2_12_4+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	6574080
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET xfft4096K212_1AMD
allfft	DP			OFFSET xfft4096K212_2AMD
allfft	DP			OFFSET xfft4096K212_3AMD
allfft	DP			OFFSET xfft4096K212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			256, GAP2_12_2+1
allfft	DD			255, 1, 256, 1, 1, 0
	DD	77910000, 4194304, 0.382,	6574080
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET xfft4096K112_1AMD
	DP			OFFSET xfft4096K112_2AMD
	DP			OFFSET xfft4096K112_3AMD
	DP			OFFSET xfft4096K112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, GAP2_12_1+1
	DD			255, 1, 256, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096K211_1AMD
allfft	DP			OFFSET xfft4096K211_2AMD
allfft	DP			OFFSET xfft4096K211_3AMD
allfft	DP			OFFSET xfft4096K211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_2+1
allfft	DD			511, 1, 512, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096K111_1AMD
allfft	DP			OFFSET xfft4096K111_2AMD
allfft	DP			OFFSET xfft4096K111_3AMD
allfft	DP			OFFSET xfft4096K111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_1+1
allfft	DD			511, 1, 512, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096K011_1AMD
allfft	DP			OFFSET xfft4096K011_2AMD
allfft	DP			OFFSET xfft4096K011_3AMD
allfft	DP			OFFSET xfft4096K011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_0+1
allfft	DD			511, 1, 512, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	5724160
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET xfft4096K413_1, OFFSET xfft4096K413_2
allfft	DP			OFFSET xfft4096K413_3, OFFSET xfft4096K413_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			128, GAP2_13_4+1
allfft	DD			127, 1, 128, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	5724160
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET xfft4096K213_1, OFFSET xfft4096K213_2
allfft	DP			OFFSET xfft4096K213_3, OFFSET xfft4096K213_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			128, GAP2_13_2+1
allfft	DD			127, 1, 128, 1, 1, 0
	DD	77910000, 4194304, 0.382,	6574080
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			131072+31*128	;; scratch area size
	DP			OFFSET xfft4096K412_1, OFFSET xfft4096K412_2
	DP			OFFSET xfft4096K412_3, OFFSET xfft4096K412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, GAP2_12_4+1
	DD			255, 1, 256, 1, 1, 0
	DD	77910000, 4194304, 0.382,	6574080
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			65536+31*128	;; scratch area size
	DP			OFFSET xfft4096K212_1, OFFSET xfft4096K212_2
	DP			OFFSET xfft4096K212_3, OFFSET xfft4096K212_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, GAP2_12_2+1
	DD			255, 1, 256, 1, 1, 0
	DD	77910000, 4194304, 0.382,	6574080
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET xfft4096K112_1, OFFSET xfft4096K112_2
	DP			OFFSET xfft4096K112_3, OFFSET xfft4096K112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, GAP2_12_1+1
	DD			255, 1, 256, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096K211_1, OFFSET xfft4096K211_2
allfft	DP			OFFSET xfft4096K211_3, OFFSET xfft4096K211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_2+1
allfft	DD			511, 1, 512, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096K111_1, OFFSET xfft4096K111_2
allfft	DP			OFFSET xfft4096K111_3, OFFSET xfft4096K111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_1+1
allfft	DD			511, 1, 512, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	8617984
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096K011_1, OFFSET xfft4096K011_2
allfft	DP			OFFSET xfft4096K011_3, OFFSET xfft4096K011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_0+1
allfft	DD			511, 1, 512, 1, 1, 0
allfft	DD	77910000, 4194304, 0.382,	10782720
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET xfft4096K010_1, OFFSET xfft4096K010_2
allfft	DP			OFFSET xfft4096K010_3, OFFSET xfft4096K010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			1024, GAP2_10_0+1
allfft	DD			1023, 1, 1024, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	7045376
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET xfft5M413_1AMD, OFFSET xfft5M413_2AMD
allfft	DP			OFFSET xfft5M413_3AMD, OFFSET xfft5M413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			160, GAP2_13_4+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	7045376
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET xfft5M213_1AMD, OFFSET xfft5M213_2AMD
allfft	DP			OFFSET xfft5M213_3AMD, OFFSET xfft5M213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			160, GAP2_13_2+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	7045376
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			20480+19*128	;; scratch area size
allfft	DP			OFFSET xfft5M113_1AMD, OFFSET xfft5M113_2AMD
allfft	DP			OFFSET xfft5M113_3AMD, OFFSET xfft5M113_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			160, GAP2_13_1+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
	DD	96830000, 5242880, 0.485,	8167040
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			163840+39*128	;; scratch area size
	DP			OFFSET xfft5M412_1AMD, OFFSET xfft5M412_2AMD
	DP			OFFSET xfft5M412_3AMD, OFFSET xfft5M412_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			320, GAP2_12_4+1
	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	8167040
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET xfft5M212_1AMD, OFFSET xfft5M212_2AMD
allfft	DP			OFFSET xfft5M212_3AMD, OFFSET xfft5M212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			320, GAP2_12_2+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
	DD	96830000, 5242880, 0.485,	8167040
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			40960+39*128	;; scratch area size
	DP			OFFSET xfft5M112_1AMD, OFFSET xfft5M112_2AMD
	DP			OFFSET xfft5M112_3AMD, OFFSET xfft5M112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			320, GAP2_12_1+1
	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	10756224
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft5M111_1AMD, OFFSET xfft5M111_2AMD
allfft	DP			OFFSET xfft5M111_3AMD, OFFSET xfft5M111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			640, GAP2_11_1+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	10756224
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft5M011_1AMD, OFFSET xfft5M011_2AMD
allfft	DP			OFFSET xfft5M011_3AMD, OFFSET xfft5M011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			640, GAP2_11_0+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	7045376
allfft	DD			1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			81920+19*128	;; scratch area size
allfft	DP			OFFSET xfft5M413_1, OFFSET xfft5M413_2
allfft	DP			OFFSET xfft5M413_3, OFFSET xfft5M413_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			160, GAP2_13_4+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	7045376
allfft	DD			1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			40960+19*128	;; scratch area size
allfft	DP			OFFSET xfft5M213_1, OFFSET xfft5M213_2
allfft	DP			OFFSET xfft5M213_3, OFFSET xfft5M213_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			160, GAP2_13_2+1
allfft	DD			159, 1, 128*2048+32, 1, 1, 0
	DD	96830000, 5242880, 0.485,	8167040
	DD			1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			163840+39*128	;; scratch area size
	DP			OFFSET xfft5M412_1, OFFSET xfft5M412_2
	DP			OFFSET xfft5M412_3, OFFSET xfft5M412_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			320, GAP2_12_4+1
	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	8167040
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET xfft5M212_1, OFFSET xfft5M212_2
allfft	DP			OFFSET xfft5M212_3, OFFSET xfft5M212_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			320, GAP2_12_2+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
	DD	96830000, 5242880, 0.485,	8167040
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			40960+39*128	;; scratch area size
	DP			OFFSET xfft5M112_1, OFFSET xfft5M112_2
	DP			OFFSET xfft5M112_3, OFFSET xfft5M112_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			320, GAP2_12_1+1
	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	8167040
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			40960+39*128	;; scratch area size
allfft	DP			OFFSET xfft5M012_1, OFFSET xfft5M012_2
allfft	DP			OFFSET xfft5M012_3, OFFSET xfft5M012_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			320, GAP2_12_0+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	10756224
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft5M111_1, OFFSET xfft5M111_2
allfft	DP			OFFSET xfft5M111_3, OFFSET xfft5M111_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			640, GAP2_11_1+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	96830000, 5242880, 0.485,	10756224
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft5M011_1, OFFSET xfft5M011_2
allfft	DP			OFFSET xfft5M011_3, OFFSET xfft5M011_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			640, GAP2_11_0+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	115300000, 6291456, 0.668,	8364416
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET xfft6M413_1AMD, OFFSET xfft6M413_2AMD
allfft	DP			OFFSET xfft6M413_3AMD, OFFSET xfft6M413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			192, GAP2_13_4+1
allfft	DD			191, 1, 192, 1, 1, 0
allfft	DD	115300000, 6291456, 0.668,	8364416
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET xfft6M213_1AMD, OFFSET xfft6M213_2AMD
allfft	DP			OFFSET xfft6M213_3AMD, OFFSET xfft6M213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			192, GAP2_13_2+1
allfft	DD			191, 1, 192, 1, 1, 0
allfft	DD	115300000, 6291456, 0.668,	8364416
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET xfft6M113_1AMD, OFFSET xfft6M113_2AMD
allfft	DP			OFFSET xfft6M113_3AMD, OFFSET xfft6M113_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			192, GAP2_13_1+1
allfft	DD			191, 1, 192, 1, 1, 0
	DD	115300000, 6291456, 0.668,	9756288
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET xfft6M412_1AMD, OFFSET xfft6M412_2AMD
	DP			OFFSET xfft6M412_3AMD, OFFSET xfft6M412_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, GAP2_12_4+1
	DD			383, 1, 384, 1, 1, 0
allfft	DD	115300000, 6291456, 0.668,	9756288
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET xfft6M212_1AMD, OFFSET xfft6M212_2AMD
allfft	DP			OFFSET xfft6M212_3AMD, OFFSET xfft6M212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			384, GAP2_12_2+1
allfft	DD			383, 1, 384, 1, 1, 0
	DD	115300000, 6291456, 0.668,	9756288
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET xfft6M112_1AMD, OFFSET xfft6M112_2AMD
	DP			OFFSET xfft6M112_3AMD, OFFSET xfft6M112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, GAP2_12_1+1
	DD			383, 1, 384, 1, 1, 0
allfft	DD	115300000, 6291456, 0.668,	12886272
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft6M111_1AMD, OFFSET xfft6M111_2AMD
allfft	DP			OFFSET xfft6M111_3AMD, OFFSET xfft6M111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			768, GAP2_11_1+1
allfft	DD			767, 1, 768, 1, 1, 0
allfft	DD	115300000, 6291456, 0.668,	12886272
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft6M011_1AMD, OFFSET xfft6M011_2AMD
allfft	DP			OFFSET xfft6M011_3AMD, OFFSET xfft6M011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			768, GAP2_11_0+1
allfft	DD			767, 1, 768, 1, 1, 0
allfft	DD	115300000, 6291456, 0.668,	8364416
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET xfft6M413_1, OFFSET xfft6M413_2
allfft	DP			OFFSET xfft6M413_3, OFFSET xfft6M413_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			192, GAP2_13_4+1
allfft	DD			191, 1, 192, 1, 1, 0
allfft	DD	115300000, 6291456, 0.668,	8364416
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET xfft6M213_1, OFFSET xfft6M213_2
allfft	DP			OFFSET xfft6M213_3, OFFSET xfft6M213_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			192, GAP2_13_2+1
allfft	DD			191, 1, 192, 1, 1, 0
	DD	115300000, 6291456, 0.668,	9756288
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET xfft6M412_1, OFFSET xfft6M412_2
	DP			OFFSET xfft6M412_3, OFFSET xfft6M412_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, GAP2_12_4+1
	DD			383, 1, 384, 1, 1, 0
	DD	115300000, 6291456, 0.668,	9756288
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			98304+47*128	;; scratch area size
	DP			OFFSET xfft6M212_1, OFFSET xfft6M212_2
	DP			OFFSET xfft6M212_3, OFFSET xfft6M212_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, GAP2_12_2+1
	DD			383, 1, 384, 1, 1, 0
	DD	115300000, 6291456, 0.668,	9756288
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET xfft6M112_1, OFFSET xfft6M112_2
	DP			OFFSET xfft6M112_3, OFFSET xfft6M112_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, GAP2_12_1+1
	DD			383, 1, 384, 1, 1, 0
allfft	DD	115300000, 6291456, 0.668,	12886272
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft6M111_1, OFFSET xfft6M111_2
allfft	DP			OFFSET xfft6M111_3, OFFSET xfft6M111_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			768, GAP2_11_1+1
allfft	DD			767, 1, 768, 1, 1, 0
allfft	DD	115300000, 6291456, 0.668,	12886272
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft6M011_1, OFFSET xfft6M011_2
allfft	DP			OFFSET xfft6M011_3, OFFSET xfft6M011_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			768, GAP2_11_0+1
allfft	DD			767, 1, 768, 1, 1, 0
allfft	DD	134200000, 7340032, 0.886,	9683456
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			114688+27*128	;; scratch area size
allfft	DP			OFFSET xfft7M413_1AMD, OFFSET xfft7M413_2AMD
allfft	DP			OFFSET xfft7M413_3AMD, OFFSET xfft7M413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			224, GAP2_13_4+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	134200000, 7340032, 0.886,	9683456
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET xfft7M213_1AMD, OFFSET xfft7M213_2AMD
allfft	DP			OFFSET xfft7M213_3AMD, OFFSET xfft7M213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			224, GAP2_13_2+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	134200000, 7340032, 0.886,	9683456
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			28672+27*128	;; scratch area size
allfft	DP			OFFSET xfft7M113_1AMD, OFFSET xfft7M113_2AMD
allfft	DP			OFFSET xfft7M113_3AMD, OFFSET xfft7M113_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			224, GAP2_13_1+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
	DD	134200000, 7340032, 0.886,	11345536
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			229376+55*128	;; scratch area size
	DP			OFFSET xfft7M412_1AMD, OFFSET xfft7M412_2AMD
	DP			OFFSET xfft7M412_3AMD, OFFSET xfft7M412_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			448, GAP2_12_4+1
	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	134200000, 7340032, 0.886,	11345536
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			114688+55*128	;; scratch area size
allfft	DP			OFFSET xfft7M212_1AMD, OFFSET xfft7M212_2AMD
allfft	DP			OFFSET xfft7M212_3AMD, OFFSET xfft7M212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			448, GAP2_12_2+1
allfft	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
	DD	134200000, 7340032, 0.886,	11345536
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			57344+55*128	;; scratch area size
	DP			OFFSET xfft7M112_1AMD, OFFSET xfft7M112_2AMD
	DP			OFFSET xfft7M112_3AMD, OFFSET xfft7M112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			448, GAP2_12_1+1
	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	134200000, 7340032, 0.886,	15016320
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET xfft7M111_1AMD, OFFSET xfft7M111_2AMD
allfft	DP			OFFSET xfft7M111_3AMD, OFFSET xfft7M111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			896, GAP2_11_1+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	134200000, 7340032, 0.886,	15016320
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET xfft7M011_1AMD, OFFSET xfft7M011_2AMD
allfft	DP			OFFSET xfft7M011_3AMD, OFFSET xfft7M011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			896, GAP2_11_0+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	134200000, 7340032, 0.886,	9683456
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			114688+27*128	;; scratch area size
allfft	DP			OFFSET xfft7M413_1, OFFSET xfft7M413_2
allfft	DP			OFFSET xfft7M413_3, OFFSET xfft7M413_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			224, GAP2_13_4+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
allfft	DD	134200000, 7340032, 0.886,	9683456
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			57344+27*128	;; scratch area size
allfft	DP			OFFSET xfft7M213_1, OFFSET xfft7M213_2
allfft	DP			OFFSET xfft7M213_3, OFFSET xfft7M213_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			224, GAP2_13_2+1
allfft	DD			223, 1, (128*2048+64)*2048+32, 1, 1, 0
	DD	134200000, 7340032, 0.886,	11345536
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			229376+55*128	;; scratch area size
	DP			OFFSET xfft7M412_1, OFFSET xfft7M412_2
	DP			OFFSET xfft7M412_3, OFFSET xfft7M412_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			448, GAP2_12_4+1
	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
	DD	134200000, 7340032, 0.886,	11345536
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			114688+55*128	;; scratch area size
	DP			OFFSET xfft7M212_1, OFFSET xfft7M212_2
	DP			OFFSET xfft7M212_3, OFFSET xfft7M212_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			448, GAP2_12_2+1
	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
	DD	134200000, 7340032, 0.886,	11345536
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			57344+55*128	;; scratch area size
	DP			OFFSET xfft7M112_1, OFFSET xfft7M112_2
	DP			OFFSET xfft7M112_3, OFFSET xfft7M112_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			448, GAP2_12_1+1
	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	134200000, 7340032, 0.886,	15016320
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET xfft7M111_1, OFFSET xfft7M111_2
allfft	DP			OFFSET xfft7M111_3, OFFSET xfft7M111_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			896, GAP2_11_1+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	134200000, 7340032, 0.886,	15016320
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET xfft7M011_1, OFFSET xfft7M011_2
allfft	DP			OFFSET xfft7M011_3, OFFSET xfft7M011_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			896, GAP2_11_0+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET xfft8M413_1AMD, OFFSET xfft8M413_2AMD
allfft	DP			OFFSET xfft8M413_3AMD, OFFSET xfft8M413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			256, GAP2_13_4+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET xfft8M213_1AMD, OFFSET xfft8M213_2AMD
allfft	DP			OFFSET xfft8M213_3AMD, OFFSET xfft8M213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			256, GAP2_13_2+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft8M113_1AMD, OFFSET xfft8M113_2AMD
allfft	DP			OFFSET xfft8M113_3AMD, OFFSET xfft8M113_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			256, GAP2_13_1+1
allfft	DD			255, 1, 256, 1, 1, 0
	DD	153400000, 8388608, 1.042,	12926976
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET xfft8M412_1AMD, OFFSET xfft8M412_2AMD
	DP			OFFSET xfft8M412_3AMD, OFFSET xfft8M412_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, GAP2_12_4+1
	DD			511, 1, 512, 1, 1, 0
allfft	DD	153400000, 8388608, 1.042,	12926976
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET xfft8M212_1AMD, OFFSET xfft8M212_2AMD
allfft	DP			OFFSET xfft8M212_3AMD, OFFSET xfft8M212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			512, GAP2_12_2+1
allfft	DD			511, 1, 512, 1, 1, 0
	DD	153400000, 8388608, 1.042,	12926976
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET xfft8M112_1AMD, OFFSET xfft8M112_2AMD
	DP			OFFSET xfft8M112_3AMD, OFFSET xfft8M112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, GAP2_12_1+1
	DD			511, 1, 512, 1, 1, 0
allfft	DD	153400000, 8388608, 1.042,	17129472
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET xfft8M111_1AMD, OFFSET xfft8M111_2AMD
allfft	DP			OFFSET xfft8M111_3AMD, OFFSET xfft8M111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			1024, GAP2_11_1+1
allfft	DD			1023, 1, 1024, 1, 1, 0
allfft	DD	153400000, 8388608, 1.042,	17129472
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET xfft8M011_1AMD, OFFSET xfft8M011_2AMD
allfft	DP			OFFSET xfft8M011_3AMD, OFFSET xfft8M011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			1024, GAP2_11_0+1
allfft	DD			1023, 1, 1024, 1, 1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET xfft8M413_1, OFFSET xfft8M413_2
allfft	DP			OFFSET xfft8M413_3, OFFSET xfft8M413_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			256, GAP2_13_4+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET xfft8M213_1, OFFSET xfft8M213_2
allfft	DP			OFFSET xfft8M213_3, OFFSET xfft8M213_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			256, GAP2_13_2+1
allfft	DD			255, 1, 256, 1, 1, 0
allfft	DD	153400000, 8388608, 1.042,	10997760
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft8M113_1, OFFSET xfft8M113_2
allfft	DP			OFFSET xfft8M113_3, OFFSET xfft8M113_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			256, GAP2_13_1+1
allfft	DD			255, 1, 256, 1, 1, 0
	DD	153400000, 8388608, 1.042,	12926976
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET xfft8M412_1, OFFSET xfft8M412_2
	DP			OFFSET xfft8M412_3, OFFSET xfft8M412_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, GAP2_12_4+1
	DD			511, 1, 512, 1, 1, 0
	DD	153400000, 8388608, 1.042,	12926976
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			131072+63*128	;; scratch area size
	DP			OFFSET xfft8M212_1, OFFSET xfft8M212_2
	DP			OFFSET xfft8M212_3, OFFSET xfft8M212_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, GAP2_12_2+1
	DD			511, 1, 512, 1, 1, 0
	DD	153400000, 8388608, 1.042,	12926976
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET xfft8M112_1, OFFSET xfft8M112_2
	DP			OFFSET xfft8M112_3, OFFSET xfft8M112_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, GAP2_12_1+1
	DD			511, 1, 512, 1, 1, 0
allfft	DD	153400000, 8388608, 1.042,	17129472
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET xfft8M111_1, OFFSET xfft8M111_2
allfft	DP			OFFSET xfft8M111_3, OFFSET xfft8M111_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			1024, GAP2_11_1+1
allfft	DD			1023, 1, 1024, 1, 1, 0
allfft	DD	153400000, 8388608, 1.042,	17129472
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET xfft8M011_1, OFFSET xfft8M011_2
allfft	DP			OFFSET xfft8M011_3, OFFSET xfft8M011_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			1024, GAP2_11_0+1
allfft	DD			1023, 1, 1024, 1, 1, 0
	DD	190700000, 10485760, 1.100,	13639296
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			163840+39*128	;; scratch area size
	DP			OFFSET xfft10M413_1AMD, OFFSET xfft10M413_2AMD
	DP			OFFSET xfft10M413_3AMD, OFFSET xfft10M413_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			320, GAP2_13_4+1
	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	190700000, 10485760, 1.100,	13639296
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET xfft10M213_1AMD, OFFSET xfft10M213_2AMD
allfft	DP			OFFSET xfft10M213_3AMD, OFFSET xfft10M213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			320, GAP2_13_2+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
	DD	190700000, 10485760, 1.100,	13639296
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			40960+39*128	;; scratch area size
	DP			OFFSET xfft10M113_1AMD, OFFSET xfft10M113_2AMD
	DP			OFFSET xfft10M113_3AMD, OFFSET xfft10M113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			320, GAP2_13_1+1
	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			163840+79*128	;; scratch area size
allfft	DP			OFFSET xfft10M212_1AMD, OFFSET xfft10M212_2AMD
allfft	DP			OFFSET xfft10M212_3AMD, OFFSET xfft10M212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			640, GAP2_12_2+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft10M112_1AMD, OFFSET xfft10M112_2AMD
allfft	DP			OFFSET xfft10M112_3AMD, OFFSET xfft10M112_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			640, GAP2_12_1+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft10M012_1AMD, OFFSET xfft10M012_2AMD
allfft	DP			OFFSET xfft10M012_3AMD, OFFSET xfft10M012_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			640, GAP2_12_0+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
	DD	190700000, 10485760, 1.100,	13639296
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			163840+39*128	;; scratch area size
	DP			OFFSET xfft10M413_1, OFFSET xfft10M413_2
	DP			OFFSET xfft10M413_3, OFFSET xfft10M413_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			320, GAP2_13_4+1
	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	190700000, 10485760, 1.100,	13639296
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			81920+39*128	;; scratch area size
allfft	DP			OFFSET xfft10M213_1, OFFSET xfft10M213_2
allfft	DP			OFFSET xfft10M213_3, OFFSET xfft10M213_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			320, GAP2_13_2+1
allfft	DD			319, 1, 256*2048+64, 1, 1, 0
	DD	190700000, 10485760, 1.100,	13639296
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			40960+39*128	;; scratch area size
	DP			OFFSET xfft10M113_1, OFFSET xfft10M113_2
	DP			OFFSET xfft10M113_3, OFFSET xfft10M113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			320, GAP2_13_1+1
	DD			319, 1, 256*2048+64, 1, 1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			163840+79*128	;; scratch area size
allfft	DP			OFFSET xfft10M212_1, OFFSET xfft10M212_2
allfft	DP			OFFSET xfft10M212_3, OFFSET xfft10M212_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			640, GAP2_12_2+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft10M112_1, OFFSET xfft10M112_2
allfft	DP			OFFSET xfft10M112_3, OFFSET xfft10M112_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			640, GAP2_12_1+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	190700000, 10485760, 1.100,	16113792
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft10M012_1, OFFSET xfft10M012_2
allfft	DP			OFFSET xfft10M012_3, OFFSET xfft10M012_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			640, GAP2_12_0+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
	DD	227300000, 12582912, 1.400,	16277120
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET xfft12M413_1AMD, OFFSET xfft12M413_2AMD
	DP			OFFSET xfft12M413_3AMD, OFFSET xfft12M413_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			384, GAP2_13_4+1
	DD			383, 1, 384, 1, 1, 0
allfft	DD	227300000, 12582912, 1.400,	16277120
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET xfft12M213_1AMD, OFFSET xfft12M213_2AMD
allfft	DP			OFFSET xfft12M213_3AMD, OFFSET xfft12M213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			384, GAP2_13_2+1
allfft	DD			383, 1, 384, 1, 1, 0
	DD	227300000, 12582912, 1.400,	16277120
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET xfft12M113_1AMD, OFFSET xfft12M113_2AMD
	DP			OFFSET xfft12M113_3AMD, OFFSET xfft12M113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			384, GAP2_13_1+1
	DD			383, 1, 384, 1, 1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			196608+95*128	;; scratch area size
allfft	DP			OFFSET xfft12M212_1AMD, OFFSET xfft12M212_2AMD
allfft	DP			OFFSET xfft12M212_3AMD, OFFSET xfft12M212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			768, GAP2_12_2+1
allfft	DD			767, 1, 768, 1, 1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft12M112_1AMD, OFFSET xfft12M112_2AMD
allfft	DP			OFFSET xfft12M112_3AMD, OFFSET xfft12M112_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			768, GAP2_12_1+1
allfft	DD			767, 1, 768, 1, 1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft12M012_1AMD, OFFSET xfft12M012_2AMD
allfft	DP			OFFSET xfft12M012_3AMD, OFFSET xfft12M012_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			768, GAP2_12_0+1
allfft	DD			767, 1, 768, 1, 1, 0
	DD	227300000, 12582912, 1.400,	16277120
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET xfft12M413_1, OFFSET xfft12M413_2
	DP			OFFSET xfft12M413_3, OFFSET xfft12M413_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			384, GAP2_13_4+1
	DD			383, 1, 384, 1, 1, 0
	DD	227300000, 12582912, 1.400,	16277120
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			98304+47*128	;; scratch area size
	DP			OFFSET xfft12M213_1, OFFSET xfft12M213_2
	DP			OFFSET xfft12M213_3, OFFSET xfft12M213_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			384, GAP2_13_2+1
	DD			383, 1, 384, 1, 1, 0
	DD	227300000, 12582912, 1.400,	16277120
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET xfft12M113_1, OFFSET xfft12M113_2
	DP			OFFSET xfft12M113_3, OFFSET xfft12M113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			384, GAP2_13_1+1
	DD			383, 1, 384, 1, 1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			196608+95*128	;; scratch area size
allfft	DP			OFFSET xfft12M212_1, OFFSET xfft12M212_2
allfft	DP			OFFSET xfft12M212_3, OFFSET xfft12M212_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			768, GAP2_12_2+1
allfft	DD			767, 1, 768, 1, 1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft12M112_1, OFFSET xfft12M112_2
allfft	DP			OFFSET xfft12M112_3, OFFSET xfft12M112_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			768, GAP2_12_1+1
allfft	DD			767, 1, 768, 1, 1, 0
allfft	DD	227300000, 12582912, 1.400,	19292416
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft12M012_1, OFFSET xfft12M012_2
allfft	DP			OFFSET xfft12M012_3, OFFSET xfft12M012_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			768, GAP2_12_0+1
allfft	DD			767, 1, 768, 1, 1, 0
	DD	264600000, 14680064, 1.700,	18914944
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			229376+55*128	;; scratch area size
	DP			OFFSET xfft14M413_1AMD, OFFSET xfft14M413_2AMD
	DP			OFFSET xfft14M413_3AMD, OFFSET xfft14M413_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			448, GAP2_13_4+1
	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	264600000, 14680064, 1.700,	18914944
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			114688+55*128	;; scratch area size
allfft	DP			OFFSET xfft14M213_1AMD, OFFSET xfft14M213_2AMD
allfft	DP			OFFSET xfft14M213_3AMD, OFFSET xfft14M213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			448, GAP2_13_2+1
allfft	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
	DD	264600000, 14680064, 1.700,	18914944
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			57344+55*128	;; scratch area size
	DP			OFFSET xfft14M113_1AMD, OFFSET xfft14M113_2AMD
	DP			OFFSET xfft14M113_3AMD, OFFSET xfft14M113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			448, GAP2_13_1+1
	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			229376+111*128	;; scratch area size
allfft	DP			OFFSET xfft14M212_1AMD, OFFSET xfft14M212_2AMD
allfft	DP			OFFSET xfft14M212_3AMD, OFFSET xfft14M212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			896, GAP2_12_2+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET xfft14M112_1AMD, OFFSET xfft14M112_2AMD
allfft	DP			OFFSET xfft14M112_3AMD, OFFSET xfft14M112_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			896, GAP2_12_1+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET xfft14M012_1AMD, OFFSET xfft14M012_2AMD
allfft	DP			OFFSET xfft14M012_3AMD, OFFSET xfft14M012_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			896, GAP2_12_0+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
	DD	264600000, 14680064, 1.700,	18914944
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			229376+55*128	;; scratch area size
	DP			OFFSET xfft14M413_1, OFFSET xfft14M413_2
	DP			OFFSET xfft14M413_3, OFFSET xfft14M413_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			448, GAP2_13_4+1
	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
	DD	264600000, 14680064, 1.700,	18914944
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			114688+55*128	;; scratch area size
	DP			OFFSET xfft14M213_1, OFFSET xfft14M213_2
	DP			OFFSET xfft14M213_3, OFFSET xfft14M213_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			448, GAP2_13_2+1
	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
	DD	264600000, 14680064, 1.700,	18914944
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			57344+55*128	;; scratch area size
	DP			OFFSET xfft14M113_1, OFFSET xfft14M113_2
	DP			OFFSET xfft14M113_3, OFFSET xfft14M113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			448, GAP2_13_1+1
	DD			447, 1, (256*2048+128)*2048+64, 1, 1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			229376+111*128	;; scratch area size
allfft	DP			OFFSET xfft14M212_1, OFFSET xfft14M212_2
allfft	DP			OFFSET xfft14M212_3, OFFSET xfft14M212_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			896, GAP2_12_2+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET xfft14M112_1, OFFSET xfft14M112_2
allfft	DP			OFFSET xfft14M112_3, OFFSET xfft14M112_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			896, GAP2_12_1+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	264600000, 14680064, 1.700,	22471040
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET xfft14M012_1, OFFSET xfft14M012_2
allfft	DP			OFFSET xfft14M012_3, OFFSET xfft14M012_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			896, GAP2_12_0+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
	DD	302600000, 16777216, 2.100,	21544960
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET xfft16M413_1AMD, OFFSET xfft16M413_2AMD
	DP			OFFSET xfft16M413_3AMD, OFFSET xfft16M413_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			512, GAP2_13_4+1
	DD			511, 1, 512, 1, 1, 0
allfft	DD	302600000, 16777216, 2.100,	21544960
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET xfft16M213_1AMD, OFFSET xfft16M213_2AMD
allfft	DP			OFFSET xfft16M213_3AMD, OFFSET xfft16M213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			512, GAP2_13_2+1
allfft	DD			511, 1, 512, 1, 1, 0
	DD	302600000, 16777216, 2.100,	21544960
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET xfft16M113_1AMD, OFFSET xfft16M113_2AMD
	DP			OFFSET xfft16M113_3AMD, OFFSET xfft16M113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			512, GAP2_13_1+1
	DD			511, 1, 512, 1, 1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			262144+127*128	;; scratch area size
allfft	DP			OFFSET xfft16M212_1AMD, OFFSET xfft16M212_2AMD
allfft	DP			OFFSET xfft16M212_3AMD, OFFSET xfft16M212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			1024, GAP2_12_2+1
allfft	DD			1023, 1, 1024, 1, 1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET xfft16M112_1AMD, OFFSET xfft16M112_2AMD
allfft	DP			OFFSET xfft16M112_3AMD, OFFSET xfft16M112_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			1024, GAP2_12_1+1
allfft	DD			1023, 1, 1024, 1, 1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET xfft16M012_1AMD, OFFSET xfft16M012_2AMD
allfft	DP			OFFSET xfft16M012_3AMD, OFFSET xfft16M012_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			1024, GAP2_12_0+1
allfft	DD			1023, 1, 1024, 1, 1, 0
	DD	302600000, 16777216, 2.100,	21544960
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET xfft16M413_1, OFFSET xfft16M413_2
	DP			OFFSET xfft16M413_3, OFFSET xfft16M413_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			512, GAP2_13_4+1
	DD			511, 1, 512, 1, 1, 0
	DD	302600000, 16777216, 2.100,	21544960
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			131072+63*128	;; scratch area size
	DP			OFFSET xfft16M213_1, OFFSET xfft16M213_2
	DP			OFFSET xfft16M213_3, OFFSET xfft16M213_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			512, GAP2_13_2+1
	DD			511, 1, 512, 1, 1, 0
	DD	302600000, 16777216, 2.100,	21544960
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET xfft16M113_1, OFFSET xfft16M113_2
	DP			OFFSET xfft16M113_3, OFFSET xfft16M113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			512, GAP2_13_1+1
	DD			511, 1, 512, 1, 1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			262144+127*128	;; scratch area size
allfft	DP			OFFSET xfft16M212_1, OFFSET xfft16M212_2
allfft	DP			OFFSET xfft16M212_3, OFFSET xfft16M212_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			1024, GAP2_12_2+1
allfft	DD			1023, 1, 1024, 1, 1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET xfft16M112_1, OFFSET xfft16M112_2
allfft	DP			OFFSET xfft16M112_3, OFFSET xfft16M112_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			1024, GAP2_12_1+1
allfft	DD			1023, 1, 1024, 1, 1, 0
allfft	DD	302600000, 16777216, 2.100,	25632768
allfft	DD			512*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET xfft16M012_1, OFFSET xfft16M012_2
allfft	DP			OFFSET xfft16M012_3, OFFSET xfft16M012_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			1024, GAP2_12_0+1
allfft	DD			1023, 1, 1024, 1, 1, 0
allfft	DD	376100000, 20971520, 2.700,	26828928
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			327680+79*128	;; scratch area size
allfft	DP			OFFSET xfft20M413_1AMD, OFFSET xfft20M413_2AMD
allfft	DP			OFFSET xfft20M413_3AMD, OFFSET xfft20M413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			640, GAP2_13_4+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
	DD	376100000, 20971520, 2.700,	26828928
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			163840+79*128	;; scratch area size
	DP			OFFSET xfft20M213_1AMD, OFFSET xfft20M213_2AMD
	DP			OFFSET xfft20M213_3AMD, OFFSET xfft20M213_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			640, GAP2_13_2+1
	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	376100000, 20971520, 2.700,	26828928
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft20M113_1AMD, OFFSET xfft20M113_2AMD
allfft	DP			OFFSET xfft20M113_3AMD, OFFSET xfft20M113_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			640, GAP2_13_1+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
	DD	376100000, 20971520, 2.700,	26828928
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			81920+79*128	;; scratch area size
	DP			OFFSET xfft20M013_1AMD, OFFSET xfft20M013_2AMD
	DP			OFFSET xfft20M013_3AMD, OFFSET xfft20M013_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			640, GAP2_13_0+1
	DD			639, 1, 512*2048+128, 1, 1, 0
	DD	376100000, 20971520, 2.700,	26828928
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			163840+79*128	;; scratch area size
	DP			OFFSET xfft20M213_1, OFFSET xfft20M213_2
	DP			OFFSET xfft20M213_3, OFFSET xfft20M213_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			640, GAP2_13_2+1
	DD			639, 1, 512*2048+128, 1, 1, 0
	DD	376100000, 20971520, 2.700,	26828928
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			81920+79*128	;; scratch area size
	DP			OFFSET xfft20M113_1, OFFSET xfft20M113_2
	DP			OFFSET xfft20M113_3, OFFSET xfft20M113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			640, GAP2_13_1+1
	DD			639, 1, 512*2048+128, 1, 1, 0
allfft	DD	376100000, 20971520, 2.700,	26828928
allfft	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			81920+79*128	;; scratch area size
allfft	DP			OFFSET xfft20M013_1, OFFSET xfft20M013_2
allfft	DP			OFFSET xfft20M013_3, OFFSET xfft20M013_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			640, GAP2_13_0+1
allfft	DD			639, 1, 512*2048+128, 1, 1, 0
	DD	448000000, 25165824, 3.300,	32104704
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			196608+95*128	;; scratch area size
	DP			OFFSET xfft24M213_1AMD, OFFSET xfft24M213_2AMD
	DP			OFFSET xfft24M213_3AMD, OFFSET xfft24M213_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			768, GAP2_13_2+1
	DD			767, 1, 768, 1, 1, 0
allfft	DD	448000000, 25165824, 3.300,	32104704
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft24M113_1AMD, OFFSET xfft24M113_2AMD
allfft	DP			OFFSET xfft24M113_3AMD, OFFSET xfft24M113_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			768, GAP2_13_1+1
allfft	DD			767, 1, 768, 1, 1, 0
	DD	448000000, 25165824, 3.300,	32104704
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET xfft24M013_1AMD, OFFSET xfft24M013_2AMD
	DP			OFFSET xfft24M013_3AMD, OFFSET xfft24M013_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			768, GAP2_13_0+1
	DD			767, 1, 768, 1, 1, 0
	DD	448000000, 25165824, 3.300,	32104704
	DD			2048*65536+2*2	;; Flags, min_l2_cache, clm
	DD			196608+95*128	;; scratch area size
	DP			OFFSET xfft24M213_1, OFFSET xfft24M213_2
	DP			OFFSET xfft24M213_3, OFFSET xfft24M213_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			768, GAP2_13_2+1
	DD			767, 1, 768, 1, 1, 0
	DD	448000000, 25165824, 3.300,	32104704
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET xfft24M113_1, OFFSET xfft24M113_2
	DP			OFFSET xfft24M113_3, OFFSET xfft24M113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			768, GAP2_13_1+1
	DD			767, 1, 768, 1, 1, 0
	DD	448000000, 25165824, 3.300,	32104704
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET xfft24M013_1, OFFSET xfft24M013_2
	DP			OFFSET xfft24M013_3, OFFSET xfft24M013_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			768, GAP2_13_0+1
	DD			767, 1, 768, 1, 1, 0
	DD	521500000, 29360128, 3.900,	37380480
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			229376+111*128	;; scratch area size
	DP			OFFSET xfft28M213_1AMD, OFFSET xfft28M213_2AMD
	DP			OFFSET xfft28M213_3AMD, OFFSET xfft28M213_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			896, GAP2_13_2+1
	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	521500000, 29360128, 3.900,	37380480
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			114688+111*128	;; scratch area size
allfft	DP			OFFSET xfft28M113_1AMD, OFFSET xfft28M113_2AMD
allfft	DP			OFFSET xfft28M113_3AMD, OFFSET xfft28M113_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			896, GAP2_13_1+1
allfft	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
	DD	521500000, 29360128, 3.900,	37380480
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			114688+111*128	;; scratch area size
	DP			OFFSET xfft28M013_1AMD, OFFSET xfft28M013_2AMD
	DP			OFFSET xfft28M013_3AMD, OFFSET xfft28M013_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			896, GAP2_13_0+1
	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
	DD	521500000, 29360128, 3.900,	37380480
	DD			2048*65536+2*2	;; Flags, min_l2_cache, clm
	DD			229376+111*128	;; scratch area size
	DP			OFFSET xfft28M213_1, OFFSET xfft28M213_2
	DP			OFFSET xfft28M213_3, OFFSET xfft28M213_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			896, GAP2_13_2+1
	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
	DD	521500000, 29360128, 3.900,	37380480
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			114688+111*128	;; scratch area size
	DP			OFFSET xfft28M113_1, OFFSET xfft28M113_2
	DP			OFFSET xfft28M113_3, OFFSET xfft28M113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			896, GAP2_13_1+1
	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
	DD	521500000, 29360128, 3.900,	37380480
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			114688+111*128	;; scratch area size
	DP			OFFSET xfft28M013_1, OFFSET xfft28M013_2
	DP			OFFSET xfft28M013_3, OFFSET xfft28M013_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			896, GAP2_13_0+1
	DD			895, 1, (512*2048+256)*2048+128, 1, 1, 0
allfft	DD	596000000, 33554432, 4.500,	42639360
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			262144+127*128	;; scratch area size
allfft	DP			OFFSET xfft32M213_1AMD, OFFSET xfft32M213_2AMD
allfft	DP			OFFSET xfft32M213_3AMD, OFFSET xfft32M213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			1024, GAP2_13_2+1
allfft	DD			1023, 1, 1024, 1, 1, 0
	DD	596000000, 33554432, 4.500,	42639360
	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET xfft32M113_1AMD, OFFSET xfft32M113_2AMD
	DP			OFFSET xfft32M113_3AMD, OFFSET xfft32M113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			1024, GAP2_13_1+1
	DD			1023, 1, 1024, 1, 1, 0
	DD	596000000, 33554432, 4.500,	42639360
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET xfft32M013_1AMD, OFFSET xfft32M013_2AMD
	DP			OFFSET xfft32M013_3AMD, OFFSET xfft32M013_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			1024, GAP2_13_0+1
	DD			1023, 1, 1024, 1, 1, 0
	DD	596000000, 33554432, 4.500,	42639360
	DD			2048*65536+2*2	;; Flags, min_l2_cache, clm
	DD			262144+127*128	;; scratch area size
	DP			OFFSET xfft32M213_1, OFFSET xfft32M213_2
	DP			OFFSET xfft32M213_3, OFFSET xfft32M213_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			1024, GAP2_13_2+1
	DD			1023, 1, 1024, 1, 1, 0
	DD	596000000, 33554432, 4.500,	42639360
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET xfft32M113_1, OFFSET xfft32M113_2
	DP			OFFSET xfft32M113_3, OFFSET xfft32M113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			1024, GAP2_13_1+1
	DD			1023, 1, 1024, 1, 1, 0
	DD	596000000, 33554432, 4.500,	42639360
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET xfft32M013_1, OFFSET xfft32M013_2
	DP			OFFSET xfft32M013_3, OFFSET xfft32M013_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			1024, GAP2_13_0+1
	DD			1023, 1, 1024, 1, 1, 0
	DD	0
xjmptablep DD	739,	32,	0.00000111,	1152
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft32p_1, OFFSET xfft32p_2
	DP			OFFSET xfft32p_3, OFFSET xfft32p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			4, 4
	DD			2, 1, 1, 1, 1, 0
	DD	1095,	48,	0.00000144,	1920
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft48p_1, OFFSET xfft48p_2
	DP			OFFSET xfft48p_3, OFFSET xfft48p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			6, 6
	DD			3, 1, 1, 1, 1, 0
	DD	1465,	64,	0.00000178,	2432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft64p_1, OFFSET xfft64p_2
	DP			OFFSET xfft64p_3, OFFSET xfft64p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			8, 8
	DD			4, 1, 1, 1, 1, 0
	DD	2173,	96,	0.00000259,	3328
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft96p_1, OFFSET xfft96p_2
	DP			OFFSET xfft96p_3, OFFSET xfft96p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			12, 12
	DD			6, 3, 1, 1, 1, 0
	DD	2897,	128,	0.00000319,	4352
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft128p_1, OFFSET xfft128p_2
	DP			OFFSET xfft128p_3, OFFSET xfft128p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			16, 16
	DD			8, 4, 1, 2, 1, 0
	DD	4295,	192,	0.00000542,	6784
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft192p_1, OFFSET xfft192p_2
	DP			OFFSET xfft192p_3, OFFSET xfft192p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			24, 24
	DD			12, 6, 3, 3, 1, 0
	DD	5729,	256,	0.00000691,	8576
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft256p_1, OFFSET xfft256p_2
	DP			OFFSET xfft256p_3, OFFSET xfft256p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			32, 32
	DD			16, 8, 4, 4, 1, 0
	DD	8493,	384,	0.0000111,	13312
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft384p_1, OFFSET xfft384p_2
	DP			OFFSET xfft384p_3, OFFSET xfft384p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			48, 48
	DD			24, 12, 3, 6, 1, 0
	DD	11319,	512,	0.0000143,	17152
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft512p_1, OFFSET xfft512p_2
	DP			OFFSET xfft512p_3, OFFSET xfft512p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			64, 64
	DD			32, 16, 4, 2*256+8, 1, 0
	DD	16779,	768,	0.0000260,	26624
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft768p_1, OFFSET xfft768p_2
	DP			OFFSET xfft768p_3, OFFSET xfft768p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			96, 96
	DD			48, 24, 3*256+6, 3*256+12, 1, 0
	DD	22381,	1024,	0.0000349,	34304
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft1024p_1, OFFSET xfft1024p_2
	DP			OFFSET xfft1024p_3, OFFSET xfft1024p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			128, 128
	DD			64, 32, 4*256+8, 4*256+16, 1, 0
	DD	33189,	1536,	0.0000601,	52992
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft1536p_1, OFFSET xfft1536p_2
	DP			OFFSET xfft1536p_3, OFFSET xfft1536p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			192, 192
	DD			96, 48, 12, 6*256+24, 3, 0
	DD	44221,	2048,	0.0000773,	68608
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft2048p_1, OFFSET xfft2048p_2
	DP			OFFSET xfft2048p_3, OFFSET xfft2048p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			256, 256
	DD			128, 64, 16, 8*256+32, 2*256+4, 0
	DD	65519,	3072,	0.000131,	106112
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft3072p_1, OFFSET xfft3072p_2
	DP			OFFSET xfft3072p_3, OFFSET xfft3072p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			384, 384
	DD			192, 96, 6*256+24, 12*256+48, 3*256+3, 0
	DD	87271,	4096,	0.000172,	137216
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft4096p_1, OFFSET xfft4096p_2
	DP			OFFSET xfft4096p_3, OFFSET xfft4096p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			512, 512
	DD			256, 128, 8*256+32, 16*256+64, 4*256+4, 0
	DD	129600,	6144,	0.000291,	211968
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft6144p_1, OFFSET xfft6144p_2
	DP			OFFSET xfft6144p_3, OFFSET xfft6144p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			768, 768
	DD			384, 192, 12*256+48, 24*256+96, 6*256+3, 0
	DD	172400,	8192,	0.000395,	274432
	DD			0		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft8192p_1, OFFSET xfft8192p_2
	DP			OFFSET xfft8192p_3, OFFSET xfft8192p_4
	DP			OFFSET xprctab1	;; Table of add/sub/norm procs
	DD			0		;; FFT levels done in pass2
	DD			1024, 1024
	DD			512, 256, 16*256+64, 32*256+128, 8*256+4, 0
	DD	254800,	12288,	0.000626,	64896
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft12Kp_1AMD, OFFSET xfft12Kp_2AMD
	DP			OFFSET xfft12Kp_3AMD, OFFSET xfft12Kp_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			12, GAP2_8_4+1
	DD			12, 1, 12, 1, 1, 0
	DD	254800,	12288,	0.000626,	64896
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft12Kp_1, OFFSET xfft12Kp_2
	DP			OFFSET xfft12Kp_3, OFFSET xfft12Kp_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			12, GAP2_8_4+1
	DD			12, 1, 12, 1, 1, 0
	DD	339100,	16384,	0.000857,	84224
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft16Kp_1AMD, OFFSET xfft16Kp_2AMD
	DP			OFFSET xfft16Kp_3AMD, OFFSET xfft16Kp_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			16, GAP2_8_4+1
	DD			16, 1, 16, 1, 1, 0
	DD	339100,	16384,	0.000857,	84224
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft16Kp_1, OFFSET xfft16Kp_2
	DP			OFFSET xfft16Kp_3, OFFSET xfft16Kp_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			16, GAP2_8_4+1
	DD			16, 1, 16, 1, 1, 0
	DD	503400,	24576,	0.00135,	123904
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft24Kp_1AMD, OFFSET xfft24Kp_2AMD
	DP			OFFSET xfft24Kp_3AMD, OFFSET xfft24Kp_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			24, GAP2_8_4+1
	DD			24, 1, 24, 1, 1, 0
	DD	503400,	24576,	0.00135,	123904
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft24Kp_1, OFFSET xfft24Kp_2
	DP			OFFSET xfft24Kp_3, OFFSET xfft24Kp_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			24, GAP2_8_4+1
	DD			24, 1, 24, 1, 1, 0
	DD	669600,	32768,	0.00179,	162816
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft32Kp_1AMD, OFFSET xfft32Kp_2AMD
	DP			OFFSET xfft32Kp_3AMD, OFFSET xfft32Kp_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			32, GAP2_8_4+1
	DD			32, 1, 32, 1, 1, 0
	DD	669600,	32768,	0.00179,	162816
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft32Kp_1, OFFSET xfft32Kp_2
	DP			OFFSET xfft32Kp_3, OFFSET xfft32Kp_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			8		;; FFT levels done in pass2
	DD			32, GAP2_8_4+1
	DD			32, 1, 32, 1, 1, 0
	DD	992800, 49152,	0.00303,	145280
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft48Kp_1AMD, OFFSET xfft48Kp_2AMD
	DP			OFFSET xfft48Kp_3AMD, OFFSET xfft48Kp_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			6, GAP2_11_4+1
	DD			6, 1, 6, 1, 1, 0
	DD	992800, 49152,	0.00303,	145280
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft48Kp_1, OFFSET xfft48Kp_2
	DP			OFFSET xfft48Kp_3, OFFSET xfft48Kp_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			6, GAP2_11_4+1
	DD			6, 1, 6, 1, 1, 0
	DD	1320000, 65536,	0.00404,	178560
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft64Kp_1AMD, OFFSET xfft64Kp_2AMD
	DP			OFFSET xfft64Kp_3AMD, OFFSET xfft64Kp_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			8, GAP2_11_4+1
	DD			8, 1, 8, 1, 1, 0
	DD	1320000, 65536,	0.00404,	178560
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft64Kp_1, OFFSET xfft64Kp_2
	DP			OFFSET xfft64Kp_3, OFFSET xfft64Kp_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			8, GAP2_11_4+1
	DD			8, 1, 8, 1, 1, 0
	DD	1962000, 98304, 0.00644,	245632
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft96Kp_1AMD, OFFSET xfft96Kp_2AMD
	DP			OFFSET xfft96Kp_3AMD, OFFSET xfft96Kp_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			12, GAP2_11_4+1
	DD			12, 1, 12, 1, 1, 0
	DD	1962000, 98304, 0.00644,	245632
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft96Kp_1, OFFSET xfft96Kp_2
	DP			OFFSET xfft96Kp_3, OFFSET xfft96Kp_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			12, GAP2_11_4+1
	DD			12, 1, 12, 1, 1, 0
	DD	2610000, 131072, 0.00871,	312064
	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft128Kp_1AMD, OFFSET xfft128Kp_2AMD
	DP			OFFSET xfft128Kp_3AMD, OFFSET xfft128Kp_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			16, GAP2_11_4+1
	DD			16, 1, 16, 1, 1, 0
	DD	2610000, 131072, 0.00871,	312064
	DD			2*4		;; Flags, min_l2_cache, clm
	DD			0		;; scratch area size
	DP			OFFSET xfft128Kp_1, OFFSET xfft128Kp_2
	DP			OFFSET xfft128Kp_3, OFFSET xfft128Kp_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			16, GAP2_11_4+1
	DD			16, 1, 16, 1, 1, 0
allfft	DD	3867000, 196608, 0.0136,	445952
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft192Kp_1AMD, OFFSET xfft192Kp_2AMD
allfft	DP			OFFSET xfft192Kp_3AMD, OFFSET xfft192Kp_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			24, GAP2_11_4+1
allfft	DD			24, 1, 24, 1, 1, 0
allfft	DD	3867000, 196608, 0.0136,	529536
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft192Kp410_1AMD
allfft	DP			OFFSET xfft192Kp410_2AMD
allfft	DP			OFFSET xfft192Kp410_3AMD
allfft	DP			OFFSET xfft192Kp410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			48, GAP2_10_4+1
allfft	DD			48, 1, 48, 1, 1, 0
	DD	3867000, 196608, 0.0136,	529536
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET xfft192Kp210_1AMD
	DP			OFFSET xfft192Kp210_2AMD
	DP			OFFSET xfft192Kp210_3AMD
	DP			OFFSET xfft192Kp210_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			48, GAP2_10_2+1
	DD			48, 1, 48, 1, 1, 0
allfft	DD	3867000, 196608, 0.0136,	529536
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			6144+5*128	;; scratch area size
allfft	DP			OFFSET xfft192Kp110_1AMD
allfft	DP			OFFSET xfft192Kp110_2AMD
allfft	DP			OFFSET xfft192Kp110_3AMD
allfft	DP			OFFSET xfft192Kp110_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			48, GAP2_10_1+1
allfft	DD			48, 1, 48, 1, 1, 0
allfft	DD	3867000, 196608, 0.0136,	445952
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft192Kp_1, OFFSET xfft192Kp_2
allfft	DP			OFFSET xfft192Kp_3, OFFSET xfft192Kp_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			24, GAP2_11_4+1
allfft	DD			24, 1, 24, 1, 1, 0
allfft	DD	3867000, 196608, 0.0136,	529536
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft192Kp410_1, OFFSET xfft192Kp410_2
allfft	DP			OFFSET xfft192Kp410_3, OFFSET xfft192Kp410_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			48, GAP2_10_4+1
allfft	DD			48, 1, 48, 1, 1, 0
allfft	DD	3867000, 196608, 0.0136,	529536
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			12288+5*128	;; scratch area size
allfft	DP			OFFSET xfft192Kp210_1, OFFSET xfft192Kp210_2
allfft	DP			OFFSET xfft192Kp210_3, OFFSET xfft192Kp210_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			48, GAP2_10_2+1
allfft	DD			48, 1, 48, 1, 1, 0
	DD	3867000, 196608, 0.0136,	529536
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			6144+5*128	;; scratch area size
	DP			OFFSET xfft192Kp110_1, OFFSET xfft192Kp110_2
	DP			OFFSET xfft192Kp110_3, OFFSET xfft192Kp110_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			48, GAP2_10_1+1
	DD			48, 1, 48, 1, 1, 0
allfft	DD	5146000, 262144, 0.0179,	579072
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft256Kp_1AMD, OFFSET xfft256Kp_2AMD
allfft	DP			OFFSET xfft256Kp_3AMD, OFFSET xfft256Kp_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			32, GAP2_11_4+1
allfft	DD			32, 1, 32, 1, 1, 0
allfft	DD	5146000, 262144, 0.0179,	697344
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET xfft256Kp410_1AMD
allfft	DP			OFFSET xfft256Kp410_2AMD
allfft	DP			OFFSET xfft256Kp410_3AMD
allfft	DP			OFFSET xfft256Kp410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			64, GAP2_10_4+1
allfft	DD			64, 1, 64, 1, 1, 0
allfft	DD	5146000, 262144, 0.0179,	697344
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			16384+7*128	;; scratch area size
allfft	DP			OFFSET xfft256Kp210_1AMD
allfft	DP			OFFSET xfft256Kp210_2AMD
allfft	DP			OFFSET xfft256Kp210_3AMD
allfft	DP			OFFSET xfft256Kp210_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			64, GAP2_10_2+1
allfft	DD			64, 1, 64, 1, 1, 0
	DD	5146000, 262144, 0.0179,	697344
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			8192+7*128	;; scratch area size
	DP			OFFSET xfft256Kp110_1AMD
	DP			OFFSET xfft256Kp110_2AMD
	DP			OFFSET xfft256Kp110_3AMD
	DP			OFFSET xfft256Kp110_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			64, GAP2_10_1+1
	DD			64, 1, 64, 1, 1, 0
allfft	DD	5146000, 262144, 0.0179,	579072
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft256Kp_1, OFFSET xfft256Kp_2
allfft	DP			OFFSET xfft256Kp_3, OFFSET xfft256Kp_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			32, GAP2_11_4+1
allfft	DD			32, 1, 32, 1, 1, 0
	DD	5146000, 262144, 0.0179,	697344
	DD			256*65536+2*4	;; Flags, min_l2_cache, clm
	DD			32768+7*128	;; scratch area size
	DP			OFFSET xfft256Kp410_1, OFFSET xfft256Kp410_2
	DP			OFFSET xfft256Kp410_3, OFFSET xfft256Kp410_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			64, GAP2_10_4+1
	DD			64, 1, 64, 1, 1, 0
	DD	5146000, 262144, 0.0179,	697344
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			16384+7*128	;; scratch area size
	DP			OFFSET xfft256Kp210_1, OFFSET xfft256Kp210_2
	DP			OFFSET xfft256Kp210_3, OFFSET xfft256Kp210_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			64, GAP2_10_2+1
	DD			64, 1, 64, 1, 1, 0
allfft	DD	5146000, 262144, 0.0179,	697344
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			8192+7*128	;; scratch area size
allfft	DP			OFFSET xfft256Kp110_1, OFFSET xfft256Kp110_2
allfft	DP			OFFSET xfft256Kp110_3, OFFSET xfft256Kp110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			64, GAP2_10_1+1
allfft	DD			64, 1, 64, 1, 1, 0
allfft	DD	7635000, 393216, 0.0273,	846976
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft384Kp411_1AMD
allfft	DP			OFFSET xfft384Kp411_2AMD
allfft	DP			OFFSET xfft384Kp411_3AMD
allfft	DP			OFFSET xfft384Kp411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			48, GAP2_11_4+1
allfft	DD			48, 1, 48, 1, 1, 0
allfft	DD	7635000, 393216, 0.0273,	846976
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft384Kp_1AMD, OFFSET xfft384Kp_2AMD
allfft	DP			OFFSET xfft384Kp_3AMD, OFFSET xfft384Kp_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			48, GAP2_11_4+1
allfft	DD			48, 1, 48, 1, 1, 0
allfft	DD	7635000, 393216, 0.0273,	1036288
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET xfft384Kp410_1AMD
allfft	DP			OFFSET xfft384Kp410_2AMD
allfft	DP			OFFSET xfft384Kp410_3AMD
allfft	DP			OFFSET xfft384Kp410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			96, GAP2_10_4+1
allfft	DD			96, 1, 96, 1, 1, 0
allfft	DD	7635000, 393216, 0.0273,	1036288
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET xfft384Kp210_1AMD
allfft	DP			OFFSET xfft384Kp210_2AMD
allfft	DP			OFFSET xfft384Kp210_3AMD
allfft	DP			OFFSET xfft384Kp210_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			96, GAP2_10_2+1
allfft	DD			96, 1, 96, 1, 1, 0
	DD	7635000, 393216, 0.0273,	1036288
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET xfft384Kp110_1AMD
	DP			OFFSET xfft384Kp110_2AMD
	DP			OFFSET xfft384Kp110_3AMD
	DP			OFFSET xfft384Kp110_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			96, GAP2_10_1+1
	DD			96, 1, 96, 1, 1, 0
allfft	DD	7635000, 393216, 0.0273,	846976
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft384Kp411_1, OFFSET xfft384Kp411_2
allfft	DP			OFFSET xfft384Kp411_3, OFFSET xfft384Kp411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			48, GAP2_11_4+1
allfft	DD			48, 1, 48, 1, 1, 0
allfft	DD	7635000, 393216, 0.0273,	846976
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft384Kp_1, OFFSET xfft384Kp_2
allfft	DP			OFFSET xfft384Kp_3, OFFSET xfft384Kp_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			48, GAP2_11_4+1
allfft	DD			48, 1, 48, 1, 1, 0
	DD	7635000, 393216, 0.0273,	1036288
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			49152+11*128	;; scratch area size
	DP			OFFSET xfft384Kp410_1, OFFSET xfft384Kp410_2
	DP			OFFSET xfft384Kp410_3, OFFSET xfft384Kp410_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			96, GAP2_10_4+1
	DD			96, 1, 96, 1, 1, 0
	DD	7635000, 393216, 0.0273,	1036288
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			24576+11*128	;; scratch area size
	DP			OFFSET xfft384Kp210_1, OFFSET xfft384Kp210_2
	DP			OFFSET xfft384Kp210_3, OFFSET xfft384Kp210_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			96, GAP2_10_2+1
	DD			96, 1, 96, 1, 1, 0
allfft	DD	7635000, 393216, 0.0273,	1036288
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET xfft384Kp110_1, OFFSET xfft384Kp110_2
allfft	DP			OFFSET xfft384Kp110_3, OFFSET xfft384Kp110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			96, GAP2_10_1+1
allfft	DD			96, 1, 96, 1, 1, 0
allfft	DD	10150000, 524288, 0.0368,	1113088
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET xfft512Kp411_1AMD
allfft	DP			OFFSET xfft512Kp411_2AMD
allfft	DP			OFFSET xfft512Kp411_3AMD
allfft	DP			OFFSET xfft512Kp411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			64, GAP2_11_4+1
allfft	DD			64, 1, 64, 1, 1, 0
allfft	DD	10150000, 524288, 0.0368,	1113088
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft512Kp_1AMD, OFFSET xfft512Kp_2AMD
allfft	DP			OFFSET xfft512Kp_3AMD, OFFSET xfft512Kp_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			64, GAP2_11_4+1
allfft	DD			64, 1, 64, 1, 1, 0
allfft	DD	10150000, 524288, 0.0368,	1372160
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET xfft512Kp410_1AMD
allfft	DP			OFFSET xfft512Kp410_2AMD
allfft	DP			OFFSET xfft512Kp410_3AMD
allfft	DP			OFFSET xfft512Kp410_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			128, GAP2_10_4+1
allfft	DD			128, 1, 128, 1, 1, 0
allfft	DD	10150000, 524288, 0.0368,	1372160
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET xfft512Kp210_1AMD
allfft	DP			OFFSET xfft512Kp210_2AMD
allfft	DP			OFFSET xfft512Kp210_3AMD
allfft	DP			OFFSET xfft512Kp210_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			128, GAP2_10_2+1
allfft	DD			128, 1, 128, 1, 1, 0
	DD	10150000, 524288, 0.0368,	1372160
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET xfft512Kp110_1AMD
	DP			OFFSET xfft512Kp110_2AMD
	DP			OFFSET xfft512Kp110_3AMD
	DP			OFFSET xfft512Kp110_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			128, GAP2_10_1+1
	DD			128, 1, 128, 1, 1, 0
allfft	DD	10150000, 524288, 0.0368,	1113088
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET xfft512Kp411_1, OFFSET xfft512Kp411_2
allfft	DP			OFFSET xfft512Kp411_3, OFFSET xfft512Kp411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			64, GAP2_11_4+1
allfft	DD			64, 1, 64, 1, 1, 0
allfft	DD	10150000, 524288, 0.0368,	1113088
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			0		;; scratch area size
allfft	DP			OFFSET xfft512Kp_1, OFFSET xfft512Kp_2
allfft	DP			OFFSET xfft512Kp_3, OFFSET xfft512Kp_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			64, GAP2_11_4+1
allfft	DD			64, 1, 64, 1, 1, 0
	DD	10150000, 524288, 0.0368,	1372160
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET xfft512Kp410_1, OFFSET xfft512Kp410_2
	DP			OFFSET xfft512Kp410_3, OFFSET xfft512Kp410_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			128, GAP2_10_4+1
	DD			128, 1, 128, 1, 1, 0
	DD	10150000, 524288, 0.0368,	1372160
	DD			2*2		;; Flags, min_l2_cache, clm
	DD			32768+15*128	;; scratch area size
	DP			OFFSET xfft512Kp210_1, OFFSET xfft512Kp210_2
	DP			OFFSET xfft512Kp210_3, OFFSET xfft512Kp210_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			10		;; FFT levels done in pass2
	DD			128, GAP2_10_2+1
	DD			128, 1, 128, 1, 1, 0
allfft	DD	10150000, 524288, 0.0368,	1372160
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET xfft512Kp110_1, OFFSET xfft512Kp110_2
allfft	DP			OFFSET xfft512Kp110_3, OFFSET xfft512Kp110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			128, GAP2_10_1+1
allfft	DD			128, 1, 128, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	1285248
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp412_1AMD
allfft	DP			OFFSET xfft768Kp412_2AMD
allfft	DP			OFFSET xfft768Kp412_3AMD
allfft	DP			OFFSET xfft768Kp412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			48, GAP2_12_4+1
allfft	DD			48, 1, 48, 1, 1, 0
	DD	15040000, 786432, 0.0566,	1285248
	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET xfft768Kp212_1AMD
	DP			OFFSET xfft768Kp212_2AMD
	DP			OFFSET xfft768Kp212_3AMD
	DP			OFFSET xfft768Kp212_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			48, GAP2_12_2+1
	DD			48, 1, 48, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	1285248
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			6144+5*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp112_1AMD
allfft	DP			OFFSET xfft768Kp112_2AMD
allfft	DP			OFFSET xfft768Kp112_3AMD
allfft	DP			OFFSET xfft768Kp112_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			48, GAP2_12_1+1
allfft	DD			48, 1, 48, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	1648640
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp411_1AMD
allfft	DP			OFFSET xfft768Kp411_2AMD
allfft	DP			OFFSET xfft768Kp411_3AMD
allfft	DP			OFFSET xfft768Kp411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			96, GAP2_11_4+1
allfft	DD			96, 1, 96, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	1285248
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			24576+5*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp412_1, OFFSET xfft768Kp412_2
allfft	DP			OFFSET xfft768Kp412_3, OFFSET xfft768Kp412_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			48, GAP2_12_4+1
allfft	DD			48, 1, 48, 1, 1, 0
	DD	15040000, 786432, 0.0566,	1285248
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			12288+5*128	;; scratch area size
	DP			OFFSET xfft768Kp212_1, OFFSET xfft768Kp212_2
	DP			OFFSET xfft768Kp212_3, OFFSET xfft768Kp212_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			48, GAP2_12_2+1
	DD			48, 1, 48, 1, 1, 0
	DD	15040000, 786432, 0.0566,	1648640
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			49152+11*128	;; scratch area size
	DP			OFFSET xfft768Kp411_1, OFFSET xfft768Kp411_2
	DP			OFFSET xfft768Kp411_3, OFFSET xfft768Kp411_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			96, GAP2_11_4+1
	DD			96, 1, 96, 1, 1, 0
	DD	15040000, 786432, 0.0566,	1648640
	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
	DD			24576+11*128	;; scratch area size
	DP			OFFSET xfft768Kp211_1, OFFSET xfft768Kp211_2
	DP			OFFSET xfft768Kp211_3, OFFSET xfft768Kp211_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			96, GAP2_11_2+1
	DD			96, 1, 96, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	1648640
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp111_1, OFFSET xfft768Kp111_2
allfft	DP			OFFSET xfft768Kp111_3, OFFSET xfft768Kp111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			96, GAP2_11_1+1
allfft	DD			96, 1, 96, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	1648640
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp9111_1, OFFSET xfft768Kp9111_2
allfft	DP			OFFSET xfft768Kp9111_3, OFFSET xfft768Kp9111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			96, GAP2_11_1+1
allfft	DD			96, 1, 96, 1, 1, 0
	DD	15040000, 786432, 0.0566,	1648640
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET xfft768Kp011_1, OFFSET xfft768Kp011_2
	DP			OFFSET xfft768Kp011_3, OFFSET xfft768Kp011_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			96, GAP2_11_0+1
	DD			96, 1, 96, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	1648640
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			12288+11*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp9011_1, OFFSET xfft768Kp9011_2
allfft	DP			OFFSET xfft768Kp9011_3, OFFSET xfft768Kp9011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			96, GAP2_11_0+1
allfft	DD			96, 1, 96, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	2050176
allfft	DD			2*4		;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp410_1, OFFSET xfft768Kp410_2
allfft	DP			OFFSET xfft768Kp410_3, OFFSET xfft768Kp410_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			192, GAP2_10_4+1
allfft	DD			192, 1, 192, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	2050176
allfft	DD			2*2		;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp210_1, OFFSET xfft768Kp210_2
allfft	DP			OFFSET xfft768Kp210_3, OFFSET xfft768Kp210_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			192, GAP2_10_2+1
allfft	DD			192, 1, 192, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	2050176
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp110_1, OFFSET xfft768Kp110_2
allfft	DP			OFFSET xfft768Kp110_3, OFFSET xfft768Kp110_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			192, GAP2_10_1+1
allfft	DD			192, 1, 192, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	2050176
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp010_1, OFFSET xfft768Kp010_2
allfft	DP			OFFSET xfft768Kp010_3, OFFSET xfft768Kp010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			192, GAP2_10_0+1
allfft	DD			192, 1, 192, 1, 1, 0
allfft	DD	15040000, 786432, 0.0566,	3790464
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft768Kp18_1, OFFSET xfft768Kp18_2
allfft	DP			OFFSET xfft768Kp18_3, OFFSET xfft768Kp18_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			8		;; FFT levels done in pass2
allfft	DD			768, GAP2_8_1+1
allfft	DD			768, 1, 768, 1, 1, 0
allfft	DD	20000000, 1048576, 0.0761,	1682432
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			32768+7*128	;; scratch area size
allfft	DP			OFFSET xfft1024Kp412_1AMD
allfft	DP			OFFSET xfft1024Kp412_2AMD
allfft	DP			OFFSET xfft1024Kp412_3AMD
allfft	DP			OFFSET xfft1024Kp412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			64, GAP2_12_4+1
allfft	DD			64, 1, 64, 1, 1, 0
allfft	DD	20000000, 1048576, 0.0761,	1682432
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			16384+7*128	;; scratch area size
allfft	DP			OFFSET xfft1024Kp212_1AMD
allfft	DP			OFFSET xfft1024Kp212_2AMD
allfft	DP			OFFSET xfft1024Kp212_3AMD
allfft	DP			OFFSET xfft1024Kp212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			64, GAP2_12_2+1
allfft	DD			64, 1, 64, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	1682432
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			8192+7*128	;; scratch area size
	DP			OFFSET xfft1024Kp112_1AMD
	DP			OFFSET xfft1024Kp112_2AMD
	DP			OFFSET xfft1024Kp112_3AMD
	DP			OFFSET xfft1024Kp112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			64, GAP2_12_1+1
	DD			64, 1, 64, 1, 1, 0
allfft	DD	20000000, 1048576, 0.0761,	2181120
allfft	DD			RPFW+512*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET xfft1024Kp411_1AMD
allfft	DP			OFFSET xfft1024Kp411_2AMD
allfft	DP			OFFSET xfft1024Kp411_3AMD
allfft	DP			OFFSET xfft1024Kp411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			128, GAP2_11_4+1
allfft	DD			128, 1, 128, 1, 1, 0
allfft	DD	20000000, 1048576, 0.0761,	2181120
allfft	DD			RPFW+256*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET xfft1024Kp211_1AMD
allfft	DP			OFFSET xfft1024Kp211_2AMD
allfft	DP			OFFSET xfft1024Kp211_3AMD
allfft	DP			OFFSET xfft1024Kp211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			128, GAP2_11_2+1
allfft	DD			128, 1, 128, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	1682432
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			32768+7*128	;; scratch area size
	DP			OFFSET xfft1024Kp412_1, OFFSET xfft1024Kp412_2
	DP			OFFSET xfft1024Kp412_3, OFFSET xfft1024Kp412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			64, GAP2_12_4+1
	DD			64, 1, 64, 1, 1, 0
allfft	DD	20000000, 1048576, 0.0761,	1682432
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			16384+7*128	;; scratch area size
allfft	DP			OFFSET xfft1024Kp212_1, OFFSET xfft1024Kp212_2
allfft	DP			OFFSET xfft1024Kp212_3, OFFSET xfft1024Kp212_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			64, GAP2_12_2+1
allfft	DD			64, 1, 64, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	2181120
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET xfft1024Kp411_1, OFFSET xfft1024Kp411_2
	DP			OFFSET xfft1024Kp411_3, OFFSET xfft1024Kp411_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			128, GAP2_11_4+1
	DD			128, 1, 128, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	1682432
	DD			CELE_D*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			8192+7*128	;; scratch area size
	DP			OFFSET xfft1024Kp112_1, OFFSET xfft1024Kp112_2
	DP			OFFSET xfft1024Kp112_3, OFFSET xfft1024Kp112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			64, GAP2_12_1+1
	DD			64, 1, 64, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	2181120
	DD			WILLI*65536+2*2	;; Flags, min_l2_cache, clm
	DD			32768+15*128	;; scratch area size
	DP			OFFSET xfft1024Kp211_1, OFFSET xfft1024Kp211_2
	DP			OFFSET xfft1024Kp211_3, OFFSET xfft1024Kp211_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			128, GAP2_11_2+1
	DD			128, 1, 128, 1, 1, 0
allfft	DD	20000000, 1048576, 0.0761,	2181120
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET xfft1024Kp111_1, OFFSET xfft1024Kp111_2
allfft	DP			OFFSET xfft1024Kp111_3, OFFSET xfft1024Kp111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			128, GAP2_11_1+1
allfft	DD			128, 1, 128, 1, 1, 0
allfft	DD	20000000, 1048576, 0.0761,	2181120
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET xfft1024Kp9111_1
allfft	DP			OFFSET xfft1024Kp9111_2
allfft	DP			OFFSET xfft1024Kp9111_3
allfft	DP			OFFSET xfft1024Kp9111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			128, GAP2_11_1+1
allfft	DD			128, 1, 128, 1, 1, 0
	DD	20000000, 1048576, 0.0761,	2181120
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET xfft1024Kp011_1, OFFSET xfft1024Kp011_2
	DP			OFFSET xfft1024Kp011_3, OFFSET xfft1024Kp011_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			128, GAP2_11_0+1
	DD			128, 1, 128, 1, 1, 0
allfft	DD	20000000, 1048576, 0.0761,	2181120
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			16384+15*128	;; scratch area size
allfft	DP			OFFSET xfft1024Kp9011_1
allfft	DP			OFFSET xfft1024Kp9011_2
allfft	DP			OFFSET xfft1024Kp9011_3
allfft	DP			OFFSET xfft1024Kp9011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			128, GAP2_11_0+1
allfft	DD			128, 1, 128, 1, 1, 0
allfft	DD	29640000, 1572864, 0.125,	2480128
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET xfft1536Kp412_1AMD
allfft	DP			OFFSET xfft1536Kp412_2AMD
allfft	DP			OFFSET xfft1536Kp412_3AMD
allfft	DP			OFFSET xfft1536Kp412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			96, GAP2_12_4+1
allfft	DD			96, 1, 96, 1, 1, 0
allfft	DD	29640000, 1572864, 0.125,	2480128
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET xfft1536Kp212_1AMD
allfft	DP			OFFSET xfft1536Kp212_2AMD
allfft	DP			OFFSET xfft1536Kp212_3AMD
allfft	DP			OFFSET xfft1536Kp212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			96, GAP2_12_2+1
allfft	DD			96, 1, 96, 1, 1, 0
	DD	29640000, 1572864, 0.125,	2480128
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET xfft1536Kp112_1AMD
	DP			OFFSET xfft1536Kp112_2AMD
	DP			OFFSET xfft1536Kp112_3AMD
	DP			OFFSET xfft1536Kp112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			96, GAP2_12_1+1
	DD			96, 1, 96, 1, 1, 0
allfft	DD	29640000, 1572864, 0.125,	3252352
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536Kp411_1AMD
allfft	DP			OFFSET xfft1536Kp411_2AMD
allfft	DP			OFFSET xfft1536Kp411_3AMD
allfft	DP			OFFSET xfft1536Kp411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_4+1
allfft	DD			192, 1, 192, 1, 1, 0
allfft	DD	29640000, 1572864, 0.125,	3252352
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536Kp211_1AMD
allfft	DP			OFFSET xfft1536Kp211_2AMD
allfft	DP			OFFSET xfft1536Kp211_3AMD
allfft	DP			OFFSET xfft1536Kp211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_2+1
allfft	DD			192, 1, 192, 1, 1, 0
allfft	DD	29640000, 1572864, 0.125,	3252352
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536Kp111_1AMD
allfft	DP			OFFSET xfft1536Kp111_2AMD
allfft	DP			OFFSET xfft1536Kp111_3AMD
allfft	DP			OFFSET xfft1536Kp111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_1+1
allfft	DD			192, 1, 192, 1, 1, 0
allfft	DD	29640000, 1572864, 0.125,	3252352
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536Kp411_1, OFFSET xfft1536Kp411_2
allfft	DP			OFFSET xfft1536Kp411_3, OFFSET xfft1536Kp411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_4+1
allfft	DD			192, 1, 192, 1, 1, 0
	DD	29640000, 1572864, 0.125,	2480128
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			49152+23*128	;; scratch area size
	DP			OFFSET xfft1536Kp412_1, OFFSET xfft1536Kp412_2
	DP			OFFSET xfft1536Kp412_3, OFFSET xfft1536Kp412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			96, GAP2_12_4+1
	DD			96, 1, 96, 1, 1, 0
allfft	DD	29640000, 1572864, 0.125,	2480128
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			24576+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536Kp212_1, OFFSET xfft1536Kp212_2
allfft	DP			OFFSET xfft1536Kp212_3, OFFSET xfft1536Kp212_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			96, GAP2_12_2+1
allfft	DD			96, 1, 96, 1, 1, 0
	DD	29640000, 1572864, 0.125,	2480128
	DD			CELE_D*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			12288+23*128	;; scratch area size
	DP			OFFSET xfft1536Kp112_1, OFFSET xfft1536Kp112_2
	DP			OFFSET xfft1536Kp112_3, OFFSET xfft1536Kp112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			96, GAP2_12_1+1
	DD			96, 1, 96, 1, 1, 0
allfft	DD	29640000, 1572864, 0.125,	3252352
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET xfft1536Kp211_1, OFFSET xfft1536Kp211_2
allfft	DP			OFFSET xfft1536Kp211_3, OFFSET xfft1536Kp211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			192, GAP2_11_2+1
allfft	DD			192, 1, 192, 1, 1, 0
	DD	29640000, 1572864, 0.125,	3252352
	DD			WILLI*65536+2*1	;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET xfft1536Kp111_1, OFFSET xfft1536Kp111_2
	DP			OFFSET xfft1536Kp111_3, OFFSET xfft1536Kp111_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			192, GAP2_11_1+1
	DD			192, 1, 192, 1, 1, 0
	DD	29640000, 1572864, 0.125,	3252352
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET xfft1536Kp011_1, OFFSET xfft1536Kp011_2
	DP			OFFSET xfft1536Kp011_3, OFFSET xfft1536Kp011_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			192, GAP2_11_0+1
	DD			192, 1, 192, 1, 1, 0
allfft	DD	39390000, 2097152, 0.169,	3274752
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET xfft2048Kp412_1AMD
allfft	DP			OFFSET xfft2048Kp412_2AMD
allfft	DP			OFFSET xfft2048Kp412_3AMD
allfft	DP			OFFSET xfft2048Kp412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			128, GAP2_12_4+1
allfft	DD			128, 1, 128, 1, 1, 0
allfft	DD	39390000, 2097152, 0.169,	3274752
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET xfft2048Kp212_1AMD
allfft	DP			OFFSET xfft2048Kp212_2AMD
allfft	DP			OFFSET xfft2048Kp212_3AMD
allfft	DP			OFFSET xfft2048Kp212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			128, GAP2_12_2+1
allfft	DD			128, 1, 128, 1, 1, 0
	DD	39390000, 2097152, 0.169,	3274752
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET xfft2048Kp112_1AMD
	DP			OFFSET xfft2048Kp112_2AMD
	DP			OFFSET xfft2048Kp112_3AMD
	DP			OFFSET xfft2048Kp112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			128, GAP2_12_1+1
	DD			128, 1, 128, 1, 1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048Kp411_1AMD
allfft	DP			OFFSET xfft2048Kp411_2AMD
allfft	DP			OFFSET xfft2048Kp411_3AMD
allfft	DP			OFFSET xfft2048Kp411_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_4+1
allfft	DD			256, 1, 256, 1, 1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			RPFW+512*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048Kp211_1AMD
allfft	DP			OFFSET xfft2048Kp211_2AMD
allfft	DP			OFFSET xfft2048Kp211_3AMD
allfft	DP			OFFSET xfft2048Kp211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_2+1
allfft	DD			256, 1, 256, 1, 1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048Kp111_1AMD
allfft	DP			OFFSET xfft2048Kp111_2AMD
allfft	DP			OFFSET xfft2048Kp111_3AMD
allfft	DP			OFFSET xfft2048Kp111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_1+1
allfft	DD			256, 1, 256, 1, 1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048Kp011_1AMD
allfft	DP			OFFSET xfft2048Kp011_2AMD
allfft	DP			OFFSET xfft2048Kp011_3AMD
allfft	DP			OFFSET xfft2048Kp011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_0+1
allfft	DD			256, 1, 256, 1, 1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048Kp411_1, OFFSET xfft2048Kp411_2
allfft	DP			OFFSET xfft2048Kp411_3, OFFSET xfft2048Kp411_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_4+1
allfft	DD			256, 1, 256, 1, 1, 0
	DD	39390000, 2097152, 0.169,	3274752
	DD			512*65536+2*4	;; Flags, min_l2_cache, clm
	DD			65536+15*128	;; scratch area size
	DP			OFFSET xfft2048Kp412_1, OFFSET xfft2048Kp412_2
	DP			OFFSET xfft2048Kp412_3, OFFSET xfft2048Kp412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			128, GAP2_12_4+1
	DD			128, 1, 128, 1, 1, 0
allfft	DD	39390000, 2097152, 0.169,	3274752
allfft	DD			256*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET xfft2048Kp212_1, OFFSET xfft2048Kp212_2
allfft	DP			OFFSET xfft2048Kp212_3, OFFSET xfft2048Kp212_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			128, GAP2_12_2+1
allfft	DD			128, 1, 128, 1, 1, 0
	DD	39390000, 2097152, 0.169,	3274752
	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET xfft2048Kp112_1, OFFSET xfft2048Kp112_2
	DP			OFFSET xfft2048Kp112_3, OFFSET xfft2048Kp112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			128, GAP2_12_1+1
	DD			128, 1, 128, 1, 1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048Kp211_1, OFFSET xfft2048Kp211_2
allfft	DP			OFFSET xfft2048Kp211_3, OFFSET xfft2048Kp211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_2+1
allfft	DD			256, 1, 256, 1, 1, 0
allfft	DD	39390000, 2097152, 0.169,	4317184
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			32768+31*128	;; scratch area size
allfft	DP			OFFSET xfft2048Kp111_1, OFFSET xfft2048Kp111_2
allfft	DP			OFFSET xfft2048Kp111_3, OFFSET xfft2048Kp111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			256, GAP2_11_1+1
allfft	DD			256, 1, 256, 1, 1, 0
	DD	39390000, 2097152, 0.169,	4317184
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET xfft2048Kp011_1, OFFSET xfft2048Kp011_2
	DP			OFFSET xfft2048Kp011_3, OFFSET xfft2048Kp011_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			11		;; FFT levels done in pass2
	DD			256, GAP2_11_0+1
	DD			256, 1, 256, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	4143104
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp413_1AMD
allfft	DP			OFFSET xfft3072Kp413_2AMD
allfft	DP			OFFSET xfft3072Kp413_3AMD
allfft	DP			OFFSET xfft3072Kp413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			96, GAP2_13_4+1
allfft	DD			96, 1, 96, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	4143104
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			24576+11*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp213_1AMD
allfft	DP			OFFSET xfft3072Kp213_2AMD
allfft	DP			OFFSET xfft3072Kp213_3AMD
allfft	DP			OFFSET xfft3072Kp213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			96, GAP2_13_2+1
allfft	DD			96, 1, 96, 1, 1, 0
	DD	58410000, 3145728, 0.289,	4143104
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			12288+11*128	;; scratch area size
	DP			OFFSET xfft3072Kp113_1AMD
	DP			OFFSET xfft3072Kp113_2AMD
	DP			OFFSET xfft3072Kp113_3AMD
	DP			OFFSET xfft3072Kp113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			96, GAP2_13_1+1
	DD			96, 1, 96, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	4870272
allfft	DD			RPFW+2*4	;; Flags, min_l2_cache, clm
allfft	DD			98304+23*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp412_1AMD
allfft	DP			OFFSET xfft3072Kp412_2AMD
allfft	DP			OFFSET xfft3072Kp412_3AMD
allfft	DP			OFFSET xfft3072Kp412_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			192, GAP2_12_4+1
allfft	DD			192, 1, 192, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	4870272
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			49152+23*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp212_1AMD
allfft	DP			OFFSET xfft3072Kp212_2AMD
allfft	DP			OFFSET xfft3072Kp212_3AMD
allfft	DP			OFFSET xfft3072Kp212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			192, GAP2_12_2+1
allfft	DD			192, 1, 192, 1, 1, 0
	DD	58410000, 3145728, 0.289,	4870272
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET xfft3072Kp112_1AMD
	DP			OFFSET xfft3072Kp112_2AMD
	DP			OFFSET xfft3072Kp112_3AMD
	DP			OFFSET xfft3072Kp112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			192, GAP2_12_1+1
	DD			192, 1, 192, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp211_1AMD
allfft	DP			OFFSET xfft3072Kp211_2AMD
allfft	DP			OFFSET xfft3072Kp211_3AMD
allfft	DP			OFFSET xfft3072Kp211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_2+1
allfft	DD			384, 1, 384, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			RPFW+256*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp111_1AMD
allfft	DP			OFFSET xfft3072Kp111_2AMD
allfft	DP			OFFSET xfft3072Kp111_3AMD
allfft	DP			OFFSET xfft3072Kp111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_1+1
allfft	DD			384, 1, 384, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp011_1AMD
allfft	DP			OFFSET xfft3072Kp011_2AMD
allfft	DP			OFFSET xfft3072Kp011_3AMD
allfft	DP			OFFSET xfft3072Kp011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_0+1
allfft	DD			384, 1, 384, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	4143104
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			49152+11*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp413_1, OFFSET xfft3072Kp413_2
allfft	DP			OFFSET xfft3072Kp413_3, OFFSET xfft3072Kp413_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			96, GAP2_13_4+1
allfft	DD			96, 1, 96, 1, 1, 0
	DD	58410000, 3145728, 0.289,	4870272
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			98304+23*128	;; scratch area size
	DP			OFFSET xfft3072Kp412_1, OFFSET xfft3072Kp412_2
	DP			OFFSET xfft3072Kp412_3, OFFSET xfft3072Kp412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			192, GAP2_12_4+1
	DD			192, 1, 192, 1, 1, 0
	DD	58410000, 3145728, 0.289,	4870272
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			49152+23*128	;; scratch area size
	DP			OFFSET xfft3072Kp212_1, OFFSET xfft3072Kp212_2
	DP			OFFSET xfft3072Kp212_3, OFFSET xfft3072Kp212_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			192, GAP2_12_2+1
	DD			192, 1, 192, 1, 1, 0
	DD	58410000, 3145728, 0.289,	4870272
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			24576+23*128	;; scratch area size
	DP			OFFSET xfft3072Kp112_1, OFFSET xfft3072Kp112_2
	DP			OFFSET xfft3072Kp112_3, OFFSET xfft3072Kp112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			192, GAP2_12_1+1
	DD			192, 1, 192, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp211_1, OFFSET xfft3072Kp211_2
allfft	DP			OFFSET xfft3072Kp211_3, OFFSET xfft3072Kp211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_2+1
allfft	DD			384, 1, 384, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			256*65536+2*1	;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp111_1, OFFSET xfft3072Kp111_2
allfft	DP			OFFSET xfft3072Kp111_3, OFFSET xfft3072Kp111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_1+1
allfft	DD			384, 1, 384, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	6459392
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			49152+47*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp011_1, OFFSET xfft3072Kp011_2
allfft	DP			OFFSET xfft3072Kp011_3, OFFSET xfft3072Kp011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			384, GAP2_11_0+1
allfft	DD			384, 1, 384, 1, 1, 0
allfft	DD	58410000, 3145728, 0.289,	8132736
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft3072Kp010_1, OFFSET xfft3072Kp010_2
allfft	DP			OFFSET xfft3072Kp010_3, OFFSET xfft3072Kp010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			768, GAP2_10_0+1
allfft	DD			768, 1, 768, 1, 1, 0
allfft	DD	77700000, 4194304, 0.425,	5462016
allfft	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
allfft	DD			65536+15*128	;; scratch area size
allfft	DP			OFFSET xfft4096Kp413_1AMD
allfft	DP			OFFSET xfft4096Kp413_2AMD
allfft	DP			OFFSET xfft4096Kp413_3AMD
allfft	DP			OFFSET xfft4096Kp413_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			128, GAP2_13_4+1
allfft	DD			128, 1, 128, 1, 1, 0
allfft	DD	77700000, 4194304, 0.425,	5462016
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			32768+15*128	;; scratch area size
allfft	DP			OFFSET xfft4096Kp213_1AMD
allfft	DP			OFFSET xfft4096Kp213_2AMD
allfft	DP			OFFSET xfft4096Kp213_3AMD
allfft	DP			OFFSET xfft4096Kp213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			128, GAP2_13_2+1
allfft	DD			128, 1, 128, 1, 1, 0
	DD	77700000, 4194304, 0.425,	5462016
	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			16384+15*128	;; scratch area size
	DP			OFFSET xfft4096Kp113_1AMD
	DP			OFFSET xfft4096Kp113_2AMD
	DP			OFFSET xfft4096Kp113_3AMD
	DP			OFFSET xfft4096Kp113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			128, GAP2_13_1+1
	DD			128, 1, 128, 1, 1, 0
allfft	DD	77700000, 4194304, 0.425,	6459392
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			65536+31*128	;; scratch area size
allfft	DP			OFFSET xfft4096Kp212_1AMD
allfft	DP			OFFSET xfft4096Kp212_2AMD
allfft	DP			OFFSET xfft4096Kp212_3AMD
allfft	DP			OFFSET xfft4096Kp212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			256, GAP2_12_2+1
allfft	DD			256, 1, 256, 1, 1, 0
	DD	77700000, 4194304, 0.425,	6459392
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET xfft4096Kp112_1AMD
	DP			OFFSET xfft4096Kp112_2AMD
	DP			OFFSET xfft4096Kp112_3AMD
	DP			OFFSET xfft4096Kp112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, GAP2_12_1+1
	DD			256, 1, 256, 1, 1, 0
allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096Kp211_1AMD
allfft	DP			OFFSET xfft4096Kp211_2AMD
allfft	DP			OFFSET xfft4096Kp211_3AMD
allfft	DP			OFFSET xfft4096Kp211_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_2+1
allfft	DD			512, 1, 512, 1, 1, 0
allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			RPFW+512*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096Kp111_1AMD
allfft	DP			OFFSET xfft4096Kp111_2AMD
allfft	DP			OFFSET xfft4096Kp111_3AMD
allfft	DP			OFFSET xfft4096Kp111_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_1+1
allfft	DD			512, 1, 512, 1, 1, 0
allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096Kp011_1AMD
allfft	DP			OFFSET xfft4096Kp011_2AMD
allfft	DP			OFFSET xfft4096Kp011_3AMD
allfft	DP			OFFSET xfft4096Kp011_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_0+1
allfft	DD			512, 1, 512, 1, 1, 0
	DD	77700000, 4194304, 0.425,	6459392
	DD			1024*65536+2*4	;; Flags, min_l2_cache, clm
	DD			131072+31*128	;; scratch area size
	DP			OFFSET xfft4096Kp412_1, OFFSET xfft4096Kp412_2
	DP			OFFSET xfft4096Kp412_3, OFFSET xfft4096Kp412_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, GAP2_12_4+1
	DD			256, 1, 256, 1, 1, 0
	DD	77700000, 4194304, 0.425,	6459392
	DD			512*65536+2*2	;; Flags, min_l2_cache, clm
	DD			65536+31*128	;; scratch area size
	DP			OFFSET xfft4096Kp212_1, OFFSET xfft4096Kp212_2
	DP			OFFSET xfft4096Kp212_3, OFFSET xfft4096Kp212_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, GAP2_12_2+1
	DD			256, 1, 256, 1, 1, 0
	DD	77700000, 4194304, 0.425,	6459392
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			32768+31*128	;; scratch area size
	DP			OFFSET xfft4096Kp112_1, OFFSET xfft4096Kp112_2
	DP			OFFSET xfft4096Kp112_3, OFFSET xfft4096Kp112_4
	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			256, GAP2_12_1+1
	DD			256, 1, 256, 1, 1, 0
allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096Kp211_1, OFFSET xfft4096Kp211_2
allfft	DP			OFFSET xfft4096Kp211_3, OFFSET xfft4096Kp211_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_2+1
allfft	DD			512, 1, 512, 1, 1, 0
allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096Kp111_1, OFFSET xfft4096Kp111_2
allfft	DP			OFFSET xfft4096Kp111_3, OFFSET xfft4096Kp111_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_1+1
allfft	DD			512, 1, 512, 1, 1, 0
allfft	DD	77700000, 4194304, 0.425,	8589312
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET xfft4096Kp011_1, OFFSET xfft4096Kp011_2
allfft	DP			OFFSET xfft4096Kp011_3, OFFSET xfft4096Kp011_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			11		;; FFT levels done in pass2
allfft	DD			512, GAP2_11_0+1
allfft	DD			512, 1, 512, 1, 1, 0
allfft	DD	77700000, 4194304, 0.425,	4*6459392
allfft	DD			2*1	;; Flags, min_l2_cache, clm
allfft	DD			131072+127*128	;; scratch area size
allfft	DP			OFFSET xfft4096Kp010_1, OFFSET xfft4096Kp010_2
allfft	DP			OFFSET xfft4096Kp010_3, OFFSET xfft4096Kp010_4
allfft	DP			OFFSET xprctab2	;; Table of add/sub/norm procs
allfft	DD			10		;; FFT levels done in pass2
allfft	DD			1024, GAP2_10_0+1
allfft	DD			1024, 1, 1024, 1, 1, 0
	DD	115000000, 6291456, 0.668,	9650176
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET xfft6Mp412_1AMD, OFFSET xfft6Mp412_2AMD
	DP			OFFSET xfft6Mp412_3AMD, OFFSET xfft6Mp412_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, GAP2_12_4+1
	DD			384, 1, 384, 1, 1, 0
allfft	DD	115000000, 6291456, 0.668,	9650176
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET xfft6Mp212_1AMD, OFFSET xfft6Mp212_2AMD
allfft	DP			OFFSET xfft6Mp212_3AMD, OFFSET xfft6Mp212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			384, GAP2_12_2+1
allfft	DD			384, 1, 384, 1, 1, 0
	DD	115000000, 6291456, 0.668,	9650176
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET xfft6Mp112_1AMD, OFFSET xfft6Mp112_2AMD
	DP			OFFSET xfft6Mp112_3AMD, OFFSET xfft6Mp112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, GAP2_12_1+1
	DD			384, 1, 384, 1, 1, 0
	DD	115000000, 6291456, 0.668,	9650176
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET xfft6Mp412_1, OFFSET xfft6Mp412_2
	DP			OFFSET xfft6Mp412_3, OFFSET xfft6Mp412_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, GAP2_12_4+1
	DD			384, 1, 384, 1, 1, 0
	DD	115000000, 6291456, 0.668,	9650176
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			98304+47*128	;; scratch area size
	DP			OFFSET xfft6Mp212_1, OFFSET xfft6Mp212_2
	DP			OFFSET xfft6Mp212_3, OFFSET xfft6Mp212_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, GAP2_12_2+1
	DD			384, 1, 384, 1, 1, 0
	DD	115000000, 6291456, 0.668,	9650176
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET xfft6Mp112_1, OFFSET xfft6Mp112_2
	DP			OFFSET xfft6Mp112_3, OFFSET xfft6Mp112_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			384, GAP2_12_1+1
	DD			384, 1, 384, 1, 1, 0
	DD	153100000, 8388608, 1.042,	12828672
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET xfft8Mp412_1AMD, OFFSET xfft8Mp412_2AMD
	DP			OFFSET xfft8Mp412_3AMD, OFFSET xfft8Mp412_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, GAP2_12_4+1
	DD			512, 1, 512, 1, 1, 0
allfft	DD	153100000, 8388608, 1.042,	12828672
allfft	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET xfft8Mp212_1AMD, OFFSET xfft8Mp212_2AMD
allfft	DP			OFFSET xfft8Mp212_3AMD, OFFSET xfft8Mp212_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			12		;; FFT levels done in pass2
allfft	DD			512, GAP2_12_2+1
allfft	DD			512, 1, 512, 1, 1, 0
	DD	153100000, 8388608, 1.042,	12828672
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET xfft8Mp112_1AMD, OFFSET xfft8Mp112_2AMD
	DP			OFFSET xfft8Mp112_3AMD, OFFSET xfft8Mp112_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, GAP2_12_1+1
	DD			512, 1, 512, 1, 1, 0
allfft	DD	153100000, 8388608, 1.042,	3*12828672
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			131072+31*128	;; scratch area size
allfft	DP			OFFSET xfft8Mp413_1, OFFSET xfft8Mp413_2
allfft	DP			OFFSET xfft8Mp413_3, OFFSET xfft8Mp413_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			256, GAP2_13_4+1
allfft	DD			256, 1, 256, 1, 1, 0
	DD	153100000, 8388608, 1.042,	12828672
	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET xfft8Mp412_1, OFFSET xfft8Mp412_2
	DP			OFFSET xfft8Mp412_3, OFFSET xfft8Mp412_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, GAP2_12_4+1
	DD			512, 1, 512, 1, 1, 0
	DD	153100000, 8388608, 1.042,	12828672
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			131072+63*128	;; scratch area size
	DP			OFFSET xfft8Mp212_1, OFFSET xfft8Mp212_2
	DP			OFFSET xfft8Mp212_3, OFFSET xfft8Mp212_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, GAP2_12_2+1
	DD			512, 1, 512, 1, 1, 0
	DD	153100000, 8388608, 1.042,	12828672
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET xfft8Mp112_1, OFFSET xfft8Mp112_2
	DP			OFFSET xfft8Mp112_3, OFFSET xfft8Mp112_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			12		;; FFT levels done in pass2
	DD			512, GAP2_12_1+1
	DD			512, 1, 512, 1, 1, 0
	DD	226800000, 12582912, 1.400,	16031744
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			196608+47*128	;; scratch area size
	DP			OFFSET xfft12Mp413_1AMD
	DP			OFFSET xfft12Mp413_2AMD
	DP			OFFSET xfft12Mp413_3AMD
	DP			OFFSET xfft12Mp413_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			384, GAP2_13_4+1
	DD			384, 1, 384, 1, 1, 0
allfft	DD	226800000, 12582912, 1.400,	16031744
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			98304+47*128	;; scratch area size
allfft	DP			OFFSET xfft12Mp213_1AMD
allfft	DP			OFFSET xfft12Mp213_2AMD
allfft	DP			OFFSET xfft12Mp213_3AMD
allfft	DP			OFFSET xfft12Mp213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			384, GAP2_13_2+1
allfft	DD			384, 1, 384, 1, 1, 0
	DD	226800000, 12582912, 1.400,	16031744
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET xfft12Mp113_1AMD
	DP			OFFSET xfft12Mp113_2AMD
	DP			OFFSET xfft12Mp113_3AMD
	DP			OFFSET xfft12Mp113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			384, GAP2_13_1+1
	DD			384, 1, 384, 1, 1, 0
allfft	DD	226800000, 12582912, 1.400,	16031744
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			196608+47*128	;; scratch area size
allfft	DP			OFFSET xfft12Mp413_1, OFFSET xfft12Mp413_2
allfft	DP			OFFSET xfft12Mp413_3, OFFSET xfft12Mp413_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			384, GAP2_13_4+1
allfft	DD			384, 1, 384, 1, 1, 0
	DD	226800000, 12582912, 1.400,	16031744
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			98304+47*128	;; scratch area size
	DP			OFFSET xfft12Mp213_1, OFFSET xfft12Mp213_2
	DP			OFFSET xfft12Mp213_3, OFFSET xfft12Mp213_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			384, GAP2_13_2+1
	DD			384, 1, 384, 1, 1, 0
	DD	226800000, 12582912, 1.400,	16031744
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			49152+47*128	;; scratch area size
	DP			OFFSET xfft12Mp113_1, OFFSET xfft12Mp113_2
	DP			OFFSET xfft12Mp113_3, OFFSET xfft12Mp113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			384, GAP2_13_1+1
	DD			384, 1, 384, 1, 1, 0
	DD	301900000, 16777216, 2.100,	21307392
	DD			RPFW+1024*65536+2*4 ;; Flags, min_l2_cache, clm
	DD			262144+63*128	;; scratch area size
	DP			OFFSET xfft16Mp413_1AMD
	DP			OFFSET xfft16Mp413_2AMD
	DP			OFFSET xfft16Mp413_3AMD
	DP			OFFSET xfft16Mp413_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			512, GAP2_13_4+1
	DD			512, 1, 512, 1, 1, 0
allfft	DD	301900000, 16777216, 2.100,	21307392
allfft	DD			RPFW+2*2	;; Flags, min_l2_cache, clm
allfft	DD			131072+63*128	;; scratch area size
allfft	DP			OFFSET xfft16Mp213_1AMD
allfft	DP			OFFSET xfft16Mp213_2AMD
allfft	DP			OFFSET xfft16Mp213_3AMD
allfft	DP			OFFSET xfft16Mp213_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			512, GAP2_13_2+1
allfft	DD			512, 1, 512, 1, 1, 0
	DD	301900000, 16777216, 2.100,	21307392
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET xfft16Mp113_1AMD
	DP			OFFSET xfft16Mp113_2AMD
	DP			OFFSET xfft16Mp113_3AMD
	DP			OFFSET xfft16Mp113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			512, GAP2_13_1+1
	DD			512, 1, 512, 1, 1, 0
allfft	DD	301900000, 16777216, 2.100,	21307392
allfft	DD			2048*65536+2*4	;; Flags, min_l2_cache, clm
allfft	DD			262144+63*128	;; scratch area size
allfft	DP			OFFSET xfft16Mp413_1, OFFSET xfft16Mp413_2
allfft	DP			OFFSET xfft16Mp413_3, OFFSET xfft16Mp413_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			512, GAP2_13_4+1
allfft	DD			512, 1, 512, 1, 1, 0
	DD	301900000, 16777216, 2.100,	21307392
	DD			1024*65536+2*2	;; Flags, min_l2_cache, clm
	DD			131072+63*128	;; scratch area size
	DP			OFFSET xfft16Mp213_1, OFFSET xfft16Mp213_2
	DP			OFFSET xfft16Mp213_3, OFFSET xfft16Mp213_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			512, GAP2_13_2+1
	DD			512, 1, 512, 1, 1, 0
	DD	301900000, 16777216, 2.100,	21307392
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			65536+63*128	;; scratch area size
	DP			OFFSET xfft16Mp113_1, OFFSET xfft16Mp113_2
	DP			OFFSET xfft16Mp113_3, OFFSET xfft16Mp113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			512, GAP2_13_1+1
	DD			512, 1, 512, 1, 1, 0
allfft	DD	301900000, 16777216, 2.100,	21307392
allfft	DD			2*1		;; Flags, min_l2_cache, clm
allfft	DD			65536+63*128	;; scratch area size
allfft	DP			OFFSET xfft16Mp013_1, OFFSET xfft16Mp013_2
allfft	DP			OFFSET xfft16Mp013_3, OFFSET xfft16Mp013_4
allfft	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			512, GAP2_13_0+1
allfft	DD			512, 1, 512, 1, 1, 0
	DD	446900000, 25165824, 3.300,	31883392
	DD			RPFW+1024*65536+2*2 ;; Flags, min_l2_cache, clm
	DD			196608+95*128	;; scratch area size
	DP			OFFSET xfft24Mp213_1AMD
	DP			OFFSET xfft24Mp213_2AMD
	DP			OFFSET xfft24Mp213_3AMD
	DP			OFFSET xfft24Mp213_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			768, GAP2_13_2+1
	DD			768, 1, 768, 1, 1, 0
allfft	DD	446900000, 25165824, 3.300,	31883392
allfft	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
allfft	DD			98304+95*128	;; scratch area size
allfft	DP			OFFSET xfft24Mp113_1AMD
allfft	DP			OFFSET xfft24Mp113_2AMD
allfft	DP			OFFSET xfft24Mp113_3AMD
allfft	DP			OFFSET xfft24Mp113_4AMD
allfft	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
allfft	DD			13		;; FFT levels done in pass2
allfft	DD			768, GAP2_13_1+1
allfft	DD			768, 1, 768, 1, 1, 0
	DD	446900000, 25165824, 3.300,	31883392
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET xfft24Mp013_1AMD
	DP			OFFSET xfft24Mp013_2AMD
	DP			OFFSET xfft24Mp013_3AMD
	DP			OFFSET xfft24Mp013_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			768, GAP2_13_0+1
	DD			768, 1, 768, 1, 1, 0
	DD	446900000, 25165824, 3.300,	31883392
	DD			2048*65536+2*2	;; Flags, min_l2_cache, clm
	DD			196608+95*128	;; scratch area size
	DP			OFFSET xfft24Mp213_1, OFFSET xfft24Mp213_2
	DP			OFFSET xfft24Mp213_3, OFFSET xfft24Mp213_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			768, GAP2_13_2+1
	DD			768, 1, 768, 1, 1, 0
	DD	446900000, 25165824, 3.300,	31883392
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET xfft24Mp113_1, OFFSET xfft24Mp113_2
	DP			OFFSET xfft24Mp113_3, OFFSET xfft24Mp113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			768, GAP2_13_1+1
	DD			768, 1, 768, 1, 1, 0
	DD	446900000, 25165824, 3.300,	31883392
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			98304+95*128	;; scratch area size
	DP			OFFSET xfft24Mp013_1, OFFSET xfft24Mp013_2
	DP			OFFSET xfft24Mp013_3, OFFSET xfft24Mp013_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			768, GAP2_13_0+1
	DD			768, 1, 768, 1, 1, 0
	DD	594600000, 33554432, 4.500,	42434560
	DD			RPFW+1024*65536+2*1 ;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET xfft32Mp113_1AMD
	DP			OFFSET xfft32Mp113_2AMD
	DP			OFFSET xfft32Mp113_3AMD
	DP			OFFSET xfft32Mp113_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			1024, GAP2_13_1+1
	DD			1024, 1, 1024, 1, 1, 0
	DD	594600000, 33554432, 4.500,	42434560
	DD			RPFW+2*1	;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET xfft32Mp013_1AMD
	DP			OFFSET xfft32Mp013_2AMD
	DP			OFFSET xfft32Mp013_3AMD
	DP			OFFSET xfft32Mp013_4AMD
	DP			OFFSET xprctab2a ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			1024, GAP2_13_0+1
	DD			1024, 1, 1024, 1, 1, 0
	DD	594600000, 33554432, 4.500,	42434560
	DD			2048*65536+2*2	;; Flags, min_l2_cache, clm
	DD			262144+127*128	;; scratch area size
	DP			OFFSET xfft32Mp213_1, OFFSET xfft32Mp213_2
	DP			OFFSET xfft32Mp213_3, OFFSET xfft32Mp213_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			1024, GAP2_13_2+1
	DD			1024, 1, 1024, 1, 1, 0
	DD	594600000, 33554432, 4.500,	42434560
	DD			1024*65536+2*1	;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET xfft32Mp113_1, OFFSET xfft32Mp113_2
	DP			OFFSET xfft32Mp113_3, OFFSET xfft32Mp113_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			1024, GAP2_13_1+1
	DD			1024, 1, 1024, 1, 1, 0
	DD	594600000, 33554432, 4.500,	42434560
	DD			2*1		;; Flags, min_l2_cache, clm
	DD			131072+127*128	;; scratch area size
	DP			OFFSET xfft32Mp013_1, OFFSET xfft32Mp013_2
	DP			OFFSET xfft32Mp013_3, OFFSET xfft32Mp013_4
	DP			OFFSET xprctab2 ;; Table of add/sub/norm procs
	DD			13		;; FFT levels done in pass2
	DD			1024, GAP2_13_0+1
	DD			1024, 1, 1024, 1, 1, 0
	DD	0

	;; Align so that other GWDATA areas are also aligned on a cache line
	align 128
_GWDATA ENDS

END

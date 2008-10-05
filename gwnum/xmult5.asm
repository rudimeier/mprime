; Copyright 2005-2007 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; This code does two passes, 13 levels on the second pass.
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
INCLUDE xfft5.mac
INCLUDE	xlucas.mac
INCLUDE xmult.mac
INCLUDE xpass1.mac
INCLUDE xpass1sc.mac

EXTRN	pass1_aux_entry_point_return:PROC
EXTRN	xgw_finish_fft:PROC
EXTRN	xgw_carries:PROC
EXTRN	xgw_finish_mult:PROC

EXTRNP	xpass2_13_levels
EXTRNP	xpass2_13_levels_np

_TEXT SEGMENT

;; Pass 2 does 13 FFT levels on 2 sets of data (2 * 2^13 complex values =
;; 2^15 doubles = 256KB).

PREFETCHING = 1
blkdst = (32*(8192+128)+GAP2_13_4)
xpass2_levels = 13

;; All the FFT routines for each FFT length

	EXPANDING = 2
;	xfft	512K
;	xfft	512Kp
;	xfft	640K
;	xfft	768K
;	xfft	768Kp
;	xfft	896K
;	xfft	1024K
;	xfft	1024Kp

IFDEF AMD
allfft	xfftclm	1280K, 4
allfft	xfftclm	1536K, 4
allfft	xfftclm	1536Kp, 4
allfft	xfftclm	1792K, 4
allfft	xfftclm	2048K, 4
allfft	xfftclm	2048Kp, 4
allfft	xfftclm	2560K, 4
allfft	xfftclm	2560K, 2
	xfftclm	2560K, 1
allfft	xfftclm	3072K, 4
allfft	xfftclm	3072K, 2
	xfftclm	3072K, 1
allfft	xfftclm	3072Kp, 4
allfft	xfftclm	3072Kp, 2
	xfftclm	3072Kp, 1
allfft	xfftclm	3584K, 4
allfft	xfftclm	3584K, 2
	xfftclm	3584K, 1
allfft	xfftclm	4096K, 4
allfft	xfftclm	4096K, 2
	xfftclm	4096K, 1
allfft	xfftclm	4096Kp, 4
allfft	xfftclm	4096Kp, 2
	xfftclm	4096Kp, 1
allfft	xfftclm	5M, 4
allfft	xfftclm	5M, 2
allfft	xfftclm	5M, 1
allfft	xfftclm	6M, 4
allfft	xfftclm	6M, 2
allfft	xfftclm	6M, 1
allfft	xfftclm	6Mp, 4
allfft	xfftclm	6Mp, 2
allfft	xfftclm	6Mp, 1
allfft	xfftclm	7M, 4
allfft	xfftclm	7M, 2
allfft	xfftclm	7M, 1
allfft	xfftclm	8M, 4
allfft	xfftclm	8M, 2
allfft	xfftclm	8M, 1
allfft	xfftclm	8Mp, 4
allfft	xfftclm	8Mp, 2
allfft	xfftclm	8Mp, 1
	xfftclm	10M, 4
allfft	xfftclm	10M, 2
	xfftclm	10M, 1
	xfftclm	12M, 4
allfft	xfftclm	12M, 2
	xfftclm	12M, 1
	xfftclm	12Mp, 4
allfft	xfftclm	12Mp, 2
	xfftclm	12Mp, 1
	xfftclm	14M, 4
allfft	xfftclm	14M, 2
	xfftclm	14M, 1
	xfftclm	16M, 4
allfft	xfftclm	16M, 2
	xfftclm	16M, 1
	xfftclm	16Mp, 4
allfft	xfftclm	16Mp, 2
	xfftclm	16Mp, 1
allfft	xfftclm	20M, 4
	xfftclm	20M, 2
allfft	xfftclm	20M, 1
	xfftclm	20M, 0
	xfftclm	24M, 2
allfft	xfftclm	24M, 1
	xfftclm	24M, 0
	xfftclm	24Mp, 2
allfft	xfftclm	24Mp, 1
	xfftclm	24Mp, 0
	xfftclm	28M, 2
allfft	xfftclm	28M, 1
	xfftclm	28M, 0
allfft	xfftclm	32M, 2
	xfftclm	32M, 1
	xfftclm	32M, 0
allfft	xfftclm	32Mp, 2
	xfftclm	32Mp, 1
	xfftclm	32Mp, 0
ELSE
allfft	xfftclm	1280K, 4
allfft	xfftclm	1536K, 4
allfft	xfftclm	1536Kp, 4
allfft	xfftclm	1792K, 4
allfft	xfftclm	2048K, 4
allfft	xfftclm	2048Kp, 4
allfft	xfftclm	2560K, 4
allfft	xfftclm	2560K, 2
allfft	xfftclm	2560K, 1
allfft	xfftclm	3072K, 4
allfft	xfftclm	3072K, 2
allfft	xfftclm	3072K, 1
allfft	xfftclm	3072Kp, 4
allfft	xfftclm	3584K, 4
allfft	xfftclm	3584K, 2
allfft	xfftclm	3584K, 1
allfft	xfftclm	4096K, 4
allfft	xfftclm	4096K, 2
allfft	xfftclm	4096K, 1
allfft	xfftclm	4096Kp, 4
allfft	xfftclm	4096Kp, 2
allfft	xfftclm	5M, 4
allfft	xfftclm	5M, 2
allfft	xfftclm	5M, 1
allfft	xfftclm	6M, 4
allfft	xfftclm	6M, 2
allfft	xfftclm	6M, 1
allfft	xfftclm	6Mp, 4
allfft	xfftclm	6Mp, 2
allfft	xfftclm	6Mp, 1
allfft	xfftclm	7M, 4
allfft	xfftclm	7M, 2
allfft	xfftclm	7M, 1
allfft	xfftclm	8M, 4
allfft	xfftclm	8M, 2
allfft	xfftclm	8M, 1
allfft	xfftclm	8Mp, 4
allfft	xfftclm	8Mp, 2
allfft	xfftclm	8Mp, 1
	xfftclm	10M, 4
allfft	xfftclm	10M, 2
	xfftclm	10M, 1
	xfftclm	12M, 4
	xfftclm	12M, 2
	xfftclm	12M, 1
allfft	xfftclm	12Mp, 4
	xfftclm	12Mp, 2
	xfftclm	12Mp, 1
	xfftclm	14M, 4
	xfftclm	14M, 2
	xfftclm	14M, 1
	xfftclm	16M, 4
	xfftclm	16M, 2
	xfftclm	16M, 1
allfft	xfftclm	16Mp, 4
	xfftclm	16Mp, 2
	xfftclm	16Mp, 1
allfft	xfftclm	16Mp, 0
	xfftclm	20M, 2
	xfftclm	20M, 1
allfft	xfftclm	20M, 0
	xfftclm	24M, 2
	xfftclm	24M, 1
	xfftclm	24M, 0
	xfftclm	24Mp, 2
	xfftclm	24Mp, 1
	xfftclm	24Mp, 0
	xfftclm	28M, 2
	xfftclm	28M, 1
	xfftclm	28M, 0
	xfftclm	32M, 2
	xfftclm	32M, 1
	xfftclm	32M, 0
	xfftclm	32Mp, 2
	xfftclm	32Mp, 1
	xfftclm	32Mp, 0
ENDIF

_TEXT	ENDS
END

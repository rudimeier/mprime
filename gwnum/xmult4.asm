; Copyright 2005-2007 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; This code does two passes, 12 levels on the second pass.
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
INCLUDE xfft4.mac
INCLUDE	xlucas.mac
INCLUDE xmult.mac
INCLUDE xpass1.mac
INCLUDE xpass1sc.mac

EXTRN	pass1_aux_entry_point_return:PROC
EXTRN	xgw_finish_fft:PROC
EXTRN	xgw_carries:PROC
EXTRN	xgw_finish_mult:PROC

EXTRNP	xpass2_12_levels
EXTRNP	xpass2_12_levels_np

_TEXT SEGMENT

;; Pass 2 does 12 FFT levels 2 sets of data (2 * 2^12 complex values =
;; 2^14 doubles = 128KB).

PREFETCHING = 1
blkdst = (16*(8192+128)+GAP2_12_4)
xpass2_levels = 12

;; All the FFT routines for each FFT length

	EXPANDING = 2
;	xfft	256K
;	xfft	256Kp
;	xfft	320K
;	xfft	384K
;	xfft	384Kp
;	xfft	448K
;	xfft	512K
;	xfft	512Kp

IFDEF AMD
allfft	xfftclm	640K, 8
allfft	xfftclm	640K, 4
	xfftclm	640K, 2
allfft	xfftclm	640K, 1
allfft	xfftclm	768K, 4
	xfftclm	768K, 2
allfft	xfftclm	768K, 1
allfft	xfftclm	768Kp, 4
	xfftclm	768Kp, 2
allfft	xfftclm	768Kp, 1
allfft	xfftclm	896K, 4
	xfftclm	896K, 2
	xfftclm	896K, 1
allfft	xfftclm	1024K, 4
	xfftclm	1024K, 2
	xfftclm	1024K, 1
allfft	xfftclm	1024Kp, 4
allfft	xfftclm	1024Kp, 2
	xfftclm	1024Kp, 1
allfft	xfftclm	1280K, 4
allfft	xfftclm	1280K, 2
	xfftclm	1280K, 1
allfft	xfftclm	1536K, 4
allfft	xfftclm	1536K, 2
	xfftclm	1536K, 1
allfft	xfftclm	1536Kp, 4
allfft	xfftclm	1536Kp, 2
	xfftclm	1536Kp, 1
allfft	xfftclm	1792K, 4
allfft	xfftclm	1792K, 2
	xfftclm	1792K, 1
allfft	xfftclm	1792K, 91
allfft	xfftclm	1792K, 0
allfft	xfftclm	1792K, 90
allfft	xfftclm	2048K, 4
allfft	xfftclm	2048K, 2
	xfftclm	2048K, 1
allfft	xfftclm	2048K, 91
allfft	xfftclm	2048K, 0
allfft	xfftclm	2048K, 90
allfft	xfftclm	2048Kp, 4
allfft	xfftclm	2048Kp, 2
	xfftclm	2048Kp, 1
allfft	xfftclm	2560K, 4
allfft	xfftclm	2560K, 2
	xfftclm	2560K, 1
allfft	xfftclm	3072K, 4
allfft	xfftclm	3072K, 2
	xfftclm	3072K, 1
	xfftclm	3072Kp, 4
allfft	xfftclm	3072Kp, 2
	xfftclm	3072Kp, 1
allfft	xfftclm	3584K, 4
allfft	xfftclm	3584K, 2
	xfftclm	3584K, 1
allfft	xfftclm	4096K, 4
allfft	xfftclm	4096K, 2
	xfftclm	4096K, 1
allfft	xfftclm	4096K, 0
allfft	xfftclm	4096Kp, 4
allfft	xfftclm	4096Kp, 2
	xfftclm	4096Kp, 1
	xfftclm	5M, 4
allfft	xfftclm	5M, 2
	xfftclm	5M, 1
allfft	xfftclm	5M, 0
	xfftclm	6M, 4
allfft	xfftclm	6M, 2
	xfftclm	6M, 1
allfft	xfftclm	6M, 0
	xfftclm	6Mp, 4
allfft	xfftclm	6Mp, 2
	xfftclm	6Mp, 1
	xfftclm	7M, 4
allfft	xfftclm	7M, 2
	xfftclm	7M, 1
allfft	xfftclm	7M, 0
	xfftclm	8M, 4
allfft	xfftclm	8M, 2
	xfftclm	8M, 1
allfft	xfftclm	8M, 0
	xfftclm	8Mp, 4
allfft	xfftclm	8Mp, 2
	xfftclm	8Mp, 1
allfft	xfftclm	10M, 2
allfft	xfftclm	10M, 1
allfft	xfftclm	10M, 0
allfft	xfftclm	12M, 2
allfft	xfftclm	12M, 1
allfft	xfftclm	12M, 0
allfft	xfftclm	12Mp, 2
allfft	xfftclm	12Mp, 1
allfft	xfftclm	12Mp, 0
allfft	xfftclm	14M, 2
allfft	xfftclm	14M, 1
allfft	xfftclm	14M, 0
allfft	xfftclm	16M, 2
allfft	xfftclm	16M, 1
allfft	xfftclm	16M, 0
allfft	xfftclm	16Mp, 2
allfft	xfftclm	16Mp, 1
allfft	xfftclm	16Mp, 0
ELSE
allfft	xfftclm	640K, 8
allfft	xfftclm	640K, 4
	xfftclm	640K, 2
allfft	xfftclm	640K, 1
allfft	xfftclm	768K, 4
	xfftclm	768K, 2
allfft	xfftclm	768K, 1
allfft	xfftclm	768Kp, 4
	xfftclm	768Kp, 2
allfft	xfftclm	896K, 4
	xfftclm	896K, 2
allfft	xfftclm	896K, 1
allfft	xfftclm	1024K, 4
allfft	xfftclm	1024K, 2
	xfftclm	1024K, 1
	xfftclm	1024Kp, 4
allfft	xfftclm	1024Kp, 2
	xfftclm	1024Kp, 1
	xfftclm	1280K, 4
allfft	xfftclm	1280K, 2
allfft	xfftclm	1280K, 1
	xfftclm	1536K, 4
allfft	xfftclm	1536K, 2
	xfftclm	1536K, 1
	xfftclm	1536Kp, 4
allfft	xfftclm	1536Kp, 2
	xfftclm	1536Kp, 1
	xfftclm	1792K, 4
allfft	xfftclm	1792K, 2
	xfftclm	1792K, 1
allfft	xfftclm	1792K, 91
allfft	xfftclm	1792K, 0
allfft	xfftclm	1792K, 90
	xfftclm	2048K, 4
	xfftclm	2048K, 2
	xfftclm	2048K, 1
allfft	xfftclm	2048K, 91
allfft	xfftclm	2048K, 0
allfft	xfftclm	2048K, 90
	xfftclm	2048Kp, 4
allfft	xfftclm	2048Kp, 2
	xfftclm	2048Kp, 1
allfft	xfftclm	2048Kp, 91
allfft	xfftclm	2048Kp, 0
allfft	xfftclm	2048Kp, 90
	xfftclm	2560K, 4
allfft	xfftclm	2560K, 2
	xfftclm	2560K, 1
	xfftclm	3072K, 4
allfft	xfftclm	3072K, 2
	xfftclm	3072K, 1
	xfftclm	3072Kp, 4
	xfftclm	3072Kp, 2
	xfftclm	3072Kp, 1
	xfftclm	3584K, 4
	xfftclm	3584K, 2
	xfftclm	3584K, 1
	xfftclm	4096K, 4
	xfftclm	4096K, 2
	xfftclm	4096K, 1
allfft	xfftclm	4096K, 0
	xfftclm	4096Kp, 4
	xfftclm	4096Kp, 2
	xfftclm	4096Kp, 1
allfft	xfftclm	4096Kp, 0
	xfftclm	5M, 4
allfft	xfftclm	5M, 2
	xfftclm	5M, 1
allfft	xfftclm	5M, 0
	xfftclm	6M, 4
	xfftclm	6M, 2
	xfftclm	6M, 1
allfft	xfftclm	6M, 0
	xfftclm	6Mp, 4
	xfftclm	6Mp, 2
	xfftclm	6Mp, 1
allfft	xfftclm	6Mp, 0
	xfftclm	7M, 4
	xfftclm	7M, 2
	xfftclm	7M, 1
allfft	xfftclm	7M, 0
	xfftclm	8M, 4
	xfftclm	8M, 2
	xfftclm	8M, 1
allfft	xfftclm	8M, 0
	xfftclm	8Mp, 4
	xfftclm	8Mp, 2
	xfftclm	8Mp, 1
allfft	xfftclm	8Mp, 0
allfft	xfftclm	10M, 2
allfft	xfftclm	10M, 1
allfft	xfftclm	10M, 0
allfft	xfftclm	12M, 2
allfft	xfftclm	12M, 1
allfft	xfftclm	12M, 0
allfft	xfftclm	12Mp, 2
allfft	xfftclm	12Mp, 1
allfft	xfftclm	12Mp, 0
allfft	xfftclm	14M, 2
allfft	xfftclm	14M, 1
allfft	xfftclm	14M, 0
allfft	xfftclm	16M, 2
allfft	xfftclm	16M, 1
allfft	xfftclm	16M, 0
allfft	xfftclm	16Mp, 2
allfft	xfftclm	16Mp, 1
allfft	xfftclm	16Mp, 0
ENDIF

_TEXT	ENDS
END

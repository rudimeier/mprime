; Copyright 2001-2007 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; FFT sizes of 40K and above are supported.
; This code does two passes, 11 levels on the second pass.
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
INCLUDE xfft3.mac
INCLUDE	xlucas.mac
INCLUDE xmult.mac
INCLUDE xpass1.mac
INCLUDE xpass1sc.mac

EXTRN	pass1_aux_entry_point_return:PROC
EXTRN	xgw_finish_fft:PROC
EXTRN	xgw_carries:PROC
EXTRN	xgw_finish_mult:PROC

EXTRNP	xpass2_11_levels
EXTRNP	xpass2_11_levels_np

_TEXT SEGMENT

;; Pass 2 does 11 FFT levels 2 sets of data (2 * 2^11 complex values =
;; 2^13 doubles = 64KB).

PREFETCHING = 1
blkdst = (8*(8192+128)+GAP2_11_4)
xpass2_levels = 11

;; All the FFT routines for each FFT length

	EXPANDING = 2
	xfft	40K
	xfft	48K
	xfft	48Kp
	xfft	56K
	xfft	64K
	xfft	64Kp
	xfft	80K
	xfft	96K
	xfft	96Kp
	xfft	112K
	xfft	128K
	xfft	128Kp
allfft	xfft	160K
allfft	xfft	192K
	xfft	192Kp
allfft	xfft	224K
allfft	xfft	256K
	xfft	256Kp
allfft	xfft	320K
allfft	xfft	384K
	xfft	384Kp
allfft	xfft	448K
allfft	xfft	512K
	xfft	512Kp

IFDEF AMD
allfft	xfftclm	320K, 4
allfft	xfftclm	384K, 4
allfft	xfftclm	384Kp, 4
allfft	xfftclm	448K, 4
allfft	xfftclm	512K, 4
allfft	xfftclm	512K, 2
allfft	xfftclm	512K, 1
allfft	xfftclm	512Kp, 4
allfft	xfftclm	640K, 8
allfft	xfftclm	640K, 4
allfft	xfftclm	640K, 2
allfft	xfftclm	640K, 1
allfft	xfftclm	768K, 4
allfft	xfftclm	768K, 2
allfft	xfftclm	768K, 1
allfft	xfftclm	768K, 0
allfft	xfftclm	768Kp, 4
allfft	xfftclm	768Kp, 2
allfft	xfftclm	896K, 4
allfft	xfftclm	896K, 2
allfft	xfftclm	896K, 1
allfft	xfftclm	896K, 0
allfft	xfftclm	1024K, 4
allfft	xfftclm	1024K, 2
allfft	xfftclm	1024K, 1
allfft	xfftclm	1024K, 0
allfft	xfftclm	1024Kp, 4
allfft	xfftclm	1024Kp, 2
allfft	xfftclm	1280K, 4
allfft	xfftclm	1280K, 2
allfft	xfftclm	1280K, 1
allfft	xfftclm	1280K, 91			; No prefetching in pass 2
allfft	xfftclm	1280K, 0
allfft	xfftclm	1280K, 90			; No prefetching in pass 2
allfft	xfftclm	1536K, 4
allfft	xfftclm	1536K, 2
allfft	xfftclm	1536K, 1
allfft	xfftclm	1536K, 91			; No prefetching in pass 2
allfft	xfftclm	1536K, 0
allfft	xfftclm	1536K, 90			; No prefetching in pass 2
allfft	xfftclm	1536Kp, 4
allfft	xfftclm	1536Kp, 2
allfft	xfftclm	1536Kp, 1
allfft	xfftclm	1792K, 4
allfft	xfftclm	1792K, 2
allfft	xfftclm	1792K, 1
allfft	xfftclm	1792K, 91			; No prefetching in pass 2
allfft	xfftclm	1792K, 0
allfft	xfftclm	1792K, 90			; No prefetching in pass 2
allfft	xfftclm	2048K, 4
allfft	xfftclm	2048K, 2
allfft	xfftclm	2048K, 1
allfft	xfftclm	2048K, 91			; No prefetching in pass 2
allfft	xfftclm	2048K, 0
allfft	xfftclm	2048K, 90			; No prefetching in pass 2
allfft	xfftclm	2048Kp, 4
allfft	xfftclm	2048Kp, 2
allfft	xfftclm	2048Kp, 1
allfft	xfftclm	2048Kp, 0
allfft	xfftclm	2560K, 2
allfft	xfftclm	2560K, 1
allfft	xfftclm	2560K, 0
allfft	xfftclm	3072K, 2
allfft	xfftclm	3072K, 1
allfft	xfftclm	3072K, 0
allfft	xfftclm	3072Kp, 2
allfft	xfftclm	3072Kp, 1
allfft	xfftclm	3072Kp, 0
allfft	xfftclm	3584K, 2
allfft	xfftclm	3584K, 1
allfft	xfftclm	3584K, 0
allfft	xfftclm	4096K, 2
allfft	xfftclm	4096K, 1
allfft	xfftclm	4096K, 0
allfft	xfftclm	4096Kp, 2
allfft	xfftclm	4096Kp, 1
allfft	xfftclm	4096Kp, 0
allfft	xfftclm	5M, 1
allfft	xfftclm	5M, 0
allfft	xfftclm	6M, 1
allfft	xfftclm	6M, 0
allfft	xfftclm	6Mp, 1
allfft	xfftclm	6Mp, 0
allfft	xfftclm	7M, 1
allfft	xfftclm	7M, 0
allfft	xfftclm	8M, 1
allfft	xfftclm	8M, 0
allfft	xfftclm	8Mp, 1
allfft	xfftclm	8Mp, 0
ELSE
allfft	xfftclm	320K, 4
allfft	xfftclm	384K, 4
allfft	xfftclm	384Kp, 4
allfft	xfftclm	448K, 4
allfft	xfftclm	512K, 4
allfft	xfftclm	512K, 2
allfft	xfftclm	512K, 1
allfft	xfftclm	512Kp, 4
allfft	xfftclm	640K, 8
	xfftclm	640K, 4
allfft	xfftclm	640K, 2
	xfftclm	640K, 1
allfft	xfftclm	640K, 91			; No prefetching in pass 2
	xfftclm	768K, 4
	xfftclm	768K, 2
allfft	xfftclm	768K, 1
allfft	xfftclm	768K, 91			; No prefetching in pass 2
allfft	xfftclm	768K, 0
allfft	xfftclm	768K, 90			; No prefetching in pass 2
	xfftclm	768Kp, 4
	xfftclm	768Kp, 2
allfft	xfftclm	768Kp, 1
allfft	xfftclm	768Kp, 91			; No prefetching in pass 2
	xfftclm	768Kp, 0
allfft	xfftclm	768Kp, 90			; No prefetching in pass 2
	xfftclm	896K, 4
	xfftclm	896K, 2
allfft	xfftclm	896K, 1
allfft	xfftclm	896K, 91			; No prefetching in pass 2
	xfftclm	896K, 0
allfft	xfftclm	896K, 90			; No prefetching in pass 2
	xfftclm	1024K, 4
	xfftclm	1024K, 2
allfft	xfftclm	1024K, 1
allfft	xfftclm	1024K, 91			; No prefetching in pass 2
	xfftclm	1024K, 0
allfft	xfftclm	1024K, 90			; No prefetching in pass 2
	xfftclm	1024Kp, 4
	xfftclm	1024Kp, 2
allfft	xfftclm	1024Kp, 1
allfft	xfftclm	1024Kp, 91			; No prefetching in pass 2
	xfftclm	1024Kp, 0
allfft	xfftclm	1024Kp, 90			; No prefetching in pass 2
allfft	xfftclm	1280K, 4
	xfftclm	1280K, 2
allfft	xfftclm	1280K, 1
allfft	xfftclm	1280K, 91			; No prefetching in pass 2
	xfftclm	1280K, 0
allfft	xfftclm	1280K, 90			; No prefetching in pass 2
allfft	xfftclm	1536K, 4
allfft	xfftclm	1536K, 2
	xfftclm	1536K, 1
allfft	xfftclm	1536K, 91			; No prefetching in pass 2
	xfftclm	1536K, 0
allfft	xfftclm	1536K, 90			; No prefetching in pass 2
allfft	xfftclm	1536Kp, 4
allfft	xfftclm	1536Kp, 2
	xfftclm	1536Kp, 1
allfft	xfftclm	1536Kp, 91			; No prefetching in pass 2
	xfftclm	1536Kp, 0
allfft	xfftclm	1536Kp, 90			; No prefetching in pass 2
allfft	xfftclm	1792K, 4
allfft	xfftclm	1792K, 2
	xfftclm	1792K, 1
allfft	xfftclm	1792K, 91			; No prefetching in pass 2
	xfftclm	1792K, 0
allfft	xfftclm	1792K, 90			; No prefetching in pass 2
allfft	xfftclm	2048K, 4
allfft	xfftclm	2048K, 2
allfft	xfftclm	2048K, 1
allfft	xfftclm	2048K, 91			; No prefetching in pass 2
	xfftclm	2048K, 0
allfft	xfftclm	2048K, 90			; No prefetching in pass 2
allfft	xfftclm	2048Kp, 4
allfft	xfftclm	2048Kp, 2
allfft	xfftclm	2048Kp, 1
allfft	xfftclm	2048Kp, 91			; No prefetching in pass 2
	xfftclm	2048Kp, 0
allfft	xfftclm	2048Kp, 90			; No prefetching in pass 2
allfft	xfftclm	2560K, 2
allfft	xfftclm	2560K, 1
allfft	xfftclm	2560K, 0
allfft	xfftclm	3072K, 2
allfft	xfftclm	3072K, 1
allfft	xfftclm	3072K, 0
allfft	xfftclm	3072Kp, 2
allfft	xfftclm	3072Kp, 1
allfft	xfftclm	3072Kp, 0
allfft	xfftclm	3584K, 2
allfft	xfftclm	3584K, 1
allfft	xfftclm	3584K, 0
allfft	xfftclm	4096K, 2
allfft	xfftclm	4096K, 1
allfft	xfftclm	4096K, 0
allfft	xfftclm	4096Kp, 2
allfft	xfftclm	4096Kp, 1
allfft	xfftclm	4096Kp, 0
allfft	xfftclm	5M, 1
allfft	xfftclm	5M, 0
allfft	xfftclm	6M, 1
allfft	xfftclm	6M, 0
allfft	xfftclm	6Mp, 1
allfft	xfftclm	6Mp, 0
allfft	xfftclm	7M, 1
allfft	xfftclm	7M, 0
allfft	xfftclm	8M, 1
allfft	xfftclm	8M, 0
allfft	xfftclm	8Mp, 1
allfft	xfftclm	8Mp, 0
ENDIF

_TEXT	ENDS
END

; Copyright 2001-2007 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This routine implements a discrete weighted transform to quickly multiply
; two numbers.
;
; This code uses Pentium 4's SSE2 instructions for very fast FFTs.
; FFT sizes up to 8K are supported in a single pass.
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
INCLUDE xfft1.mac
INCLUDE	xlucas.mac
INCLUDE xmult.mac

;; All the FFT routines for each FFT length

_TEXT SEGMENT

	EXPANDING = 1
	xfft	32
	xfft	48
	xfft	64
	xfft	80
	xfft	96
	xfft	112
	xfft	128
	xfft	160
	xfft	192
	xfft	224
	xfft	256
	xfft	320
	xfft	384
	xfft	448
	xfft	512
	xfft	640
	xfft	768
	xfft	896
	xfft	1024
	xfft	1280
	xfft	1536
	xfft	1792
	xfft	2048
	xfft	2560
	xfft	3072
	xfft	3584
	xfft	4096
	xfft	5120
	xfft	6144
	xfft	7168
	xfft	8192

	xfft	32p
	xfft	48p
	xfft	64p
	xfft	96p
	xfft	128p
	xfft	192p
	xfft	256p
	xfft	384p
	xfft	512p
	xfft	768p
	xfft	1024p
	xfft	1536p
	xfft	2048p
	xfft	3072p
	xfft	4096p
	xfft	6144p
	xfft	8192p

_TEXT	ENDS
END

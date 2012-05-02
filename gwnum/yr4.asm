; Copyright 2011 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Assemble the small AVX FFTs.
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

yfft_type TEXTEQU <r4>

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE yarch.mac
INCLUDE ybasics.mac
INCLUDE ymult.mac
INCLUDE yonepass.mac
INCLUDE yr4.mac

_TEXT SEGMENT

;; Implement the small one pass FFTs, for now only assemble one-pass FFTs for the BLEND architecture

IF @INSTR(,%yarch,<BLEND>) NE 0

PREFETCHING = 0

	yonepass 32, 0
	yonepass 64, 0
	yonepass 96, 0
;	yonepass 112, 0
	yonepass 128, 0
	yonepass 160, 0
	yonepass 192, 0
;	yonepass 224, 0
	yonepass 256, 0
	yonepass 320, 0
	yonepass 384, 0
;	yonepass 448, 0
	yonepass 512, 0
	yonepass 640, 0
	yonepass 768, 0
;	yonepass 896, 0
	yonepass 1K, 0
	yonepass 1280, 0
	yonepass 1536, 0
;	yonepass 1792, 0
	yonepass 2K, 0
	yonepass 2560, 0
	yonepass 3K, 0
;	yonepass 3584, 0
	yonepass 4K, 0
	yonepass 5K, 0
	yonepass 6K, 0
;	yonepass 7168, 0
	yonepass 8K, 0
	yonepass 10K, 0
	yonepass 12K, 0
	yonepass 16K, 0
	yonepass 20K, 0
	yonepass 18K, 0
	yonepass 24K, 0
	yonepass 32K, 0

	yonepass 32, 1
	yonepass 64, 1
	yonepass 96, 1
	yonepass 128, 1
	yonepass 160, 1
	yonepass 192, 1
	yonepass 256, 1
	yonepass 320, 1
	yonepass 384, 1
	yonepass 512, 1
	yonepass 640, 1
	yonepass 768, 1
	yonepass 1K, 1
	yonepass 1280, 1
	yonepass 1536, 1
	yonepass 2K, 1
	yonepass 2560, 1
	yonepass 3K, 1
	yonepass 4K, 1
	yonepass 5K, 1
	yonepass 6K, 1
	yonepass 8K, 1
	yonepass 10K, 1
	yonepass 12K, 1
	yonepass 16K, 1
	yonepass 20K, 1
	yonepass 18K, 1
	yonepass 24K, 1
	yonepass 32K, 1
ENDIF

_TEXT	ENDS
END

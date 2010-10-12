; Copyright 1995-2010 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; Time low-level operations for optimizing macros and comparing CPUs
;

	TITLE   setup

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE	unravel.mac
INCLUDE extrn.mac
INCLUDE memory.mac
INCLUDE pfa.mac
INCLUDE lucas.mac
INCLUDE xarch.mac
INCLUDE xbasics.mac
INCLUDE xmult.mac
INCLUDE hg.mac
INCLUDE r4.mac


IFDEF X86_64
X87_CASES	EQU	0
ELSE
X87_CASES	EQU	13
ENDIF
SSE2_CASES	EQU	216

loopent	MACRO	y,z		; Create a entry in the loop entry table
	DP	&y&z
	ENDM
looptab	MACRO	y, cnt		; Create the loop entry table
	x = 0
	REPT	cnt
	loopent	y, %x
	x = x + 1
	ENDM
	ENDM

;; Macros from p4notes.doc

read1	MACRO	mem, c			;; Bytes to read
	LOCAL loop1, loop2
	cnt = mem/64		        ;; Read 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:	mov     eax, [rsi]		;; Read one cache line
	mov     eax, [rsi+4]		;; 4 bytes at a time
	mov     eax, [rsi+8]
	mov     eax, [rsi+12]
	mov     eax, [rsi+16]
	mov     eax, [rsi+20]
	mov     eax, [rsi+24]
	mov     eax, [rsi+28]
	mov     eax, [rsi+32]
	mov     eax, [rsi+36]
	mov     eax, [rsi+40]
	mov     eax, [rsi+44]
	mov     eax, [rsi+48]
	mov     eax, [rsi+52]
	mov     eax, [rsi+56]
	mov     eax, [rsi+60]
	lea     rsi, [rsi+64]			; Next cache line
	sub     ecx, 1
	jnz     loop1
	lea     rsi, [rsi-mem]		; Restore esi
	dec	edx
	jnz	loop2
	ENDM

read2	MACRO	mem, c			;; Bytes to read
	LOCAL loop1, loop2
        cnt = mem/64                    ; 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  xmm1, [rsi]			; Read one cache line
        movapd  xmm1, [rsi+16]		; 16 bytes at a time
        movapd  xmm1, [rsi+32]
        movapd  xmm1, [rsi+48]
        lea     rsi, [rsi+64]			; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

write1	MACRO	mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/64                    ; 64 bytes per iteration
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  [rsi], xmm1			; Write one cache line
        movapd  [rsi+16], xmm1		; 16 bytes at a time
        movapd  [rsi+32], xmm1
        movapd  [rsi+48], xmm1
        lea     rsi, [rsi+64]			; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

write2	MACRO	mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/(4*128)               ; 128 bytes per iteration
        dist = 64
	mov	edx, c
	sub     ebx, ebx
loop2:	mov     ecx, cnt
loop1:  movapd  [rsi+0*dist], xmm1      ; Write 8 cache lines
        movapd  [rsi+1*dist], xmm1
        movapd  [rsi+2*dist], xmm1
        movapd  [rsi+3*dist], xmm1
        movapd  [rsi+4*dist], xmm1
        movapd  [rsi+5*dist], xmm1
        movapd  [rsi+6*dist], xmm1
        movapd  [rsi+7*dist], xmm1
        lea     rsi, [rsi+16]           ; Same cache lines
        add     bl, 256/4			; 4 inner loop iterations
        jnc     loop1
        lea     rsi, [rsi-4*16+8*dist]  ; Next set of 8 cache lines
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

x4cl_empty MACRO srcreg,srcinc,d1,d2,screg,scoff
	xload	xmm0, [srcreg+0]
	xload	xmm1, [srcreg+32]
	xload	xmm2, [srcreg+d1+0]
	xload	xmm3, [srcreg+d1+32]
	xload	xmm4, [srcreg+d2+0]
	xload	xmm5, [srcreg+d2+32]
	xload	xmm6, [srcreg+d2+d1+0]
	xload	xmm7, [srcreg+d2+d1+32]
	xstore	[srcreg+0], xmm0
	xstore	[srcreg+32], xmm1
	xload	xmm0, [srcreg+0+16]
	xload	xmm1, [srcreg+32+16]
	xstore	[srcreg+16], xmm2
	xstore	[srcreg+48], xmm3
	xload	xmm2, [srcreg+d1+0+16]
	xload	xmm3, [srcreg+d1+32+16]
	xstore	[srcreg+d1], xmm4
	xstore	[srcreg+d1+16], xmm5
	xstore	[srcreg+d1+32], xmm6
	xstore	[srcreg+d1+48], xmm7
	xload	xmm4, [srcreg+d2+0+16]
	xload	xmm5, [srcreg+d2+32+16]
	xload	xmm6, [srcreg+d2+d1+0+16]
	xload	xmm7, [srcreg+d2+d1+32+16]
	xstore	[srcreg+d2], xmm0
	xstore	[srcreg+d2+16], xmm1
	xstore	[srcreg+d2+32], xmm2
	xstore	[srcreg+d2+48], xmm3
	xstore	[srcreg+d2+d1], xmm4
	xstore	[srcreg+d2+d1+16], xmm5
	xstore	[srcreg+d2+d1+32], xmm6
	xstore	[srcreg+d2+d1+48], xmm7
	bump	srcreg, srcinc
	ENDM

g4cl_empty MACRO srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg,scoff
	xload	xmm0, [srcreg+0]
	xload	xmm1, [srcreg+32]
	xload	xmm2, [srcreg+d1+0]
	xload	xmm3, [srcreg+d1+32]
	xload	xmm4, [srcreg+d2+0]
	xload	xmm5, [srcreg+d2+32]
	xload	xmm6, [srcreg+d2+d1+0]
	xload	xmm7, [srcreg+d2+d1+32]
	xstore	[dstreg+0], xmm0
	xstore	[dstreg+32], xmm1
	xload	xmm0, [srcreg+0+16]
	xload	xmm1, [srcreg+32+16]
	xstore	[dstreg+16], xmm2
	xstore	[dstreg+48], xmm3
	xload	xmm2, [srcreg+d1+0+16]
	xload	xmm3, [srcreg+d1+32+16]
	xstore	[dstreg+e1], xmm4
	xstore	[dstreg+e1+16], xmm5
	xstore	[dstreg+e1+32], xmm6
	xstore	[dstreg+e1+48], xmm7
	xload	xmm4, [srcreg+d2+0+16]
	xload	xmm5, [srcreg+d2+32+16]
	xload	xmm6, [srcreg+d2+d1+0+16]
	xload	xmm7, [srcreg+d2+d1+32+16]
	bump	srcreg, srcinc
	xstore	[dstreg+e2], xmm0
	xstore	[dstreg+e2+16], xmm1
	xstore	[dstreg+e2+32], xmm2
	xstore	[dstreg+e2+48], xmm3
	xstore	[dstreg+e2+e1], xmm4
	xstore	[dstreg+e2+e1+16], xmm5
	xstore	[dstreg+e2+e1+32], xmm6
	xstore	[dstreg+e2+e1+48], xmm7
	bump	dstreg, dstinc
	ENDM

g4cl_empty_nt MACRO srcreg,srcinc,d1,d2,dstreg,dstinc,e1,e2,screg,scoff
	xload	xmm0, [srcreg+0]
	xload	xmm1, [srcreg+32]
	xload	xmm2, [srcreg+d1+0]
	xload	xmm3, [srcreg+d1+32]
	xload	xmm4, [srcreg+d2+0]
	xload	xmm5, [srcreg+d2+32]
	xload	xmm6, [srcreg+d2+d1+0]
	xload	xmm7, [srcreg+d2+d1+32]
	movntpd	[dstreg+0], xmm0
	movntpd	[dstreg+32], xmm1
	xload	xmm0, [srcreg+0+16]
	xload	xmm1, [srcreg+32+16]
	movntpd	[dstreg+16], xmm2
	movntpd	[dstreg+48], xmm3
	xload	xmm2, [srcreg+d1+0+16]
	xload	xmm3, [srcreg+d1+32+16]
	movntpd	[dstreg+e1], xmm4
	movntpd	[dstreg+e1+16], xmm5
	movntpd	[dstreg+e1+32], xmm6
	movntpd	[dstreg+e1+48], xmm7
	xload	xmm4, [srcreg+d2+0+16]
	xload	xmm5, [srcreg+d2+32+16]
	xload	xmm6, [srcreg+d2+d1+0+16]
	xload	xmm7, [srcreg+d2+d1+32+16]
	bump	srcreg, srcinc
	movntpd	[dstreg+e2], xmm0
	movntpd	[dstreg+e2+16], xmm1
	movntpd	[dstreg+e2+32], xmm2
	movntpd	[dstreg+e2+48], xmm3
	movntpd	[dstreg+e2+e1], xmm4
	movntpd	[dstreg+e2+e1+16], xmm5
	movntpd	[dstreg+e2+e1+32], xmm6
	movntpd	[dstreg+e2+e1+48], xmm7
	bump	dstreg, dstinc
	ENDM


readwrite1 MACRO mem, c			;; Bytes to write
	LOCAL loop1, loop2
        cnt = mem/64
	mov	edx, c
loop2:	mov     ecx, cnt
loop1:  movapd  xmm0, [rsi]             ; Read one cache line
        movapd  xmm1, [rsi+16]
        movapd  xmm2, [rsi+32]
        movapd  xmm3, [rsi+48]
        subpd   xmm0, xmm0              ; Operate on the data
        pxor    xmm1, xmm1
        subpd   xmm2, xmm2
        pxor    xmm3, xmm3
        movapd  [rsi], xmm0             ; Write the cache line
        movapd  [rsi+16], xmm1
        movapd  [rsi+32], xmm2
        movapd  [rsi+48], xmm3
        lea     rsi, [rsi+64]           ; Next cache line
        sub     ecx, 1
        jnz     loop1
        lea     rsi, [rsi-mem]
	dec	edx
	jnz	loop2
	ENDM

;; Time one of the basic FFT building blocks

x87mac	MACRO	memused, memarea, ops:vararg
	LOCAL	ss0a, ss0b
	inner_iters = memarea / memused
	outer_iters = 10000 / inner_iters
	mov	eax, outer_iters
	mov	SRCARG, rdi		;; Save work buf addr
ss0a:	mov	rdi, SRCARG		;; Reload work buf addr
	lea	rsi, [rdi+4096]
	mov	ecx, inner_iters
ss0b:	disp	&ops
	lea	rsi, [rsi+memused]	;; Next source pointer
	lea	rdi, [rdi+SCD]		;; Next sine/cosine pointer
	dec	ecx
	jnz	ss0b
	dec	eax
	jnz	ss0a
	ENDM

sse2mac MACRO	lab, memused, memarea, ops:vararg
	LOCAL	ss0a, ss0b
	inner_iters = memarea / memused
	outer_iters = 10000 / inner_iters
	odd_iters = 10000 - inner_iters * outer_iters
	IF odd_iters EQ 0
	odd_iters = inner_iters
	ELSE
	outer_iters = outer_iters + 1
	ENDIF
lab:	mov	rbx, 0			;; Offset for some sse2 macros (s.b. non-zero for with_mult macros)
	mov	rbp, 524288+256		;; Offset for mulf sse2 macros
	mov	eax, outer_iters
	mov	ecx, odd_iters
	mov	SRCARG, rdi		;; Save work buf addr
ss0a:	mov	rdi, SRCARG		;; Reload work buf addr (sincos data)
	lea	rsi, [rdi+4096]		;; Source & dest ptr
	lea	rdx, [rsi+524288+256]	;; Destination for "g" macros
align 16
ss0b:	&ops
IF memused NE 192
;	lea	rdi, [rdi+2*XMM_SCD]	;; Next sine/cosine pointer
ELSE
;	lea	rdi, [rdi+XMM_SCD1]	;; Next sine/cosine pointer
ENDIF
	dec	ecx
	jnz	ss0b
	mov	ecx, inner_iters
	dec	eax
	jnz	ss0a
	jmp	exit
	ENDM
sse2macbx MACRO	lab, memused, memarea, ops:vararg
	LOCAL	ss0a, ss0b
	inner_iters = memarea / memused
	outer_iters = 10000 / inner_iters
	odd_iters = 10000 - inner_iters * outer_iters
	IF odd_iters EQ 0
	odd_iters = inner_iters
	ELSE
	outer_iters = outer_iters + 1
	ENDIF
lab:	mov	rbx, 0;;262144+128		;; Offset for some sse2 macros (s.b. non-zero for with_mult macros)
	mov	rbp, 524288+256		;; Offset for mulf sse2 macros
	mov	eax, outer_iters
	mov	ecx, odd_iters
	mov	SRCARG, rdi		;; Save work buf addr
ss0a:	mov	rdi, SRCARG		;; Reload work buf addr (sincos data)
	lea	rsi, [rdi+4096]		;; Source & dest ptr
	lea	rdx, [rsi+524288+256]	;; Destination for "g" macros
align 16
ss0b:	&ops
IF memused NE 192
;	lea	rdi, [rdi+2*XMM_SCD]	;; Next sine/cosine pointer
ELSE
;	lea	rdi, [rdi+XMM_SCD1]	;; Next sine/cosine pointer
ENDIF
	dec	ecx
	jnz	ss0b
	mov	ecx, inner_iters
	dec	eax
	jnz	ss0a
	jmp	exit
	ENDM

_TEXT	SEGMENT

x87table: looptab case, X87_CASES
sse2table: looptab sscase, SSE2_CASES

; gwtimeit (asm_data)
;	Time a mini benchmark
; Windows 32-bit (_gwtimeit)
; Linux 32-bit (gwtimeit)
;	Parameter asm_data = [esp+4]
; Windows 64-bit (gwtimeit)
;	Parameter asm_data = rcx
; Linux 64-bit (gwtimeit)
;	Parameter asm_data = rdi

PROCFL	gwtimeit
	ad_prolog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

	mov	rdi, PPTR [AD_BASE]	; Load work buf address
	mov	rsi, rdi

	mov	edx, DPTR [AD_BASE+8]	; Load n (which test to run)

	sub	ebx, ebx		; Clear registers
	sub	ecx, ecx
	sub	ebp, ebp

	mov	eax, X87_CASES
	cmp	edx, -1			; -1 = get num x87 cases
	je	exit
	mov	eax, SSE2_CASES
	cmp	edx, -2			; -2 = get num sse2 cases
	je	exit

	cmp	edx, 1000		; Tests above 1000 are SSE2 code
	jl	short x87

	subpd	xmm0, xmm0		; Clear XMM registers
	subpd	xmm1, xmm1
	subpd	xmm2, xmm2
	subpd	xmm3, xmm3
	subpd	xmm4, xmm4
	subpd	xmm5, xmm5
	subpd	xmm6, xmm6
	subpd	xmm7, xmm7
IFDEF X86_64
	subpd	xmm8, xmm8
	subpd	xmm9, xmm9
	subpd	xmm10, xmm10
	subpd	xmm11, xmm11
	subpd	xmm12, xmm12
	subpd	xmm13, xmm13
	subpd	xmm14, xmm14
	subpd	xmm15, xmm15
ENDIF

; Jump to desired test case

	sub	rdx, 1000
	mov	rax, OFFSET sse2table
	mov	rax, [rax+rdx*SZPTR]; Get address of test to execute
	jmp	rax
x87:	mov	rax, OFFSET x87table
	mov	rax, [rax+rdx*SZPTR]; Get address of test to execute
	jmp	rax

; Time the "do-nothing" case.

IFNDEF X86_64
case0:
	jmp	exit

; Time the loop of "do-nothing" case.  1000 iterations.

case1:
	mov	ecx, 1000
c1a:	dec	ecx
	jnz	c1a
	jmp	exit

; This code reads a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will read from the L1 cache
; only, 96KB will read from the L2 cache only, and 2MB will test reading
; from main memory.

case2:	read1	4096, 1000	; Read 4KB
	jmp	exit

case3:	read1	96*1024, 100	; Read 96KB
	jmp	exit

case4:	read1	2048*1024, 2	; Read 2MB
	jmp	exit

case5:	x87mac	64, 4096, eight_reals_fft, 8, 16, 32
	jmp	exit

case6:	x87mac	64, 100000, eight_reals_fft, 8, 16, 32
	jmp	exit

case7:	x87mac	64, 4096, eight_reals_unfft, 8, 16, 32
	jmp	exit

case8:	x87mac	64, 100000, eight_reals_unfft, 8, 16, 32
	jmp	exit

case9:	x87mac	64, 4096, four_complex_fft, 8, 16, 32
	jmp	exit

case10:	x87mac	64, 100000, four_complex_fft, 8, 16, 32
	jmp	exit

case11:	x87mac	64, 4096, four_complex_unfft, 8, 16, 32
	jmp	exit

case12:	x87mac	64, 100000, four_complex_unfft, 8, 16, 32
	jmp	exit
ENDIF

; This code reads a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will read from the L1 cache
; only, 96KB will read from the L2 cache only, and 32MB will test reading
; from main memory.

sscase0:
	read2	4096, 1000	; Read 4KB
	jmp	exit

sscase1:
	read2	96*1024, 100	; Read 96KB
	jmp	exit

sscase2:
	read2	32768*1024, 2	; Read 32MB
	jmp	exit

; This code writes a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 32MB will test writing
; to main memory.

sscase3:
	write1	4096, 1000	; Write 4KB
	jmp	exit

sscase4:
	write1	96*1024, 100	; Write 96KB
	jmp	exit

sscase5:
	write1	32768*1024, 2	; Write 32MB
	jmp	exit

; This code writes a block of memory non-contiguously.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 32MB will test writing
; to main memory.

sscase6:
	write2	4096, 1000	; Read 4KB
	jmp	exit

sscase7:
	write2	96*1024, 100	; Read 96KB
	jmp	exit

sscase8:
	write2	32768*1024, 2	; Read 32MB
	jmp	exit

; This code reads & writes a block of memory.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 32MB will test writing
; to main memory.

sscase9:
	readwrite1	4096, 1000	; Read 4KB
	jmp	exit

sscase10:
	readwrite1	96*1024, 100	; Read 96KB
	jmp	exit

sscase11:
	readwrite1	32768*1024, 2	; Read 32MB
	jmp	exit

; Time ~10000 iterations of the SSE2 macros in L1 and L2 caches

	sse2mac sscase12, 128, 4096, x2cl_eight_reals_fft rsi, 2*64, 64
	sse2mac sscase13, 128, 100000, x2cl_eight_reals_fft rsi, 2*64, 64
	sse2mac sscase14, 128, 4096, x2cl_eight_reals_first_fft rsi, 2*64, 64
	sse2mac sscase15, 128, 100000, x2cl_eight_reals_first_fft rsi, 2*64, 64
	sse2mac sscase16, 128, 4096, x2cl_eight_reals_fft_2 rsi, 2*64, 64
	sse2mac sscase17, 128, 100000, x2cl_eight_reals_fft_2 rsi, 2*64, 64
	sse2mac sscase18, 128, 4096, x2cl_eight_reals_fft_1 rsi, 2*64, 64
	sse2mac sscase19, 128, 100000, x2cl_eight_reals_fft_1 rsi, 2*64, 64
	sse2mac sscase20, 128, 4096, s2cl_eight_reals_first_fft rsi, 2*64, 64
	sse2mac sscase21, 128, 100000, s2cl_eight_reals_first_fft rsi, 2*64, 64
	sse2mac sscase22, 128, 4096, s2cl_eight_reals_fft_1 rsi, 2*64, 64
	sse2mac sscase23, 128, 100000, s2cl_eight_reals_fft_1 rsi, 2*64, 64
	sse2mac sscase24, 128, 4096, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	sse2mac sscase25, 128, 100000, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	sse2mac sscase26, 128, 4096, s2cl_eight_reals_fft_2_final rsi, 2*64, 64
	sse2mac sscase27, 128, 100000, s2cl_eight_reals_fft_2_final rsi, 2*64, 64
	sse2mac sscase28, 128, 4096, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	sse2mac sscase29, 128, 100000, s2cl_eight_reals_with_square_2 rsi, 2*64, 64
	sse2macbx sscase30, 128, 4096, s2cl_eight_reals_with_mult_2 rsi, 2*64, 64
	sse2macbx sscase31, 128, 100000, s2cl_eight_reals_with_mult_2 rsi, 2*64, 64
	sse2macbx sscase32, 128, 4096, s2cl_eight_reals_with_mulf_2 rsi, 2*64, 64
	sse2macbx sscase33, 128, 100000, s2cl_eight_reals_with_mulf_2 rsi, 2*64, 64
	sse2mac sscase34, 256, 4096, x4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase35, 256, 100000, x4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase36, 256, 4096, g4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64
	sse2mac sscase37, 256, 100000, g4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64

	sse2mac sscase38, 128, 4096, x2cl_eight_reals_unfft_2 rsi, 2*64, 64
	sse2mac sscase39, 128, 100000, x2cl_eight_reals_unfft_2 rsi, 2*64, 64
	sse2mac sscase40, 128, 4096, x2cl_half_eight_reals_unfft_2 rsi, 2*64, 64
	sse2mac sscase41, 128, 100000, x2cl_half_eight_reals_unfft_2 rsi, 2*64, 64
	sse2mac sscase42, 128, 4096, g2cl_eight_reals_unfft_2 rsi, 2*64, 64, rdx, 2*64, 64
	sse2mac sscase43, 128, 100000, g2cl_eight_reals_unfft_2 rsi, 2*64, 64, rdx, 2*64, 64
	sse2mac sscase44, 256, 4096, x4cl_eight_reals_last_unfft rsi, 4*64, 64, 2*64
	sse2mac sscase45, 256, 100000, x4cl_eight_reals_last_unfft rsi, 4*64, 64, 2*64
	sse2mac sscase46, 256, 4096, x4cl_eight_reals_unfft_2 rsi, 4*64, 64, 2*64
	sse2mac sscase47, 256, 100000, x4cl_eight_reals_unfft_2 rsi, 4*64, 64, 2*64
	sse2mac sscase48, 256, 4096, s4cl_eight_reals_unfft_1 rsi, 4*64, 64, 2*64
	sse2mac sscase49, 256, 100000, s4cl_eight_reals_unfft_1 rsi, 4*64, 64, 2*64

	sse2mac sscase50, 128, 4096, x2cl_two_complex_fft rsi, 2*64, 64, rdi
	sse2mac sscase51, 128, 100000, x2cl_two_complex_fft rsi, 2*64, 64, rdi
	sse2mac sscase52, 128, 4096, x2cl_two_complex_fft_in_place rsi, 2*64, 64, rdi
	sse2mac sscase53, 128, 100000, x2cl_two_complex_fft_in_place rsi, 2*64, 64, rdi

	sse2mac sscase54, 128, 4096, x2cl_two_complex_unfft rsi, 2*64, 64
	sse2mac sscase55, 128, 100000, x2cl_two_complex_unfft rsi, 2*64, 64

	sse2mac sscase56, 128, 4096, x2cl_four_complex_fft rsi, 2*64, 64
	sse2mac sscase57, 128, 100000, x2cl_four_complex_fft rsi, 2*64, 64
	sse2mac sscase58, 128, 4096, x2cl_four_complex_first_fft rsi, 2*64, 64
	sse2mac sscase59, 128, 100000, x2cl_four_complex_first_fft rsi, 2*64, 64
	sse2mac sscase60, 128, 4096, s2cl_four_complex_gpm_fft rsi, 2*64, 64
	sse2mac sscase61, 128, 100000, s2cl_four_complex_gpm_fft rsi, 2*64, 64
	sse2mac sscase62, 128, 4096, s2cl_four_complex_first_fft rsi, 2*64, 64
	sse2mac sscase63, 128, 100000, s2cl_four_complex_first_fft rsi, 2*64, 64
	sse2mac sscase64, 128, 4096, s2cl_four_complex_fft_final rsi, 2*64, 64
	sse2mac sscase65, 128, 100000, s2cl_four_complex_fft_final rsi, 2*64, 64
	sse2mac sscase66, 128, 4096, s2cl_four_complex_with_square rsi, 2*64, 64
	sse2mac sscase67, 128, 100000, s2cl_four_complex_with_square rsi, 2*64, 64
	sse2macbx sscase68, 128, 4096, s2cl_four_complex_with_mult rsi, 2*64, 64
	sse2macbx sscase69, 128, 100000, s2cl_four_complex_with_mult rsi, 2*64, 64
	sse2macbx sscase70, 128, 4096, s2cl_four_complex_with_mulf rsi, 2*64, 64
	sse2macbx sscase71, 128, 100000, s2cl_four_complex_with_mulf rsi, 2*64, 64
	sse2mac sscase72, 256, 4096, x4cl_four_complex_fft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase73, 256, 100000, x4cl_four_complex_fft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase74, 256, 4096, x4cl_four_complex_cpm_fft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase75, 256, 100000, x4cl_four_complex_cpm_fft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase76, 256, 4096, g4cl_four_complex_fft rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64
	sse2mac sscase77, 256, 100000, g4cl_four_complex_fft rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64
	sse2mac sscase78, 256, 4096, x4cl_four_complex_with_square rsi, 4*64, 64, 2*64
	sse2mac sscase79, 256, 100000, x4cl_four_complex_with_square rsi, 4*64, 64, 2*64
	sse2macbx sscase80, 256, 4096, x4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	sse2macbx sscase81, 256, 100000, x4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	sse2macbx sscase82, 256, 4096, x4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64
	sse2macbx sscase83, 256, 100000, x4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64

	sse2mac sscase84, 128, 4096, x2cl_four_complex_unfft rsi, 2*64, 64
	sse2mac sscase85, 128, 100000, x2cl_four_complex_unfft rsi, 2*64, 64
	sse2mac sscase86, 128, 4096, g2cl_four_complex_unfft rsi, 2*64, 64, rdx, 2*64, 64
	sse2mac sscase87, 128, 100000, g2cl_four_complex_unfft rsi, 2*64, 64, rdx, 2*64, 64
	sse2mac sscase88, 256, 4096, x4cl_four_complex_unfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase89, 256, 100000, x4cl_four_complex_unfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase90, 256, 4096, x4cl_four_complex_last_unfft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase91, 256, 100000, x4cl_four_complex_last_unfft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase92, 256, 4096, s4cl_four_complex_gpm_unfft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase93, 256, 100000, s4cl_four_complex_gpm_unfft rsi, 4*64, 64, 2*64, 0
	sse2mac sscase94, 256, 4096, x4cl_four_complex_cpm_unfft rsi, 4*64, 64, 2*64
	sse2mac sscase95, 256, 100000, x4cl_four_complex_cpm_unfft rsi, 4*64, 64, 2*64

	sse2mac sscase96, 192, 4096, x3cl_six_reals_first_fft rsi, 3*64, 64
	sse2mac sscase97, 192, 100000, x3cl_six_reals_first_fft rsi, 3*64, 64
	sse2mac sscase98, 192, 4096, g3cl_six_reals_first_fft rsi, 3*64, 64, rdx, 3*64, 64
	sse2mac sscase99, 192, 100000, g3cl_six_reals_first_fft rsi, 3*64, 64, rdx, 3*64, 64
	sse2mac sscase100, 192, 4096, s3cl_six_reals_first_fft rsi, 3*64, 64
	sse2mac sscase101, 192, 100000, s3cl_six_reals_first_fft rsi, 3*64, 64
	sse2mac sscase102, 192, 4096, x3cl_six_reals_last_unfft rsi, 3*64, 64
	sse2mac sscase103, 192, 100000, x3cl_six_reals_last_unfft rsi, 3*64, 64

	sse2mac sscase104, 192, 4096, x3cl_three_complex_first_fft rsi, 3*64, 64
	sse2mac sscase105, 192, 100000, x3cl_three_complex_first_fft rsi, 3*64, 64
	sse2mac sscase106, 192, 4096, s3cl_three_complex_first_fft rsi, 3*64, 64
	sse2mac sscase107, 192, 100000, s3cl_three_complex_first_fft rsi, 3*64, 64
	sse2mac sscase108, 192, 4096, x3cl_three_complex_last_unfft rsi, 3*64, 64
	sse2mac sscase109, 192, 100000, x3cl_three_complex_last_unfft rsi, 3*64, 64

	sse2mac sscase110, 320, 4096, x5cl_five_reals_first_fft rsi, 5*64, 64
	sse2mac sscase111, 320, 100000, x5cl_five_reals_first_fft rsi, 5*64, 64
	sse2mac sscase112, 320, 4096, g5cl_five_reals_first_fft rsi, 5*64, 64, rdx, 5*64, 64
	sse2mac sscase113, 320, 100000, g5cl_five_reals_first_fft rsi, 5*64, 64, rdx, 5*64, 64
	sse2mac sscase114, 320, 4096, s5cl_five_reals_first_fft rsi, 5*64, 64
	sse2mac sscase115, 320, 100000, s5cl_five_reals_first_fft rsi, 5*64, 64
	sse2mac sscase116, 320, 4096, x5cl_five_reals_last_unfft rsi, 5*64, 64
	sse2mac sscase117, 320, 100000, x5cl_five_reals_last_unfft rsi, 5*64, 64

	sse2mac sscase118, 448, 4096, x7cl_seven_reals_first_fft rsi, 7*64, 64
	sse2mac sscase119, 448, 100000, x7cl_seven_reals_first_fft rsi, 7*64, 64
	sse2mac sscase120, 448, 4096, g7cl_seven_reals_first_fft rsi, 7*64, 64, rdx, 7*64, 64
	sse2mac sscase121, 448, 100000, g7cl_seven_reals_first_fft rsi, 7*64, 64, rdx, 7*64, 64
	sse2mac sscase122, 448, 4096, s7cl_seven_reals_first_fft rsi, 7*64, 64
	sse2mac sscase123, 448, 100000, s7cl_seven_reals_first_fft rsi, 7*64, 64
	sse2mac sscase124, 448, 4096, x7cl_seven_reals_last_unfft rsi, 7*64, 64
	sse2mac sscase125, 448, 100000, x7cl_seven_reals_last_unfft rsi, 7*64, 64

	sse2mac sscase126, 256, 4096, x4cl_empty rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase127, 256, 100000, x4cl_empty rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase128, 256, 4096, g4cl_empty rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase129, 256, 100000, g4cl_empty rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase130, 256, 4096, g4cl_empty_nt rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase131, 256, 100000, g4cl_empty_nt rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, 0

	sse2mac sscase132, 128, 4096, r4_x2cl_four_complex_first_fft4_scratch rsi, 2*64, 64, rdi
	sse2mac sscase133, 128, 100000, r4_x2cl_four_complex_first_fft4_scratch rsi, 2*64, 64, rdi
	sse2mac sscase134, 256, 4096, r4_x4cl_four_complex_last_unfft4 rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase135, 256, 100000, r4_x4cl_four_complex_last_unfft4 rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase136, 256, 4096, r4_x4cl_four_complex_fft_final rsi, 4*64, 64, 2*64
	sse2mac sscase137, 256, 100000, r4_x4cl_four_complex_fft_final rsi, 4*64, 64, 2*64
	sse2mac sscase138, 256, 4096, r4_x4cl_four_complex_with_square rsi, 4*64, 64, 2*64
	sse2mac sscase139, 256, 100000, r4_x4cl_four_complex_with_square rsi, 4*64, 64, 2*64
	sse2macbx sscase140, 256, 4096, r4_x4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	sse2macbx sscase141, 256, 100000, r4_x4cl_four_complex_with_mult rsi, 4*64, 64, 2*64
	sse2macbx sscase142, 256, 4096, r4_x4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64
	sse2macbx sscase143, 256, 100000, r4_x4cl_four_complex_with_mulf rsi, 4*64, 64, 2*64

	sse2mac sscase144, 128, 4096, r2_x2cl_two_complex_fft rsi, 2*64, 64, rdi
	sse2mac sscase145, 128, 100000, r2_x2cl_two_complex_fft rsi, 2*64, 64, rdi
	sse2mac sscase146, 128, 4096, r2_x2cl_two_complex_unfft rsi, 2*64, 64, rdi, 0
	sse2mac sscase147, 128, 100000, r2_x2cl_two_complex_unfft rsi, 2*64, 64, rdi, 0

	sse2mac sscase148, 128, 4096, r4_x2cl_four_complex_first_djbfft_scratch rsi, 2*64, 64, rdi, rdi
	sse2mac sscase149, 128, 100000, r4_x2cl_four_complex_first_djbfft_scratch rsi, 2*64, 64, rdi, rdi
	sse2mac sscase150, 256, 4096, r4_x4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase151, 256, 100000, r4_x4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase152, 256, 4096, r4_x4cl_wpn_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase153, 256, 100000, r4_x4cl_wpn_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi
	sse2mac sscase154, 256, 4096, r4_sg4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi
	sse2mac sscase155, 256, 100000, r4_sg4cl_four_complex_djbfft rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi
	sse2macbx sscase156, 128, 4096, r4_f2cl_four_complex_djbfft rsi, 2*64, 64, rdi
	sse2macbx sscase157, 128, 100000, r4_f2cl_four_complex_djbfft rsi, 2*64, 64, rdi

	sse2mac sscase158, 256, 4096, r4_x4cl_four_complex_last_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, 0
	sse2mac sscase159, 256, 100000, r4_x4cl_four_complex_last_djbunfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, 0
	sse2mac sscase160, 256, 4096, r4_x4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase161, 256, 100000, r4_x4cl_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase162, 256, 4096, r4_x4cl_wpn_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase163, 256, 100000, r4_x4cl_wpn_four_complex_djbunfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase164, 128, 4096, r4_sg2cl_four_complex_djbunfft rsi, 2*64, 64, rdx, 2*64, 64, rdi
	sse2mac sscase165, 128, 100000, r4_sg2cl_four_complex_djbunfft rsi, 2*64, 64, rdx, 2*64, 64, rdi

	sse2mac sscase166, 256, 4096, r4_sg4cl_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi
	sse2mac sscase167, 256, 100000, r4_sg4cl_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi
	sse2mac sscase168, 128, 4096, r4_sg2cl_four_complex_unfft4 rsi, 2*64, 64, rdx, 2*64, 64, rdi
	sse2mac sscase169, 128, 100000, r4_sg2cl_four_complex_unfft4 rsi, 2*64, 64, rdx, 2*64, 64, rdi

	sse2mac sscase170, 128, 4096, r4_x2cl_eight_reals_first_fft_scratch rsi, 2*64, 64, rdi
	sse2mac sscase171, 128, 100000, r4_x2cl_eight_reals_first_fft_scratch rsi, 2*64, 64, rdi
	sse2mac sscase172, 256, 4096, r4_x4cl_eight_reals_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, rdi
	sse2mac sscase173, 256, 100000, r4_x4cl_eight_reals_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, rdi
	sse2mac sscase174, 256, 4096, r4_x4cl_wpn_eight_reals_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, rdi
	sse2mac sscase175, 256, 100000, r4_x4cl_wpn_eight_reals_four_complex_djbfft rsi, 4*64, 64, 2*64, rdi, rdi
	sse2mac sscase176, 256, 4096, r4_sg4cl_eight_reals_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, rdi
	sse2mac sscase177, 256, 100000, r4_sg4cl_eight_reals_four_complex_fft4 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi, rdi

	sse2mac sscase178, 256, 4096, r4_x4cl_eight_reals_last_unfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase179, 256, 100000, r4_x4cl_eight_reals_last_unfft rsi, 4*64, 64, 2*64, rdi, 0
	sse2mac sscase180, 256, 4096, r4_x4cl_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, 0
	sse2mac sscase181, 256, 100000, r4_x4cl_eight_reals_unfft rsi, 4*64, 64, 2*64, rdi, 0, rdi, 0
	sse2mac sscase182, 128, 4096, r4_sg2cl_eight_reals_unfft4 rsi, 2*64, 64, rdx, 2*64, 64, rdi, rdi
	sse2mac sscase183, 128, 100000, r4_sg2cl_eight_reals_unfft4 rsi, 2*64, 64, rdx, 2*64, 64, rdi, rdi

	sse2mac sscase184, 192, 4096, r3_x3cl_six_reals_three_complex_djbfft rsi, 3*64, 64, rdi, rdi
	sse2mac sscase185, 192, 100000, r3_x3cl_six_reals_three_complex_djbfft rsi, 3*64, 64, rdi, rdi
	sse2mac sscase186, 192, 4096, r3_x3cl_six_reals_unfft rsi, 3*64, 64, rdi, 0, rdi, 0
	sse2mac sscase187, 192, 100000, r3_x3cl_six_reals_unfft rsi, 3*64, 64, rdi, 0, rdi, 0

	sse2mac sscase188, 192, 4096, r3_x3cl_three_complex_djbfft rsi, 3*64, 64, rdi
	sse2mac sscase189, 192, 100000, r3_x3cl_three_complex_djbfft rsi, 3*64, 64, rdi
	sse2mac sscase190, 192, 4096, r3_x3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, 0
	sse2mac sscase191, 192, 100000, r3_x3cl_three_complex_djbunfft rsi, 3*64, 64, rdi, 0

	sse2mac sscase192, 320, 4096, r5_x5cl_20_reals_first_fft_scratch rsi, 5*64, 64, rdi
	sse2mac sscase193, 320, 100000, r5_x5cl_20_reals_first_fft_scratch rsi, 5*64, 64, rdi
	sse2mac sscase194, 640, 4096, r5_x10cl_20_reals_last_unfft rsi, 10*64, 64, rdi, 8*XMM_SCD9
	sse2mac sscase195, 640, 100000, r5_x10cl_20_reals_last_unfft rsi, 10*64, 64, rdi, 8*XMM_SCD9

	sse2mac sscase196, 320, 4096, r5_x5cl_five_complex_djbfft rsi, 5*64, 64, rdi
	sse2mac sscase197, 320, 100000, r5_x5cl_five_complex_djbfft rsi, 5*64, 64, rdi
	sse2mac sscase198, 320, 4096, r5_x5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, 0
	sse2mac sscase199, 320, 100000, r5_x5cl_five_complex_djbunfft rsi, 5*64, 64, rdi, 0

	sse2mac sscase200, 448, 4096, r7_x7cl_28_reals_first_fft_scratch rsi, 7*64, 64, rdi
	sse2mac sscase201, 448, 100000, r7_x7cl_28_reals_first_fft_scratch rsi, 7*64, 64, rdi
	sse2mac sscase202, 896, 4096, r7_x14cl_28_reals_last_unfft rsi, 14*64, 64, rdi, 8*XMM_SCD13
	sse2mac sscase203, 896, 100000, r7_x14cl_28_reals_last_unfft rsi, 14*64, 64, rdi, 8*XMM_SCD13

	sse2mac sscase204, 512, 4096, r8_x8cl_eight_complex_fft_final rsi, 8*64, 64, 2*64, 4*64
	sse2mac sscase205, 512, 100000, r8_x8cl_eight_complex_fft_final rsi, 8*64, 64, 2*64, 4*64
	sse2mac sscase206, 512, 4096, r8_x8cl_eight_complex_with_square rsi, 8*64, 64, 2*64, 4*64
	sse2mac sscase207, 512, 100000, r8_x8cl_eight_complex_with_square rsi, 8*64, 64, 2*64, 4*64
	sse2macbx sscase208, 512, 4096, r8_x8cl_eight_complex_with_mult rsi, 8*64, 64, 2*64, 4*64
	sse2macbx sscase209, 512, 100000, r8_x8cl_eight_complex_with_mult rsi, 8*64, 64, 2*64, 4*64
	sse2macbx sscase210, 512, 4096, r8_x8cl_eight_complex_with_mulf rsi, 8*64, 64, 2*64, 4*64
	sse2macbx sscase211, 512, 100000, r8_x8cl_eight_complex_with_mulf rsi, 8*64, 64, 2*64, 4*64

	sse2mac sscase212, 512, 4096, r8_sg8cl_eight_complex_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi
	sse2mac sscase213, 512, 100000, r8_sg8cl_eight_complex_fft8 rsi, 8*64, 64, 2*64, 4*64, rdx, 8*64, 64, 2*64, 4*64, rdi
	sse2mac sscase214, 256, 4096, r8_sg4cl_eight_complex_unfft8 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi
	sse2mac sscase215, 256, 100000, r8_sg4cl_eight_complex_unfft8 rsi, 4*64, 64, 2*64, rdx, 4*64, 64, 2*64, rdi

	; need 16 reals cases

; Exit the timing code

exit:	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

gwtimeit ENDP

_TEXT ENDS
END

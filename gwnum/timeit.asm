; Copyright 1995-2009 Mersenne Research, Inc.  All rights reserved
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
INCLUDE xmult.mac
INCLUDE xlucas.mac

IFDEF X86_64
X87_CASES	EQU	0
ELSE
X87_CASES	EQU	13
ENDIF
SSE2_CASES	EQU	126

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
lab:	mov	rbx, 0			;; Offset for some sse2 macros (s.b. non-zero for with_mult macros)
	mov	rbp, 524288+256		;; Offset for mulf sse2 macros
	mov	eax, outer_iters
	mov	SRCARG, rdi		;; Save work buf addr
ss0a:	mov	rdi, SRCARG		;; Reload work buf addr (sincos data)
	lea	rsi, [rdi+4096]		;; Source & dest ptr
	lea	rdx, [rsi+524288+256]	;; Destination for "g" macros
	mov	ecx, inner_iters
ss0b:	&ops
	lea	rdi, [rdi+2*XMM_SCD]	;; Next sine/cosine pointer
	dec	ecx
	jnz	ss0b
	dec	eax
	jnz	ss0a
	jmp	exit
	ENDM
sse2macbx MACRO	lab, memused, memarea, ops:vararg
	LOCAL	ss0a, ss0b
	inner_iters = memarea / memused
	outer_iters = 10000 / inner_iters
lab:	mov	rbx, 262144+128		;; Offset for some sse2 macros (s.b. non-zero for with_mult macros)
	mov	rbp, 524288+256		;; Offset for mulf sse2 macros
	mov	eax, outer_iters
	mov	SRCARG, rdi		;; Save work buf addr
ss0a:	mov	rdi, SRCARG		;; Reload work buf addr (sincos data)
	lea	rsi, [rdi+4096]		;; Source & dest ptr
	lea	rdx, [rsi+524288+256]	;; Destination for "g" macros
	mov	ecx, inner_iters
ss0b:	&ops
	lea	rdi, [rdi+2*XMM_SCD]	;; Next sine/cosine pointer
	dec	ecx
	jnz	ss0b
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

PROCFLP	gwtimeit
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
; only, 96KB will read from the L2 cache only, and 2MB will test reading
; from main memory.

sscase0:
	read2	4096, 1000	; Read 4KB
	jmp	exit

sscase1:
	read2	96*1024, 100	; Read 96KB
	jmp	exit

sscase2:
	read2	2048*1024, 2	; Read 2MB
	jmp	exit

; This code writes a contiguous block of memory.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 2MB will test writing
; to main memory.

sscase3:
	write1	4096, 1000	; Read 4KB
	jmp	exit

sscase4:
	write1	96*1024, 100	; Read 96KB
	jmp	exit

sscase5:
	write1	2048*1024, 2	; Read 2MB
	jmp	exit

; This code writes a block of memory non-contiguously.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 2MB will test writing
; to main memory.

sscase6:
	write2	4096, 1000	; Read 4KB
	jmp	exit

sscase7:
	write2	96*1024, 100	; Read 96KB
	jmp	exit

sscase8:
	write2	2048*1024, 2	; Read 2MB
	jmp	exit

; This code reads & writes a block of memory.
; Timings are done on 3 memory size.  4KB will write to the L1 cache
; only, 96KB will write to L2 cache only, and 2MB will test writing
; to main memory.

sscase9:
	readwrite1	4096, 1000	; Read 4KB
	jmp	exit

sscase10:
	readwrite1	96*1024, 100	; Read 96KB
	jmp	exit

sscase11:
	readwrite1	2048*1024, 2	; Read 2MB
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
	sse2mac sscase34, 256, 4096, x4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64
	sse2mac sscase35, 256, 100000, x4cl_eight_reals_fft_2 rsi, 4*64, 64, 2*64
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

	sse2mac sscase50, 128, 4096, x2cl_two_complex_fft rsi, 2*64, 64
	sse2mac sscase51, 128, 100000, x2cl_two_complex_fft rsi, 2*64, 64
	sse2mac sscase52, 128, 4096, x2cl_two_complex_fft_in_place rsi, 2*64, 64
	sse2mac sscase53, 128, 100000, x2cl_two_complex_fft_in_place rsi, 2*64, 64

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
	sse2mac sscase72, 256, 4096, x4cl_four_complex_fft rsi, 4*64, 64, 2*64
	sse2mac sscase73, 256, 100000, x4cl_four_complex_fft rsi, 4*64, 64, 2*64
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
	sse2mac sscase88, 256, 4096, x4cl_four_complex_unfft rsi, 4*64, 64, 2*64
	sse2mac sscase89, 256, 100000, x4cl_four_complex_unfft rsi, 4*64, 64, 2*64
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

; Exit the timing code

exit:	ad_epilog 0,0,rbx,rbp,rsi,rdi,xmm6,xmm7,xmm8,xmm9,xmm10,xmm11,xmm12,xmm13,xmm14,xmm15

ENDPP	gwtimeit

_TEXT ENDS
END

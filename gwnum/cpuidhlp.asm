; Copyright 1995-2011 Mersenne Research, Inc.  All rights reserved
; Author:  George Woltman
; Email: woltman@alum.mit.edu
;
; This file implements helper routines for the CPU identification code
; as well as other helper functions
;

	TITLE   cpuidhlp

IFNDEF X86_64
	.686
	.XMM
	.MODEL	FLAT
ENDIF

INCLUDE unravel.mac

_TEXT SEGMENT

;
; Utility routine to initialize the FPU
;

; fpu_init ()
;	Initializes the x87 FPU.  Probably no longer necessary, but it was
;	necessary a long, long, time ago.
; Windows 32-bit (_fpu_init)
; Linux 32-bit (fpu_init)
;	No parameters or return value
; Windows 64-bit (fpu_init) - leaf routine, no unwind info necessary
; Linux 64-bit (fpu_init)
;	No parameters or return value

PROCL	fpu_init
IFNDEF X86_64
	fninit
ENDIF
	ret
fpu_init ENDP

;
; Utility routine to read the time stamp counter
;

; erdtsc (&hi, &lo)
;	Read the 64-bit timestamp counter and return it in 32-bit chunks
; Windows 32-bit (_erdtsc)
; Linux 32-bit (erdtsc)
;	Parameter hi = [esp+4]
;	Parameter lo = [esp+8]
; Windows 64-bit (erdtsc) - leaf routine, no unwind info necessary
;	Parameter hi = rcx
;	Parameter lo = rdx
; Linux 64-bit (erdtsc)
;	Parameter hi = rdi
;	Parameter lo = rsi

PROCL	erdtsc
IFNDEF X86_64
	rdtsc
	mov	ecx, [esp+4]		; Address to return info
	mov	[rcx], edx		; store high 32 bits
	mov	ecx, [esp+8]		; Address to return info
	mov	[rcx], eax		; store low 32 bits
ENDIF
IFDEF WINDOWS64
	mov	r8, rdx			; Copy return info address
	rdtsc
	mov	[rcx], edx		; store high 32 bits
	mov	[r8], eax		; store low 32 bits
ENDIF
IFDEF LINUX64
	rdtsc
	mov	[rdi], edx		; store high 32 bits
	mov	[rsi], eax		; store low 32 bits
ENDIF
	ret
erdtsc	ENDP


;
; Utility routine to see if CPUID is supported
; Returns non-zero if CPUID is supported
;

; ecpuidsupport ()
;	Return true if CPU supports CPUID
; Windows 32-bit (_ecpuidsupport)
; Linux 32-bit (ecpuidsupport)
;	Return value = eax
; Windows 64-bit (ecpuidsupport) - leaf routine, no unwind info necessary
; Linux 64-bit (ecpuidsupport)
;	Return value = rax

PROCL	ecpuidsupport
IFNDEF X86_64
        pushfd				; Get original EFLAGS
	pop	eax
	mov 	ecx, eax
        xor     eax, 200000h		; Flip ID bit in EFLAGS
        push    eax			; Save new EFLAGS value on stack
        popfd				; Replace current EFLAGS value
        pushfd				; Get new EFLAGS
        pop     eax			; Store new EFLAGS in EAX
        xor     eax, ecx		; Test toggled ID bit
ELSE
	mov	rax, 1			; 64-bit CPUs support CPUID
ENDIF
	ret				; Processor supports CPUID if eax != 0
ecpuidsupport ENDP

;
; Utility routine to execute CPUID
;

; ecpuid (resptr)
;	Call cpuid a fill a structure with the eax, ebx, ecx, edx values
; Windows 32-bit (_ecpuid)
; Linux 32-bit (ecpuid)
;	Parameter resptr = [esp+4]
; Windows 64-bit (ecpuid)
;	Parameter resptr = rcx
; Linux 64-bit (ecpuid)
;	Parameter res = rdi

PROCFL	ecpuid
	ah_prolog 1,0,0,rbx,rdi
IFNDEF X86_64
	mov	rdi, [esp+push_amt+4]	; Address of data struct to return info
ENDIF
IFDEF WINDOWS64
	mov	rdi, rcx		; Address of data struct to return info
ENDIF
	mov	eax, [rdi]		; EAX argument to CPUID
	mov	ecx, [rdi+8]		; ECX argument to CPUID
	cpuid				; Perform the CPUID
	mov	[rdi], eax		; Return the results
	mov	[rdi+4], ebx		; Return the results
	mov	[rdi+8], ecx		; Return the results
	mov	[rdi+12], edx		; Return the results
	ah_epilog 1,0,0,rbx,rdi
ecpuid	ENDP

;
; Utility routine to execute XGETBV
;

; exgetbv (resptr)
;	Call xgetbx a fill a structure with the eax, ebx, ecx, edx values
; Windows 32-bit (_ecpuid)
; Linux 32-bit (ecpuid)
;	Parameter resptr = [esp+4]
; Windows 64-bit (ecpuid)
;	Parameter resptr = rcx
; Linux 64-bit (ecpuid)
;	Parameter res = rdi

PROCFL	exgetbv
	ah_prolog 1,0,0,rbx,rdi
IFNDEF X86_64
	mov	rdi, [esp+push_amt+4]	; Address of data struct to return info
ENDIF
IFDEF WINDOWS64
	mov	rdi, rcx		; Address of data struct to return info
ENDIF
	mov	ecx, [rdi+8]		; ECX argument to XGETBV
	xgetbv				; Perform the XGETBV
	mov	[rdi], eax		; Return the results
	mov	[rdi+4], ebx		; Return the results
	mov	[rdi+8], ecx		; Return the results
	mov	[rdi+12], edx		; Return the results
	ah_epilog 1,0,0,rbx,rdi
exgetbv	ENDP

;
; Utility routines that hopefully take (roughly) 100,000 and 1,000,000 clock cycles
; This implementation is based on the ROR op-code should has a reciprocal throughput
; of 1 on Intel architectures according to Agner Fog's documents.  This routine will
; be inaccurate on AMD cpus which have have a reciprocal thoughput of 1/3.
;
; NOTE:  The ROR op-code inexplicably takes fewer clock cycles than expected
; on a Pentium 4 (97600 instead of 100000).  I tried using LEA instead:
;	lea	rbx, [rbx+rbx*8]
;	lea	rcx, [rcx+rcx*8]
;	lea	rdx, [rdx+rdx*8]
;	lea	rdi, [rdi+rdi*8]
;	lea	rsi, [rsi+rsi*8]
; but this also ran too fast (83000 instead of 100000).

PROCFL	one_hundred_thousand_clocks
	ah_prolog 0,0,0,rdi,rsi
	mov	rax, 100000 / 100
biglp:	REPEAT 25
	ror	ecx, 1
	ror	edx, 1
	ror	edi, 1
	ror	esi, 1
	ENDM
	sub	rax, 1
	jnz	biglp
	ah_epilog 0,0,0,rdi,rsi
one_hundred_thousand_clocks ENDP

PROCFL	one_million_clocks
	ah_prolog 0,0,0,rdi,rsi
	mov	rax, 1000000 / 100
biglp2:	REPEAT 25
	ror	ecx, 1
	ror	edx, 1
	ror	edi, 1
	ror	esi, 1
	ENDM
	sub	rax, 1
	jnz	biglp2
	ah_epilog 0,0,0,rdi,rsi
one_million_clocks ENDP

; This was our original implementation.  Unfortunately it took 200,000 and 2,000,000
; clocks on a Pentium 4 (and 1.5 times that amount on a hyperthreaded Pentium 4).

IFDEF ORIGNAL_IMPLEMENTATION

PROCFL	one_hundred_thousand_clocks
	ah_prolog 0,0,0,xmm6,xmm7
	xorps	xmm0, xmm0
	xorps	xmm1, xmm1
	xorps	xmm2, xmm2
	xorps	xmm3, xmm3
	xorps	xmm4, xmm4
	xorps	xmm5, xmm5
	xorps	xmm6, xmm6
	xorps	xmm7, xmm7
	mov	rax, 100000 / 200
biglp:	REPEAT 25
	addps	xmm0, xmm0
	addps	xmm1, xmm1
	addps	xmm2, xmm2
	addps	xmm3, xmm3
	addps	xmm4, xmm4
	addps	xmm5, xmm5
	addps	xmm6, xmm6
	addps	xmm7, xmm7
	ENDM
	sub	rax, 1
	jnz	biglp
	ah_epilog 0,0,0,xmm6,xmm7
one_hundred_thousand_clocks ENDP

PROCFL	one_million_clocks
	ah_prolog 0,0,0,xmm6,xmm7
	xorps	xmm0, xmm0
	xorps	xmm1, xmm1
	xorps	xmm2, xmm2
	xorps	xmm3, xmm3
	xorps	xmm4, xmm4
	xorps	xmm5, xmm5
	xorps	xmm6, xmm6
	xorps	xmm7, xmm7
	mov	rax, 1000000 / 200
biglp2:	REPEAT 25
	addps	xmm0, xmm0
	addps	xmm1, xmm1
	addps	xmm2, xmm2
	addps	xmm3, xmm3
	addps	xmm4, xmm4
	addps	xmm5, xmm5
	addps	xmm6, xmm6
	addps	xmm7, xmm7
	ENDM
	sub	rax, 1
	jnz	biglp2
	ah_epilog 0,0,0,xmm6,xmm7
one_million_clocks ENDP

ENDIF

_TEXT	ENDS
END

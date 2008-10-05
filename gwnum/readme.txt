Welcome to the gwnum library.

-> how to use
	Before calling gwinit and one of the gwsetup routines, the
	global variables defined in cpuid.h must be initialized.  You can
	do this yourself or let gwinit do it for you (it will call the
	guessCpuType and guessCpuSpeed routines).

	Next, study gwnum.h for a list of functions.  Also, study prp.c
	or commonb.c or ecm.c for examples of how to use the gwnum library.

-> how to compile / how to port

Windows:
	To make the libraries, your path must already be set to the proper
	C compiler (32-bit or 64-bit).  Then execute
		nmake /f compile
	to build the 32-bit libraries and 
		nmake /f compil64
	to build the 64-bit libraries.

Linux:
	The object files have already been converted to ELF format using
	Agner Fog's objconv on a Windows machine.  Copy gwnum.a from the
	gwnum\linux directory to gwnum.

	makefile is used to compile the necessary C and C++ source files on the
	Linux machine and to finish building the gwnum library.  It adds to the
	ELF library built on the Windows machine.
	
	A C application using gwnum must link with:
		gwnum.a gwnum.ld -lpthread -lstdc++

Linux 64-bit:
	The object files have already been converted to ELF format using
	Agner Fog's objconv on a Windows machine.  Copy gwnum.a from the
	gwnum\linux64 directory to gwnum.

	make64 is used to compile the necessary C and C++ source files on the
	Linux machine and to finish building the gwnum library
	
	A C application using gwnum must link with:
		gwnum.a gwnum.ld -lpthread -lstdc++

FreeBSD:
	The object files have already been converted to ELF format using
	Agner Fog's objconv on a Windows machine.  Copy gwnum.a from the
	gwnum\linux directory to gwnum.

	make.bsd is used to compile the	necessary C and C++ source files and
	finish building the gwnum library.  It adds to the ELF library
	built on the Windows machine.

	A C application must link with:
		-lcompat gwnum.a gwnum.ld -lpthread -lstdc++

Mac OS X:
	The object files have already been converted to Mach-O format using
	Agner Fog's objconv on a Windows machine.  Copy gwnum.a from the
	gwnum\macosx directory to gwnum.

	makemac is used to compile the necessary C and C++ source files and
	finish building the gwnum library.  It adds to the MAC OS X library
	built on the Windows machine.

	A C application must link with gwnum.a -lpthread -lstdc++

-> Legal stuff

Copyright (c) 1996-2008, Mersenne Research, Inc.  All rights reserved. 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are 
met: 

(1) Redistributing source code must contain this copyright notice,
limitations, and disclaimer. 
(2) If this software is used to find Mersenne Prime numbers, then
GIMPS will be considered the discoverer of any prime numbers found
and the prize rules at http://mersenne.org/prize.htm will apply.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
  


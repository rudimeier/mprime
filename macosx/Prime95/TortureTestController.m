//
//  TortureTestController.m
//  Prime95
//
//  Created by George Woltman on 4/25/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import "TortureTestController.h"
#include "prime95.h"

@implementation TortureTestController

- (id)init
{
	if (![super initWithWindowNibName:@"TortureTest"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	int	mem, in_place_fft;

	[self setTortureType:2];

	[self setNumberOfThreads:NUM_CPUS * CPU_HYPERTHREADS];
	[self setNumberOfThreadsMax:NUM_CPUS * CPU_HYPERTHREADS];
	[self setNumberOfThreadsEnabled:(NUM_CPUS * CPU_HYPERTHREADS > 1)];
	[self setCustomSettingsEnabled:NO];
	[self setCustomMemoryEnabled:NO];
	[self setMinFFTSize:8];
	[self setMaxFFTSize:4096];

	mem = physical_memory ();
	if (mem >= 2000) {
		blendMemory = GetSuggestedMemory (1600);
		in_place_fft = FALSE;
	} else if (mem >= 500) {
		blendMemory = GetSuggestedMemory (mem - 256);
		in_place_fft = FALSE;
	} else if (mem >= 200) {
		blendMemory = GetSuggestedMemory (mem / 2);
		in_place_fft = TRUE;
	} else {
		blendMemory = 8;
		in_place_fft = TRUE;
	}
	[self setRunFFTsInPlace:in_place_fft];
	[self setMemoryToUse:blendMemory];
	[self setTimeToRunEachFFT:15];
}

- (int)tortureType
{
	return tortureType;
}

- (void)setTortureType:(int) _value
{
	if (_value == 0) {			// Small FFTs (in L2 cache)
		int	maxfft;
		[self setCustomSettingsEnabled:NO];
		[self setCustomMemoryEnabled:NO];
		[self setMinFFTSize:8];
		if (CPU_L2_CACHE_SIZE <= 128) maxfft = 8;
		else if (CPU_L2_CACHE_SIZE <= 256) maxfft = 16;
		else if (CPU_L2_CACHE_SIZE <= 512) maxfft = 32;
		else maxfft = 64;
		[self setMaxFFTSize:maxfft];
		[self setRunFFTsInPlace:YES];
		[self setMemoryToUse:0];
		[self setTimeToRunEachFFT:15];
	} else if (_value == 1) {		// Large in-place FFTs
		[self setCustomSettingsEnabled:NO];
		[self setCustomMemoryEnabled:NO];
		[self setMinFFTSize:8];
		[self setMaxFFTSize:1024];
		[self setRunFFTsInPlace:YES];
		[self setMemoryToUse:8];
		[self setTimeToRunEachFFT:15];
	} else if (_value == 2) {		// Blend
		[self setCustomSettingsEnabled:NO];
		[self setCustomMemoryEnabled:NO];
		[self setMinFFTSize:8];
		[self setMaxFFTSize:4096];
		[self setRunFFTsInPlace:NO];
		[self setMemoryToUse:blendMemory];
		[self setTimeToRunEachFFT:15];
	} else {				// Custom
		[self setCustomSettingsEnabled:YES];
		[self setCustomMemoryEnabled:!runFFTsInPlace];
	}

	tortureType = _value;
}

- (int)runFFTsInPlace
{
	return runFFTsInPlace;
}

- (void)setRunFFTsInPlace:(int) _value
{
	[self setCustomMemoryEnabled:(customMemoryEnabled && !_value)];

	runFFTsInPlace = _value;
}

@synthesize numberOfThreads;
@synthesize numberOfThreadsMax;
@synthesize numberOfThreadsEnabled;
@synthesize customSettingsEnabled;
@synthesize customMemoryEnabled;
@synthesize minFFTSize;
@synthesize maxFFTSize;
@synthesize memoryToUse;
@synthesize timeToRunEachFFT;

- (IBAction)ok:(id)sender
{
	int	mem;

	IniWriteInt (INI_FILE, "MinTortureFFT", minFFTSize);
	IniWriteInt (INI_FILE, "MaxTortureFFT", maxFFTSize);
	mem = memoryToUse / numberOfThreads;
	if (runFFTsInPlace) mem = 8;
	IniWriteInt (INI_FILE, "TortureMem", mem);
	IniWriteInt (INI_FILE, "TortureTime", timeToRunEachFFT);
	LaunchTortureTest (numberOfThreads, FALSE);

	[[self window] performClose:self];
}

@end

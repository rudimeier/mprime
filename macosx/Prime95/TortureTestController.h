//
//  TortureTestController.h
//  Prime95
//
//  Created by George Woltman on 4/25/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface TortureTestController : NSWindowController {
	int	tortureType;
	int	numberOfThreads;
	int	numberOfThreadsMax;
	int	numberOfThreadsEnabled;
	int	customSettingsEnabled;
	int	customMemoryEnabled;
	int	minFFTSize;
	int	maxFFTSize;
	int	runFFTsInPlace;
	int	memoryToUse;
	int	timeToRunEachFFT;
	int	blendMemory;
}

@property(readwrite, assign) int tortureType;
@property(readwrite, assign) int numberOfThreads;
@property(readwrite, assign) int numberOfThreadsMax;
@property(readwrite, assign) int numberOfThreadsEnabled;
@property(readwrite, assign) int customSettingsEnabled;
@property(readwrite, assign) int customMemoryEnabled;
@property(readwrite, assign) int minFFTSize;
@property(readwrite, assign) int maxFFTSize;
@property(readwrite, assign) int runFFTsInPlace;
@property(readwrite, assign) int memoryToUse;
@property(readwrite, assign) int timeToRunEachFFT;

- (void)reInit;
- (IBAction)ok:(id)sender;

@end

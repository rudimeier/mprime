//
//  TimeController.m
//  Prime95
//
//  Created by George Woltman on 4/24/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import "TimeController.h"
#include "prime95.h"

@implementation TimeController

- (id)init
{
	if (![super initWithWindowNibName:@"Time"]) return nil;
	exponentToTime = 50000000;
	exponentToTimeMax = MAX_PRIME_SSE2;
	numberOfIterations = 50;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
}

@synthesize exponentToTime;
@synthesize exponentToTimeMax;
@synthesize numberOfIterations;

- (IBAction)ok:(id)sender
{
	LaunchAdvancedTime (exponentToTime, numberOfIterations);
	[[self window] performClose:self];
}

@end

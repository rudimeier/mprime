//
//  StopController.m
//  Prime95
//
//  Created by George Woltman on 4/19/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import "StopController.h"
#include "prime95.h"

@implementation StopController

- (id)init
{
	if (![super initWithWindowNibName:@"Stop"]) return nil;
	return self;
}

- (void)windowDidLoad
{
	[self reInit];
}

- (void)reInit
{
	[self setStopAllWorkers:YES];
	[self setWorkerNumber:1];
	[self setWorkerNumberMax:NUM_WORKER_THREADS];
}

@synthesize stopAllWorkers;
@synthesize workerNumber;
@synthesize workerNumberMax;

- (IBAction)ok:(id)sender
{
	if (stopAllWorkers)
		stop_workers_for_escape ();
	else
		stop_one_worker (workerNumber - 1);
	[[self window] performClose:self];
}

@end

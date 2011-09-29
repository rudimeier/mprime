//
//  WorkerWindowsController.h
//  Prime95
//
//  Created by George Woltman on 4/26/09.
//  Copyright 2009 Mersenne Research, Inc. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface WorkerWindowsController : NSWindowController {
	int	numWorkers;
	int	numWorkersMax;
	int	numWorkersEnabled;
	int	priority;
	NSMutableArray *workerData;
	NSMutableArray *affinityValues;
	IBOutlet NSArrayController *workerDataArrayController;
	int	startupInProgress;
}

@property(readwrite, assign) int numWorkers;
@property(readwrite, assign) int numWorkersMax;
@property(readwrite, assign) int numWorkersEnabled;
@property(readwrite, assign) int priority;
@property(copy) NSArray *workerData;
@property(copy) NSArray *affinityValues;

- (void)reInit;
- (IBAction)ok:(id)sender;

@end


@interface WorkerData : NSObject
{
	NSString *workerNumber;
	int	typeOfWork;
	int	affinity;
	int	affinityEnabled;
	int	multithreading;
	int	multithreadingMax;
	int	multithreadingEnabled;
}

@property(copy) NSString *workerNumber;
@property(readwrite, assign) int typeOfWork;
@property(readwrite, assign) int affinity;
@property(readwrite, assign) int affinityEnabled;
@property(readwrite, assign) int multithreading;
@property(readwrite, assign) int multithreadingMax;
@property(readwrite, assign) int multithreadingEnabled;

@end

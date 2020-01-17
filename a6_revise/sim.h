//
//  sim.h
//  hwk6
//
//  Created by REAL CYH on 11/10/19.
//  Copyright © 2019 REAL CYH. All rights reserved.
//

#ifndef sim_h
#define sim_h

#include <stdio.h>

void RunSim (double EndTime, char* outfile);

// Schedule an event with timestamp ts, event parameters *data
void Schedule (double ts, void *data);

// This function returns the current simulation time
double CurrentTime (void);

//
// Function defined in the simulation application called by the simulation engine
//
//  Event handler function: called to process an event
void EventHandler (void *data, char* outfile);

#endif /* engine_h */

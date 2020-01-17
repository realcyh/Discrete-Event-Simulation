//
//  engine.h
//  hwk6
//
//  Created by REAL CYH on 11/10/19.
//  Copyright © 2019 REAL CYH. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "sim.h"

struct Event {
    double timestamp;        // event timestamp
    void *AppData;            // pointer to application defined event parameters
    struct Event *Next;        // priority queue pointer
};

// Simulation clock variable
double Now = 0.0;

// Future Event List
// Use an event structure as the header for the future event list (priority queue)
// Using an event struct as the header for the priority queue simplifies the code for
// inserting/removing events by eliminating the need to explicitly code special cases
// such as inserting into an empty priority queue, or removing the last event in the priority queue.
// See the Remove() and Schedule() functions below.
struct Event FEL ={-1.0, NULL, NULL};

/////////////////////////////////////////////////////////////////////////////////////////////
// Prototypes for functions used within the Simulation Engine
/////////////////////////////////////////////////////////////////////////////////////////////

// Function to print timestamps of events in event list
void PrintList (char* outfile);

// Function to remove smallest timestamped event
struct Event *Remove (void);


/////////////////////////////////////////////////////////////////////////////////////////////
//
// General Purpose Discrete Event Simulation Engine
//
/////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////
// Simulation Engine Data Structures
/////////////////////////////////////////////////////////////////////////////////////////////
//
// Data srtucture for an event; note this is designed to be independent of the application domain.
// Each event can have parameters defined within the spplication. We want the simulation engine
// not to have to know the number or types of these parameters, since they are dependent on the
// application, and we want To keep the engine independent of the application. To address this problem
// each event only contains a single parameter, a pointer to the application defined parameters (AppData).
// The simulation engine only knows it has a pointer to the event parameters, but does not know the
// structure (number of parameters or their type) of the information pointed to.
// This way the event can have application-defined information, but the simulation engine need not
// know the number or type of the application-defined parameters.
//


// Remove smallest timestamped event from FEL, return pointer to this event
// return NULL if FEL is empty
struct Event* Remove(void)
{
    struct Event* e;

    if (FEL.Next == NULL) return (NULL);
    e = FEL.Next;        // remove first event in list
    FEL.Next = e->Next;
    return (e);
}

// Print timestamps of all events in the event list (used for debugging)
void PrintList(char* outfile)
{
    FILE* fwrite;
    if ((fwrite = fopen(outfile, "a")) == NULL) {
        fprintf(stderr, "cannot open output file\n"); exit(1);
    }
    
    struct Event* p;

    fprintf(fwrite, "\nEvent List: ");
    for (p = FEL.Next; p != NULL; p = p->Next) {
        fprintf(fwrite,"%f ", p->timestamp);
    }
    fprintf(fwrite, "\n");
    fclose(fwrite);
}

// Return current simulation time
double CurrentTime(void)
{
    return (Now);
}

// Schedule new event in FEL
// queue is implemented as a timestamp ordered linear list
void Schedule(double ts, void* data)
{
    struct Event* e, * track, * dummy;

    // create event data structure and fill it in
    if ((e = malloc(sizeof(struct Event))) == NULL) exit(1);
    e->timestamp = ts;
    e->AppData = data;

    // insert into priority queue
    for (dummy = &FEL, track = FEL.Next; track != NULL; track = track->Next, dummy = dummy->Next) {
        if (track->timestamp >= e->timestamp) break;
    }
    // insert after q (before p)
    e->Next = dummy->Next;
    dummy->Next = e;
}

// Function to execute simulation up to a specified time (EndTime)
void RunSim(double EndTime, char* outfile)
{
    FILE* fwrite;
    if ((fwrite = fopen(outfile, "a")) == NULL) {
        fprintf(stderr, "cannot open output file\n"); exit(1);
    }
    struct Event* e;

    fprintf(fwrite, "Initial event list:\n");
    fclose(fwrite);
    PrintList(outfile);
    

    // Main scheduler loop
    while ((e = Remove()) != NULL) {
        Now = e->timestamp;
		if (Now > EndTime) {
			free(e);
			while (FEL.Next) {
				struct Event* p = FEL.Next;
				FEL.Next = (FEL.Next)->Next;
				free(p);
			}
			break;

		}
        EventHandler(e->AppData, outfile);
        free(e);    // it is up to the event handler to free memory for parameters
        PrintList(outfile);
    }
}


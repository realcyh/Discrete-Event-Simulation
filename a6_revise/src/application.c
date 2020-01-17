//
//  application.h
//  hwk6
//
//  Created by REAL CYH on 11/10/19.
//  Copyright © 2019 REAL CYH. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "sim.h"

#define     GENERATE    0
#define     ARRIVAL     1
#define     DEPARTURE   2

#define     MAXCOMPONENTS    100
#define     GENERATE    0

#define     GENERATOR       0
#define     QUEUE_STATION   1
#define     EXIT            2


//statistics for calculation;
int inCust = 0;//number of customers entered this system
int outCust = 0;//number fo customers exited this system
//The minimum, maximum, and average amount of time a customer remains in the system among those customers who exited during the simulation run.---------------
double minTime = 1000.0;//update in arrival events in exit;
double maxTime = 0.0;
double sumTime = 0.0;
//-----------------------
//The minimum, maximum, and average amount of the total time each customer must wait in queues in traveling through the system.--------------
double minTimeinQ = 1000.0;// update in arrival events in exit;
double maxTimeinQ = 0.0;
double sumTimeinQ = 0.0;

double timeLimit=0;

typedef struct Customer {
    double creationTime;        // time customer was created
    double exitTime;            // time customer departs system
    double enterStaTime;         // time customer enter a QS
    double totalWaitTimeinQ;   // total time customer wait in QS;
    struct Customer* next;      // pointer to next customer when it is queued in a FIFO
}Customer;

struct Component {
    int ComponentType;  // GENERATOR, QUEUE_STATION, EXIT
    void* ptrComp;         // Pointer to information on component (Generator, Exit struct, etc.)
} arrComps[MAXCOMPONENTS];

typedef struct EventData {
    int eventType;                // Type of event (GENERATE, ARRIVAL, DEPARTURE)
    struct Customer* cust;        // Arriving or departing customer; unused for GENERATE events
    int compID;                    // ID of component where customer created, arrives, or departs
}EventData;

typedef struct Generator {//struct for component G
    double geneTime;     // mean interarrival time for generated components
    int destComp;          // ID of next component customers are sent to
}Generator;

//Customers waitting in queue struct for Station component:
typedef struct FIFOQueue {
    struct Customer* first;     // pointer to first customer in queue
    struct Customer* last;      // pointer to last customer in queue
    double servetime;           // average serve time of a QS
    int* destComp;              // ID of next component customers are sent to
    double totalWaitTime;        // total waiting time
    int numberofServedCust;
    double* destPossi;
    int routingNum;
}FIFOQueue;


typedef struct Exit {        // struct for component E
    int Count;              // number of customers that exited at this component
}Exit;


double urand(void);
double randexp(double U);
double expTime(double avgTime);

void makeGenerator(int GenID, double interTime, int DestID, char* outfile);
void makeSta(int staID, double servTime, int destID[], double possibility[], int routNum, char* outfile);
void makeExit(int ExitID, char* outfile);

void Generate(struct EventData* e, char* outfile);
void Arrival(struct EventData* e, char* outfile);
void Departure(struct EventData* e, char* outfile);

void configFinished(int numofComps, char* outfile);
int inputConfig(char* config, char* outfile);

int getInCust(void);
int getOutCust(void);
double getMinTime(void);
double getMaxTime(void);
double getAvgTime(void);
double getMinTimeinQ(void);
double getMaxTimeinQ(void);
double getAvgTimeinQ(void);
double avgTimeofSta(int CompID);

double urand() {
    return ((double)rand() / (RAND_MAX));
}

double randexp(double U) {
    return   (-1) * U * (log(1.0 - urand()));
}

double expTime(double avgTime) //To generate random numbers from an exponential distribution with mean U,
{
    return randexp(avgTime);
}


void EventHandler(void* data, char* outfile)
{
    struct EventData* d;

    // coerce type so the compiler knows the type of information pointed to by the parameter data.
    d = (struct EventData*) data;
    // call an event handler based on the type of event
    if (d->eventType == GENERATE) Generate(d, outfile);
    else if (d->eventType == ARRIVAL) Arrival(d, outfile);
    else if (d->eventType == DEPARTURE) Departure(d, outfile);
    else { fprintf(stderr, "Illegal event found\n"); exit(1); }
    free(d); // Release memory for event paramters
}

void makeGenerator(int GenID, double interTime, int DestID, char* outfile)
{
    FILE* fwrite;
    if ((fwrite = fopen(outfile, "w")) == NULL) {fprintf(stderr, "cannot open output file\n"); exit(1);}
    Generator* ptrGene;
    fprintf(fwrite, "Creating Generator Component, ID=%d, Interarrival time=%f, Destination=%d\n", GenID, interTime, DestID);

    // Allocate space for component, fill in parameters
    if ((ptrGene = malloc(sizeof(struct Generator))) == NULL) {
        fprintf(stderr, "malloc error\n");
        exit(1);
    }
    ptrGene->geneTime = interTime;
    ptrGene->destComp = DestID;
    arrComps[GenID].ComponentType = GENERATOR;
    arrComps[GenID].ptrComp = ptrGene;

    // schedule initial, first generator event for this component
    EventData* geneEvent;
    if ((geneEvent = malloc(sizeof(struct EventData))) == NULL) {
        fprintf(stderr, "malloc error\n");
        exit(1);
    }
    geneEvent->eventType = GENERATE;
    geneEvent->cust = NULL;
    geneEvent->compID = GenID;
    Schedule(0, geneEvent);
    //Schedule(expTime(interTime), geneEvent);
    fclose(fwrite);
}


// Creat an Queue Station Component with its ID and average serve time. The station
// takes in customer into the queue and send head customer to component DestID
void makeSta(int staID, double servTime, int destID[], double possibility[], int routNum, char* outfile) {

    FILE* fwrite;
    if ((fwrite = fopen(outfile, "a")) == NULL) {
        fprintf(stderr, "cannot open output file\n"); exit(1);
    }
    fprintf(fwrite, "Creating Queue Station Component, ID=%d, Average serve time=%f\n", staID, servTime);

    //Allocate space for component, fill in parameters;
    FIFOQueue* queue;
    if ((queue = malloc(sizeof(struct FIFOQueue))) == NULL) {
        fprintf(stderr, "malloc error\n");
        exit(1);
    }
    queue->first = NULL;
    queue->last = NULL;
    queue->servetime = servTime;
    queue->destComp = destID;
    queue->totalWaitTime = 0.0;
    queue->destPossi = possibility;
    queue->numberofServedCust = 0;//for static 4.
    queue->routingNum = routNum;
    //Add component to component array;
    arrComps[staID].ComponentType = QUEUE_STATION;
    arrComps[staID].ptrComp = queue;
    fclose(fwrite);
}

// Create an Exit Component with identifier ExitID
void makeExit(int ExitID, char* outfile)
{
    FILE* fwrite;
    if ((fwrite = fopen(outfile, "a")) == NULL) {
        fprintf(stderr, "cannot open output file\n"); exit(1);
    }
    Exit* p;

    fprintf(fwrite, "Creating Exit Component, ID=%d\n", ExitID);
    arrComps[ExitID].ComponentType = EXIT;

    // Allocate space for component, fill in parameters
    if ((p = malloc(sizeof(struct Exit))) == NULL) {
        fprintf(stderr, "malloc error\n");
        exit(1);
    }
    p->Count = 0;
    arrComps[ExitID].ptrComp = p;
    fclose(fwrite);
}

void Generate(struct EventData* e, char* outfile)
{
    FILE* fwrite;
    if ((fwrite = fopen(outfile, "a")) == NULL) {
        fprintf(stderr, "cannot open output file\n"); exit(1);
    }

    EventData* geneEvent;
    Customer* newCust;
    double curTime;
    if (e->eventType != GENERATE) {
        fprintf(stderr, "Unexpected event type\n");
        exit(1);
    }

    if (arrComps[e->compID].ComponentType != GENERATOR) {
        fprintf(stderr, "bad componenet type\n");
        exit(1);
    }

    //Get the corresponding Generator of this event
    Generator* ptrGene = (Generator*)arrComps[e->compID].ptrComp;

    // Create a new customer
    if ((newCust = malloc(sizeof(struct Customer))) == NULL) {
        fprintf(stderr, "malloc error\n");
        exit(1);
    }

    inCust++; //a costumer goes into the system

    //Initiate the customer
    newCust->creationTime = CurrentTime();
    newCust->exitTime = 0;
    newCust->enterStaTime = 0;
    newCust->totalWaitTimeinQ = 0;
    newCust->next = NULL;

    fprintf(fwrite, "Generate new customer at %f\n", CurrentTime());

    //Schedule arrival event at component connected to generator
    if ((geneEvent = malloc(sizeof(struct EventData))) == NULL) {
        fprintf(stderr, "malloc error\n");
        exit(1);
    }
    //Put the follow ARRIVAL event after the generate into log
    geneEvent->eventType = ARRIVAL;
    geneEvent->cust = newCust;
    geneEvent->compID = ptrGene->destComp;
    curTime = CurrentTime();
    Schedule(curTime, geneEvent);

    // Schedule next generation event

    double interTime = expTime(ptrGene->geneTime);
    //double interTime = ptrGene->geneTime;       

    curTime = CurrentTime() + interTime;       // After the interarrival time, next customer is generated

    if ((geneEvent = malloc(sizeof(struct EventData))) == NULL) {
        fprintf(stderr, "malloc error\n");
        exit(1);
    }
    geneEvent->eventType = GENERATE;
    geneEvent->compID = e->compID; //a generate event happens
    Schedule(curTime, geneEvent);
    if(curTime > timeLimit) free(geneEvent);
    fprintf(fwrite, "schedule next generator at time %f\n", curTime);
    fclose(fwrite);
}

void Arrival(struct EventData* e, char* outfile) {
    FILE* fwrite;
    if ((fwrite = fopen(outfile, "a")) == NULL) {
        fprintf(stderr, "cannot open output file\n"); exit(1);
    }
    EventData* geneEvent;
    Exit* ptrExit;
    FIFOQueue* ptrFIFO;
    double ts;

    if (e->eventType != ARRIVAL) {
        fprintf(stderr, "Unexpected event type\n");
        exit(1);
    }

    if (arrComps[e->compID].ComponentType == GENERATOR) {
        fprintf(stderr, "Wrong function\n");
        exit(1);
    }

    if (arrComps[e->compID].ComponentType == EXIT) {
        fprintf(fwrite, "Processing arrival event in Exit, component ID is %d, current time is %f\n", e->compID, CurrentTime());
        ptrExit = (Exit*)arrComps[e->compID].ptrComp;
        (ptrExit->Count)++;
        outCust++;
        e->cust->exitTime = CurrentTime();      //Final exiting time.

        double waitTime = e->cust->exitTime - e->cust->creationTime;    //Time that the customer stays in the system
        if (waitTime > maxTime) { maxTime = waitTime; }
        if (waitTime < minTime) { minTime = waitTime; }
        sumTime += waitTime;

        double waitTimeinQ = e->cust->totalWaitTimeinQ;
        if (waitTimeinQ > maxTimeinQ) {maxTimeinQ = waitTimeinQ; }
        if (waitTimeinQ < minTimeinQ) {minTimeinQ = waitTimeinQ; }
        sumTimeinQ += waitTimeinQ;

        free(e->cust);                                                  //Free this space allocated for this customer
    }
    else if (arrComps[e->compID].ComponentType == QUEUE_STATION) {
        fprintf(fwrite, "Processing arrival event in Queue Station, component ID is %d, current time is %f\n", e->compID, CurrentTime());
        ptrFIFO = (FIFOQueue*)arrComps[e->compID].ptrComp;
        if (ptrFIFO->first == NULL && ptrFIFO->last == NULL) { // FIFO queue is empty
            ptrFIFO->first = e->cust;
            ptrFIFO->last = e->cust;

            e->cust->enterStaTime = CurrentTime();             //The time customer enters this queue

            //fprintf(fwrite, "CurrentTime is %f\n", e->cust->enterStaTime);

            if ((geneEvent = malloc(sizeof(struct EventData))) == NULL) { fprintf(stderr, "malloc error\n"); exit(1); } //new departure event
            geneEvent->eventType = DEPARTURE;                   //First in, first serve and departure
            geneEvent->cust = e->cust;
            geneEvent->compID = e->compID;

            double serveTime = expTime(ptrFIFO->servetime);       // Get the exp serve time
            //double serveTime = ptrFIFO->servetime;       

            ts = CurrentTime() + serveTime;
            //fprintf(fwrite, "Serve time is %f\n********average serve time is %f\n", serveTime, ptrFIFO->servetime);
            Schedule(ts, geneEvent);
			if (ts > timeLimit) free(geneEvent);
        }
        else { // fifo queue is not empty
            e->cust->enterStaTime = CurrentTime();   
            ptrFIFO->last->next = e->cust;
            ptrFIFO->last = ptrFIFO->last->next;
        }
    }
    else { fprintf(stderr, "Unexpected component type for arrival\n"); exit(1); }
    fclose(fwrite);
}


void Departure(struct EventData* e, char* outfile) {
    FILE* fwrite;
    if ((fwrite = fopen(outfile, "a")) == NULL) {
        fprintf(stderr, "cannot open output file\n"); exit(1);
    }
    EventData* geneEvent;
    EventData* nextDeparture;
    double ts;
    FIFOQueue* ptrFIFO;
    if (e->eventType != DEPARTURE) { fprintf(stderr, "Unexpected event type\n"); exit(1); }
    fprintf(fwrite, "Processing departure event in Queue Station, component ID is %d, current time is %f\n", e->compID, CurrentTime());

    ptrFIFO = (FIFOQueue*)arrComps[e->compID].ptrComp;

    if (ptrFIFO->first == ptrFIFO->last) { // only one customer in the queue
        if ((geneEvent = malloc(sizeof(struct EventData))) == NULL) { fprintf(stderr, "malloc error\n"); exit(1); }
        double rand = urand();
        double prob = 0;
        int i = 0;
        for (; ; i++) {
            prob += ptrFIFO->destPossi[i];
            if (rand < prob) break;
        }

        geneEvent->eventType = ARRIVAL;
        geneEvent->cust = ptrFIFO->first;
        geneEvent->compID = ptrFIFO->destComp[i];
        fprintf(fwrite, "Next component ID is %d\n", geneEvent->compID);

        geneEvent->cust->totalWaitTimeinQ = CurrentTime() - geneEvent->cust->enterStaTime;

        ptrFIFO->numberofServedCust++;
        ptrFIFO->totalWaitTime += (CurrentTime() - geneEvent->cust->enterStaTime);

        ts = CurrentTime();
        Schedule(ts, geneEvent);
		if (ts > timeLimit) free(geneEvent);
        ptrFIFO->first = NULL;
        ptrFIFO->last = NULL;
    }
    else { // more than one customer in fifo queue
        if ((geneEvent = malloc(sizeof(struct EventData))) == NULL) { fprintf(stderr, "malloc error\n"); exit(1); }
        double rand = urand();
        double prob = 0;
        int i = 0;
        for (; ; i++) {
            prob += ptrFIFO->destPossi[i];
            if (rand < prob) break;
        }
        geneEvent->eventType = ARRIVAL;
        geneEvent->cust = ptrFIFO->first;
        geneEvent->compID = ptrFIFO->destComp[i];
        fprintf(fwrite, "Next component ID is %d\n", geneEvent->compID);




        double serveTime = randexp(ptrFIFO->servetime);
        //double serveTime = ptrFIFO->servetime;              


        geneEvent->cust->totalWaitTimeinQ = CurrentTime() - geneEvent->cust->enterStaTime;

        ptrFIFO->numberofServedCust++;
        ptrFIFO->totalWaitTime += CurrentTime() - geneEvent->cust->enterStaTime;

        ts = CurrentTime();
        Schedule(ts, geneEvent);
		if (ts > timeLimit) free(geneEvent);
        if ((nextDeparture = malloc(sizeof(struct EventData))) == NULL) { fprintf(stderr, "malloc error\n"); exit(1); }
        ptrFIFO->first = ptrFIFO->first->next;
        nextDeparture->eventType = DEPARTURE;
        nextDeparture->cust = ptrFIFO->first;
        nextDeparture->compID = e->compID;
        ts = CurrentTime() + serveTime;
        Schedule(ts, nextDeparture);
		if (ts > timeLimit) free(nextDeparture);
    }
    fclose(fwrite);
}

//--------------------------------------------------------------------------------------------
void configFinished(int numofComps, char* outfile) {
    FILE* fwrite;
    if ((fwrite = fopen(outfile, "a")) == NULL) {
        fprintf(stderr, "cannot open output file\n"); exit(1);
    }
    //printf("Loading configuration successful\n");
    for (int i = 0; i < numofComps; ++i) {
        if (arrComps[i].ComponentType == GENERATOR) {
            Generator* ptrGene = arrComps[i].ptrComp;
            fprintf(fwrite, "\nGenerator's ID is %d, its mean interval time is %lf, id of its dest is %d\n", i, ptrGene->geneTime, ptrGene->destComp);
        }
        else if (arrComps[i].ComponentType == QUEUE_STATION) {
            FIFOQueue* ptrSta = arrComps[i].ptrComp;
            fprintf(fwrite, "Station's ID is %d, its mean interval time is %lf\n", i, ptrSta->servetime);
            for (int j = 0; j < ptrSta->routingNum; ++j) {
                fprintf(fwrite, "ID of its dest are %d, corresponding possibility is %lf\n", ptrSta->destComp[j], ptrSta->destPossi[j]);
            }
        }
        else {
            fprintf(fwrite, "Exit's ID is %d\n", i);
        }

        fprintf(fwrite, "\n");
    }
    fclose(fwrite);
}

int inputConfig(char* config, char* outfile) {
    FILE* configFile;
    if ((configFile = fopen(config, "r")) == NULL)
        fprintf(stderr, "Cannot find config file\n");

    //Get the amount of components
    int numofComps = 0;
    fscanf(configFile, "%d", &numofComps);   //get number of all Comps

    int compID;
    char typeName;    //put G into comps[0]
    double avgIntTime;

    int geneDestID;
    int routingNum;
    int* staDestID;
    double* possibility;
    //    double servTime;
    //    char statTypeNanme;

    for (int i = 0; i < numofComps; i++) {

        fscanf(configFile, "%d %c %lf", &compID, &typeName, &avgIntTime);
        if (typeName == 'G') {
            fscanf(configFile, "%d", &geneDestID);
            makeGenerator(compID, avgIntTime, geneDestID, outfile);
        }

        else if (typeName == 'E')
            makeExit(compID, outfile);

        else {
            fscanf(configFile, "%d", &routingNum);

            possibility = (double*)malloc(routingNum * sizeof(double));
            staDestID = (int*)malloc(routingNum * sizeof(int));

            for (int j = 0; j < routingNum; j++) {
                fscanf(configFile, "%lf", &possibility[j]);
            }

            for (int j = 0; j < routingNum; j++) {
                fscanf(configFile, "%d", &staDestID[j]);
            }
            makeSta(compID, avgIntTime, staDestID, possibility, routingNum, outfile);
        }

    }
    fclose(configFile);
    //printf("Loading configuration successful\n", staID, servTime);
    return numofComps;
}

//state variable functions
//static 1.
int getInCust() {
    return inCust;
}
int getOutCust() {
    return outCust;
}
//static 2.
double getMinTime() {
    return minTime;
}
double getMaxTime() {
    return maxTime;
}
double getAvgTime() {
    return (sumTime / outCust);
}
//static 3.
double getMinTimeinQ() {
    return minTimeinQ;
}
double getMaxTimeinQ() {
    return maxTimeinQ;
}
double getAvgTimeinQ() {
    return (sumTimeinQ / outCust);
}
// static 4.
double avgTimeofSta(int CompID) {
    struct FIFOQueue* pFIFO = (struct FIFOQueue*) arrComps[CompID].ptrComp;
    if (pFIFO->numberofServedCust == 0)                                    //if there is no customer served in station
        return 0;
    return (pFIFO->totalWaitTime / pFIFO->numberofServedCust);
}

void freeArr(int numofComps){
    for (int i = 1; i < numofComps - 1; i++) {
        FIFOQueue* p;
        p = (FIFOQueue*) arrComps[i].ptrComp;
        free(p->destComp);
        free(p->destPossi);
        if (p->first == NULL && p->last == NULL) {;}
        else {
            while (p->first != NULL ){
                Customer * c = p->first;
                p->first = c->next;
                free(c);
		if(p->first == p->last)
			break;
		
            }
	    
free(p->first);
			//free(p->last);
        }
    }
    for (int i = 0; i < numofComps; i++){
        free(arrComps[i].ptrComp);
    }
    
}
//typedef struct FIFOQueue {
//	struct Customer* first;     // pointer to first customer in queue
//	struct Customer* last;      // pointer to last customer in queue
//	double servetime;           // average serve time of a QS
//	int* destComp;              // ID of next component customers are sent to
//	double totalWaitTime;        // total waiting time
//	int numberofServedCust;
//	double* destPossi;
//	int routingNum;
//}FIFOQueue;

int main(int argc, char* argv[]) {
    timeLimit =  atof(argv[1]);
    char* config = argv[2];
    char* outfile = argv[3];

    //timeLimit = 240;
    //char* config = "config.txt";
    //char* outfile = "output.txt";

    int numofComps = inputConfig(config, outfile);
    configFinished(numofComps, outfile);

    RunSim(timeLimit, outfile);

    FILE* fwrite;
    if ((fwrite = fopen(outfile, "a")) == NULL) {
        fprintf(stderr, "cannot open output file\n"); exit(1);
    }
    fprintf(fwrite, "\n*************** Data ***************\n");
    fprintf(fwrite, "number of customers entering the sysytem = %d\n", getInCust());
    fprintf(fwrite, "number of customers exiting the system = %d\n", getOutCust());
    fprintf(fwrite, "min time in system = %f, max time in system = %f, average time in system = %f\n", getMinTime(), getMaxTime(), getAvgTime());
    fprintf(fwrite, "min time in queue = %f, max time in queue = %f, average time in queue = %f\n", getMinTimeinQ(), getMaxTimeinQ(), getAvgTimeinQ());
    for (int i=1; i<numofComps-1; i++) {
        fprintf(fwrite, "avgTimeofSta(%d) = %f\n", i, avgTimeofSta(i));
    }
    fclose(fwrite);

	//freeArr(numofComps);

    return 0;
}



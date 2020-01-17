Author          : Hongyu Su and Yuhang Chen
Created         : November 5 , 2019
Last Modified   : November 5, 2019

Affiliation          : Georgia Institute of Technology


Description
-------------
engine.c
This part of the program is independent of the simulation model. It includes a data structure called the future event list (FEL) that is a priority queue that contains the set of events that have been scheduled, but have not yet been processed. The simulation engine holds the main event processing loop which repeatedly (1) removes the smallest timestamped event from the FEL, and (2) calls a function (defined in the simulation application, discussed next) to simulate that event. The loop continues processing events until some termination condition is met, e.g., simulation time reaches a certain value. The simulation maintains a clock variable indicating how far the simulation has advanced in simulation time. It also includes a function called by the simulation application to schedule a new event. The main event processing loop updates the clock variable in each iteration of the loop to the timestamp of the event it just removed from the FEL.

application.c
This part of the program contains code to model the physical system. It includes a set of state variables that represent the current state of the system being modeled. It also includes one or more functions or event handler procedures, one for each type of event modeled by the simulation. The event handler procedures collectively model the operation of the system being simulated. To develop the simulation application, one must define the state variables and the different types of events that are modeled, and implement a function for each different event type.

sim.h
Application independent simulation engine interface.



Installation
------------

For Part one to install, simply run

    gcc -std=c99 -lm application.c engine.c -o cpssim

Note: Note: The -std=c99 flag is required on some C compilers 
to allow variable declarations at arbitrary locations in a function.
And -lm is required because the program includes math.h


Execution
----------

Assuming your executable is called "cpssim", run it typing  

    ./cpssim 240 config.txt output.txt
    
The program should take the following command line parameters in sequence

The program has three command line parameters.
The first one is the time of simulation.
The second one is the name of the configuration file.
The third one is the name of the output file.

 



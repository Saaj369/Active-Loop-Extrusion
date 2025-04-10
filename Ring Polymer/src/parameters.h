// Global parameter definitions for the simulation.
// Includes tunable simulation constants and system setup variables.

# include <cmath>
# include <iostream>

#ifndef PARAMETERS_H
#define PARAMETERS_H

#define PI 3.141592653589793

const int NUM_PARTICLES = 500;
const double TEMPERATURE = 1.0;

const double SPRING_CONSTANT = 1000.0;
const double SPRING_L = std::pow(2.0,1.0/6);
// const float RING_DIAMETER = NUM_PARTICLES * SPRING_L / PI;
const double RING_DIAMETER = SPRING_L/std::sin(PI/NUM_PARTICLES);

const double CUTTOFF_RADIUS = std::pow(2.0,1.0/6);
const double dt = 0.00005;
/* 
steps = 200,000,000
dissociationtime = 200,000
DISSOCIATION_TIME_INTERVAL = 10;
*/

// If utilising Neighbour List
const double NL_cutoff_sq = (CUTTOFF_RADIUS+0.5)*(CUTTOFF_RADIUS+0.5);
const double SKIN = sqrt(NL_cutoff_sq) - CUTTOFF_RADIUS;

// If simulating SMC
const int NUM_SMC = 75;
const double K_SMC_POLY = 100;
const double K_SMCs = 1000;
const double REE_FORCE = 500;
const double DISSOCIATION_TIME_INTERVAL = 10;
const int DISSOCIATION_STEP_INTERVAL = DISSOCIATION_TIME_INTERVAL/dt;


//Debug
const double ideal_sim_max_distance = RING_DIAMETER*PI*0.5;

#endif
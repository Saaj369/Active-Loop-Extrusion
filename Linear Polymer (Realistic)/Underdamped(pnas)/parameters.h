# include <cmath>
# include <iostream>
# include <string>
#include <cstddef>

#ifndef DNA_PARAMETERS_H
#define DNA_PARAMETERS_H

#define PI 3.141592653589793

const int NUM_PARTICLES = 1000;
const int NUM_CTCF = 146; // 146
const double h = 5; //5, 2
const double Mu = 4.4; //5.3, 2.3

const double SPRING_CONSTANT = 1000.0; //1000
const double SPRING_L = 1; //std::pow(2.0,1.0/6)

const double CUTTOFF_RADIUS = std::pow(2.0,1.0/6);
const double CUTTOFF_RADIUS2 = std::pow(2.0,1.0/3.0);
const double dt = 0.01; //0.00005
const double dth = 0.5 * dt;
const double C = (2 * dt / (2 + dt) );
const double A = ((2 - dt) / (2 + dt));
const double B = sqrt(0.5 * dt);
const double dtau = 10;
const int MC_interval =   dtau/dt; //10

// If utilising Neighbour List
const double NL_cutoff_sq = (CUTTOFF_RADIUS+0.7)*(CUTTOFF_RADIUS+0.7);
const double SKIN = sqrt(NL_cutoff_sq) - CUTTOFF_RADIUS;
const double CSKIN = SKIN / (2*sqrt(3));

// If simulating SMC
const double K_SMC_POLY = 100; //100

//Probabilities

const double Un_pob = exp(-h); // exp(-h) 0
const double B_prob = exp(-Mu-h); //exp(-Mu-h) 0.4
const double Free_jump_prob = 0.5; // 0.4
const double POc_jump_prob = 0.5; // 0.4

//Debug
// const double ideal_sim_max_distance = RING_DIAMETER*PI*0.5;
const int BD_RAND_SEED =  88284;
const int binding_RAND_SEED = 96012;
const int unbind_RAND_SEED =  4318;
const int choice_RAND_SEED =  94874;
const int rGenFreq = 1;

// const int BD_RAND_SEED =  82901;
// const int binding_RAND_SEED = 20518;
// const int unbind_RAND_SEED =  937548;
// const int choice_RAND_SEED =  642329;
// const int rGenFreq = 1;

// color coding array

const std::string color_code[5] = {
    "255 0 0",   // Red
    "0 0 255",   // Blue
    "0 255 0",   // Green
    "255 0 255", // Magenta
    "255 255 0"  // Yellow
};


// Binary writing parameters
size_t BUFFER_SIZE = NUM_PARTICLES * (sizeof(int) + 3 * sizeof(float) + 3 * sizeof(uint8_t));

#endif
# include <cmath>
# include <iostream>
# include <string>
#include <cstddef>

#ifndef DNA_PARAMETERS_H
#define DNA_PARAMETERS_H

#define PI 3.141592653589793

const int NUM_PARTICLES = 1000;
const double h = 5;
const double Mu = 5.3; // 5.3,  5.0, 4.7, 4.4, 4.0
const int N_tad = 1;
const int TAD[N_tad] = {200};

//Probabilities

const double Un_pob = exp(-h);
const double B_prob = exp(-Mu-h);
const double Free_jump_prob = 0.5;
const double POc_jump_prob = 0.5;

const int BD_RAND_SEED = 23918;
const int binding_RAND_SEED = 903847;
const int unbind_RAND_SEED = 152338;
const int choice_RAND_SEED = 928014;


#endif
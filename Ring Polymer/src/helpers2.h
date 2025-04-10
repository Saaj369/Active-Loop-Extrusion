// Additional helper declarations for extended MD functionality.
// Complements helpers.h with more specialized utilities.

# include <array>
# include <cstring>
# include <iostream>
# include "parameters.h"
# include "helpers.h"
using namespace std;
typedef std::array<std::array<double,3>, NUM_PARTICLES> c2D_t;
typedef std::array<int, NUM_PARTICLES*NUM_PARTICLES> NList_t;
typedef std::array< std::array<int,2>, NUM_PARTICLES> NLindex2d_t;

#ifndef HELPERS2_H
#define HELPERS2_H

void initialise(c2D_t &POS, c2D_t &VEL);
void initialiseBD(c2D_t &POS);
void eval_force(c2D_t &NF, c2D_t &POS, double &PE);
void eval_force_using_NL(c2D_t &NF, c2D_t &POS, double &PE, NList_t &NL, NLindex2d_t &NLsize);
void velocity_verlet_step(c2D_t &POS, c2D_t &VEL, c2D_t &F, c2D_t &NF, double &PE);
void velocity_verlet_step_using_NL(c2D_t &POS,c2D_t &POS0, c2D_t &VEL, c2D_t &F, c2D_t &NF, double &PE,NList_t &NL, NLindex2d_t &NLsize);
void velocity_verlet_BD(c2D_t &POS, c2D_t &F);
void velocity_verlet_BD_using_NL(c2D_t &POS,c2D_t &F);
void store_position_data(c2D_t &POS, string file_name, int index);
void store_KE_PE_momentum(double &PE, c2D_t &VEL,string file_name, int index);

#endif
 // Header file for general helper functions used in polymer MD simulation.
// Includes utility and math functions for position, force, etc.

 # include <iostream>
 # include <array>
 # include "parameters.h"

#ifndef HELPERS_H
#define HELPERS_H

using namespace std;
typedef std::array<std::array<double,3>, NUM_PARTICLES> c2D_t;
typedef std::array<double, 3> vec_t;
typedef std::array<int, NUM_PARTICLES*NUM_PARTICLES> NList_t;
typedef std::array< std::array<int,2>, NUM_PARTICLES> NLindex2d_t;
typedef std::array<std::array<double, 3>,2*NUM_SMC> s2D_t;

// using vector = vec_t;

vec_t difference(vec_t &v1, vec_t &v2);
vec_t multiply(vec_t &v1, double c);
double magnitude(vec_t &v);
double dot(vec_t &v1, vec_t &v2);

vec_t spring_force(vec_t &r12, double r);
vec_t lj_force(vec_t &r12, double r);

void neighbour_list(NList_t &NL, NLindex2d_t &NLsize,c2D_t &POS);
bool check_if_NList_requires_update(c2D_t &POS0, c2D_t &POS);
void copy_POS_data(c2D_t &POS,c2D_t &POS0);
void copy_SPOS_data(s2D_t &SPOS,s2D_t &SPOS0);


#endif
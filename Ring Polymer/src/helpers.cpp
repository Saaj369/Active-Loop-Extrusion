 // Implementation of helper functions for molecular dynamics simulation.
// Supports core MD logic like force calculation, integration, etc.

 # include <iostream>
 # include <array>
 # include <cmath>
 # include "parameters.h"
 # include "helpers.h"

using namespace std;
typedef std::array<std::array<double, 3>, NUM_PARTICLES> c2D_t;
typedef std::array<double, 3> vec_t;
typedef std::array<int, NUM_PARTICLES*NUM_PARTICLES> NList_t;
typedef std::array< std::array<int,2>, NUM_PARTICLES> NLindex2d_t;
typedef std::array<std::array<double, 3>,2*NUM_SMC> s2D_t;
// using vector = vec_t;

vec_t difference(vec_t &v1, vec_t &v2){
    vec_t v = {};
    for (int i = 0; i < 3; i++){
        v[i] = v1[i] - v2[i];
    }
    return v;
}

vec_t multiply(vec_t &v1, double c){
    vec_t v = {};
    for(int i = 0; i < 3; i++){
        v[i] = v1[i]*c;
    }
    return v;
}

double magnitude(vec_t &v){
    double mag = 0;
    for(int i = 0; i < 3; i++){
        mag += v[i]*v[i];
    }
    mag = sqrt(mag);
    return mag;
}

double dot(vec_t &v1, vec_t &v2){
    double dot = 0;
    for( int i = 0; i < 3; i++){
        dot += v1[i]*v2[i];
    }
    return dot;
}

vec_t spring_force(vec_t &r12, double r){
    double factor = -(SPRING_CONSTANT/TEMPERATURE) * (r - SPRING_L)/r;
    vec_t springForce =  multiply(r12,factor);
    
    return springForce;
}

vec_t lj_force(vec_t &r12, double r){
    double factor = (24/TEMPERATURE)*(2*pow(r, -14) - pow(r,-8));
    vec_t ljForce =  multiply(r12,factor);

    return ljForce;
}

void neighbour_list(NList_t &NL, NLindex2d_t &NLsize,c2D_t &POS){
    int NL_index = 0;
    for(int i = 0; i < NUM_PARTICLES-1; i++){
        NLsize[i][0] = NL_index; // starting for ith particle 
        for(int j = i+1; j < NUM_PARTICLES; j++){
            vec_t rij = difference(POS[i],POS[j]); 
            double dist_sq = dot(rij,rij);
            if( dist_sq < NL_cutoff_sq){
                NL[NL_index] = j;
                NL_index += 1;
            }
        }
        NLsize[i][1] = NL_index; // ending for ith particle
    }
}

bool check_if_NList_requires_update(c2D_t &POS0, c2D_t &POS){
    double maxDisp = 0.0;
    for(int i = 0; i < NUM_PARTICLES; i++){
        for(int j = 0; j < 3; j++){
            maxDisp = max(maxDisp, abs(POS[i][j] - POS0[i][j]));
        }
    }
    maxDisp = 2.0 * sqrt(3.0*maxDisp*maxDisp);
    bool update = (maxDisp > SKIN);
    return update;
}

void copy_POS_data(c2D_t &POS,c2D_t &POS0){
    for(int i = 0; i < NUM_PARTICLES; i++){
        for(int j = 0; j < 3; j++){
            POS0[i][j] = POS[i][j];
        }
    }
}

void copy_SPOS_data(s2D_t &SPOS,s2D_t &SPOS0){
    for(int i = 0; i < 2*NUM_SMC; i++){
        for(int j = 0; j < 3; j++){
            SPOS0[i][j] = SPOS[i][j];
        }
    }
}
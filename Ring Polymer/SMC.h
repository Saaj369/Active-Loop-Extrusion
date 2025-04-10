// Header file defining structures and interfaces for SMC logic.
// Part of the SMC (Structural Maintenance of Chromosomes) model.

# include <array>
# include <cstring>
# include <iostream>
#include "../src/parameters.h"
#include "../src/helpers.h"

#ifndef SMC_H
#define SMC_H

typedef std::array<double, 3> vec_t;
typedef std::array<std::array<double, 3>, NUM_PARTICLES> c2D_t;
typedef std::array<int, NUM_PARTICLES * NUM_PARTICLES> NList_t;
typedef std::array<std::array<int, 2>, NUM_PARTICLES> NLindex2d_t;
typedef std::array<std::array<double, 3>,2*NUM_SMC> s2D_t;
typedef std::array<int,2*NUM_SMC> links_t;

void attach_SMCs(c2D_t &POS, s2D_t &SPOS, links_t &Slinks,std::array<int,NUM_PARTICLES> &sites_availability);
void eval_force_SMCs(s2D_t &SPOS, s2D_t &SF, links_t &Slinks);
void inject_focre_SMC_poly(c2D_t &POS, c2D_t &F,s2D_t &SPOS, s2D_t &SF, links_t &Slinks );
void velocity_verlet_SMC(s2D_t &SPOS,s2D_t &SF);
void random_detach_and_attach(s2D_t &SPOS, c2D_t &POS, links_t &Slinks, std::array<int,NUM_PARTICLES> &sites_availability);
void check_proximity_and_take_action(c2D_t &POS, s2D_t &SPOS, links_t &Slinks,links_t &if_Slinks_blocked);
void save_SMC_polymer_xyz(std::ofstream& out, c2D_t &POS, s2D_t &SPOS, int step);
void save_SMC_links(std::ofstream& out, links_t &Slinks, int step);
void write_config(string configFile , int steps, int equisteps, int interval);
bool validate_sim(c2D_t &POS, s2D_t &SPOS, int &index, string &msg);
void save_debug_data(links_t &Slinks, links_t &if_Slinks_blocked, s2D_t &SF, c2D_t &F, std::array<int,NUM_PARTICLES> &sites_availability, string fileSMC, string fileC, int step);


#endif
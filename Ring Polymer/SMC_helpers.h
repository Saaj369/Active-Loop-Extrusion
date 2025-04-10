// Helper function declarations for SMC dynamics.
// Used by SMC.cpp and the main simulation file.

# include <array>
# include <cstring>
# include <iostream>
#include "../src/parameters.h"
#include "../src/helpers.h"

#ifndef SMC_HELPERS_H
#define SMC_HELPERS_H

typedef std::array<double, 3> vec_t;
typedef std::array<std::array<double, 3>, NUM_PARTICLES> c2D_t;
typedef std::array<int, NUM_PARTICLES * NUM_PARTICLES> NList_t;
typedef std::array<std::array<int, 2>, NUM_PARTICLES> NLindex2d_t;
typedef std::array<std::array<double, 3>,2*NUM_SMC> s2D_t;
typedef std::array<int,2*NUM_SMC> links_t;
typedef std::array<int,NUM_PARTICLES> states_t;

void calculate_sites_availability(links_t &Slinks, states_t &sites_availability);
int ranndom_site_chooser_for_smcs(std::array<int,NUM_PARTICLES> &sites_availability);
vec_t springforce_smc_polymer(vec_t &r12, double r);
vec_t springforce_smcs(vec_t &r12, double r);
void update_if_Slinks_blocked(links_t &if_Slinks_blocked, links_t &Slinks);
bool is_leftlink_blocked(int link, links_t &Slinks, links_t &if_Slinks_blocked);
bool is_rightlink_blocked(int link, links_t &Slinks, links_t &if_Slinks_blocked);
bool check_and_update_blockness(int index, links_t &Slinks, links_t &if_Slinks_blocked);
void display_progress(int totalSteps, int finishedSteps, string msg);


#endif
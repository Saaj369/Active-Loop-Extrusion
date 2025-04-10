// Helper function implementations for SMC-related computations.
// Contains utility functions used across the SMC module.

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <cstring>
#include <unistd.h>
#include <array>
#include <algorithm>
#include <iomanip>
#include <cstdlib> // For std::exit
#include "../src/parameters.h"
#include "../src/helpers.h"
#include "../src/helpers2.h"

using namespace std;
typedef std::array<double, 3> vec_t;
typedef std::array<std::array<double, 3>, NUM_PARTICLES> c2D_t;
typedef std::array<int, NUM_PARTICLES * NUM_PARTICLES> NList_t;
typedef std::array<std::array<int, 2>, NUM_PARTICLES> NLindex2d_t;
typedef std::array<std::array<double, 3>, 2 * NUM_SMC> s2D_t;
typedef std::array<int, 2 * NUM_SMC> links_t;
typedef std::array<int, NUM_PARTICLES> states_t;


/* This function updates sites_Availability array based on the Slinks. */
void calculate_sites_availability(links_t &Slinks, states_t &sites_availability)
{
    fill(&sites_availability[0], &sites_availability[0] + NUM_PARTICLES, 0);
    for (int i = 0; i < 2 * NUM_SMC; i++)
    {
        sites_availability[Slinks[i]] = 1; // Blocked monomer/bead
    }
}
//-----------------------------------------------------------------------------------

/*
This function takes array for sites availability and generates a random site for SMCs to be hooked
*/
int ranndom_site_chooser_for_smcs(std::array<int,NUM_PARTICLES> &sites_availability){
    //Calculating total available sites
    auto available_sites = [](std::array<int,NUM_PARTICLES> &sites_availability) -> int{
        int Nsites = 0;
        for(int i = 0; i < NUM_PARTICLES; i++){
            if(sites_availability[i] + sites_availability[(i+1)%NUM_PARTICLES] == 0){Nsites++;}}
        return Nsites;
    };
    // finding out index of site
    auto index_of_site = [](std::array<int,NUM_PARTICLES> &sites_availability,int rndsite) -> int{
        for(int i = 0 ; i < NUM_PARTICLES; i++){
            if(sites_availability[i] + sites_availability[(i+1)%NUM_PARTICLES] == 0){
                rndsite--;
                if(rndsite == 0){
                    return i;
                }}}
        cout<<"Something went wrong in attach_SMCs function"<<endl;
        return -1;
        };
    // Initialising random number generator    
    random_device rd{};
    mt19937 gen{rd()};
    int NO_OF_SITES = available_sites(sites_availability);
    if(NO_OF_SITES==0){
        // cout<<"No sites available for SMCs to attach"<<endl;
        return -1;
    }
    else{
        uniform_int_distribution<> dis(1,NO_OF_SITES); // generating a random site number
        int random_site = dis(gen);
        int index = index_of_site(sites_availability,random_site); // got the index of the site
        return index;
    }
}
//--------------------------------------------------------------------------------------------------

/* This function returns spring force between SMC monomer and Ring Polymer monomer */
vec_t springforce_smc_polymer(vec_t &r12, double r)
{
    double factor = -(K_SMC_POLY / TEMPERATURE);
    vec_t springForce = multiply(r12, factor);

    return springForce;
}
//-------------------------------------------------------------------------------------------------


/* This function returns spring force between SMCs only */
vec_t springforce_smcs(vec_t &r12, double r)
{
    double factor = -(K_SMCs / TEMPERATURE);
    vec_t springForce = multiply(r12, factor);

    return springForce;
}
//-------------------------------------------------------------------------------------------------


/* This function updates the if_Slinks_blocked to latest state of the system */
void update_if_Slinks_blocked(links_t &if_Slinks_blocked, links_t &Slinks)
{
    for(int i = 0; i < 2*NUM_SMC; i++){
        for(int j = 0; j < 2*NUM_SMC; j++){
            if(i%2==0){
                if(Slinks[i] == (Slinks[j]+1)%NUM_PARTICLES){
                    if_Slinks_blocked[i] = 1;
                    break;
                }
                else{
                    if_Slinks_blocked[i] = 0;
                }
            }
            else{
                if(Slinks[i] == (Slinks[j]-1+NUM_PARTICLES)%NUM_PARTICLES){
                    if_Slinks_blocked[i] = 1;
                    break;
                }
                else{
                    if_Slinks_blocked[i] = 0;
                }
            }
        }
    }
}
//-------------------------------------------------------------------------------------------


/* This function is reponsible for checking blockness and updating*/
bool check_and_update_blockness(int index, links_t &Slinks, links_t &if_Slinks_blocked)
{
    if(index%2==0){
        for(int k = 0; k < 2*NUM_SMC; k++){
            if(Slinks[index] == (Slinks[k]+1)%NUM_PARTICLES){
                if_Slinks_blocked[index] = 1;
                if(k%2!=0){
                    if_Slinks_blocked[k] = 1;
                }
                return true;
            }
        }
    }
    else{
        for(int k = 0; k < 2*NUM_SMC; k++){
            if(Slinks[index] == (Slinks[k]-1+NUM_PARTICLES)%NUM_PARTICLES){
                if_Slinks_blocked[index] = 1;
                if(k%2==0){
                    if_Slinks_blocked[k] = 1;
                }
                return true;
            }
        }
    }
    return false;
}
//-------------------------------------------------------------------------------------------


/* A function to display progress in terminal */
void display_progress(int totalSteps, int finishedSteps, string msg)
{
    cout << "\033[2J\033[H";
    cout.flush();
    // Display progress
    cout << finishedSteps << "/" << totalSteps << "------->" << msg << 100.0f * finishedSteps / totalSteps << " %" << "\r";
    cout.flush();
}
//------------------------------------------------------------------------------------------
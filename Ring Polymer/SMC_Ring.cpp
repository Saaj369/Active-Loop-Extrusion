// Main file for running SMC Ring Polymer simulations.
// Integrates SMC logic and molecular dynamics of the polymer.
// Entry point of the simulation.

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <cstring>
#include <unistd.h>
#include <chrono>
#include <iomanip>
#include <array>
#include <cstdlib> // For std::exit
#include "../src/parameters.h"
#include "../src/helpers.h"
#include "../src/helpers2.h"
#include "SMC.h"
#include "SMC_helpers.h"

using namespace std;
using namespace std::chrono;
typedef std::array<std::array<double, 3>, NUM_PARTICLES> c2D_t;
typedef std::array<int, NUM_PARTICLES * NUM_PARTICLES> NList_t;
typedef std::array<std::array<int, 2>, NUM_PARTICLES> NLindex2d_t;
typedef std::array<std::array<double, 3>,2*NUM_SMC> s2D_t;
typedef std::array<int,2*NUM_SMC> links_t;


void smc_simulateBDwithNL(int steps, int equiSteps, int data_entry_interval)
{
    // Polymer data containers
    c2D_t POS = {}; // To store position data in present step
    c2D_t POS0 = {}; // To store position data when neighbourlist gets initialised, this is used to decide weather to update NList
    c2D_t F = {}; // To store force data in present step

    // c2D_t prevPOS = {}; // (for debugging)  To store position data of previous step
    // s2D_t prevSPOS = {}; // (for debugginh) To store position data of SMCs for previous step

    s2D_t SPOS = {}; // To store data of SMCs in present step
    s2D_t SF = {}; // To store forces on SMCs in present step
    links_t Slinks = {}; // To store the location where each SMC monomer is attached
    links_t if_Slinks_blocked = {}; // To store information if a particular monomer of a SMC is blocked or not.
    std::array<int,NUM_PARTICLES> sites_availability = {}; // To store sites availability info on ring polymer

    NList_t NL = {};         // Neighbour List
    NLindex2d_t NLsize = {}; // Neighbour List size
    double PE = 0; // NO use of it, but few function require this as an argument

    // filenames with directory to save data
    string configFile = "dataSMC/config.txt";
    string POSfile = "dataSMC/SMCRingPOS.xyz";
    string SMClinks = "dataSMC/SMClinks.txt";
    string ePOSfile = "dataSMC/ePOS.xyz";
    // string VELfile = "dataSMC/SMCRingVEL.txt";
    // string datafile = "dataSMC/SMCDCdata_NL.txt";
    // string debugfile = "dataSMC/debug.xyz";
    // string debugSMC = "dataSMC/debugSMC.txt";
    // string debugC = "dataSMC/debugC.txt";


    // information tags for various simulation state
    string msg_init = "(Equilibrating)";
    string msg_simulate = "(Simulation Running)";

    // Writing a Configuration file beforehand
    write_config(configFile,steps,equiSteps,data_entry_interval);


    // Initialisation
    initialiseBD(POS); // Initialising a Ring Polymer
    neighbour_list(NL, NLsize, POS); // Initialising a neighbour list for the ring polymer
    eval_force_using_NL(F, POS, PE, NL, NLsize); // Evaluating force on ring polymer (excluding SMCs)
    copy_POS_data(POS, POS0); // Copying current position to POS0 because we will be updating NList.
   

    // ******************Velocity Verlet for equilibirum*********************8
    ofstream eout(ePOSfile, ios::app);
    for (int i = 0; i < equiSteps / data_entry_interval; i++)
    {
        save_SMC_polymer_xyz(eout, POS,SPOS,i);
        display_progress(equiSteps / data_entry_interval, i + 1, msg_init);
        for (int j = 0; j < data_entry_interval; j++)
        {
            // ####### Validating if the simulation haven't broken #######
            int index; string msg;
            // if(validate_sim(POS,SPOS,index,msg)){
            //     std::cerr << "Simulation Broke "<<" Particle : "<<msg<<index<<" Step: "<<i*data_entry_interval+j<< std::endl;
            //     std::exit(EXIT_FAILURE);
            // }
            //----------------------------------------------------------------

            // ####### Handling Detachment and Attachment #######
            // if((i*data_entry_interval + j + 1)%DISSOCIATION_STEP_INTERVAL == 0){
            //     random_detach_and_attach(SPOS, POS, Slinks, sites_availability);
            // }
            //-----------------------------------------------------------------

            // ####### Checking if Neighbour list needs to be updated #######
            if (check_if_NList_requires_update(POS0, POS)){
                neighbour_list(NL, NLsize, POS);
                copy_POS_data(POS, POS0);
            }
            //----------------------------------------------------------------

            eval_force_using_NL(F, POS, PE, NL, NLsize); // Evaluating force on RingPolymer (Excluding SMCs)
            eval_force_SMCs(SPOS,SF,Slinks);             // Evaluating forces on SMCs (Excluding Ring Polymer)
            // inject_focre_SMC_poly(POS,F,SPOS,SF,Slinks); // Now Injecting forces on Ring due to SMCs and vice versa
            // velocity_verlet_SMC(SPOS,SF);                // Updating positions on SMCs
            velocity_verlet_BD_using_NL(POS,F);          // Updating force on Ring Polymer
            // check_proximity_and_take_action(POS,SPOS,Slinks, if_Slinks_blocked); // Checking if adjacent beads on ring polymer is closer to SMCs and if yes, updating links
        }
    }
    eout.close();

    attach_SMCs(POS,SPOS,Slinks,sites_availability); //Attaching SMCs to the Ring Polymer
    update_if_Slinks_blocked(if_Slinks_blocked,Slinks); // Updating the info about if Slinks are blocked

    // ***************Velocity verlet with data reading*******************
     ofstream out(POSfile, ios::app); // Creating a file object which will be used to write stuff
    //  ofstream outLinks(SMClinks, ios::app); // Creating file object to write Slink data

    for (int i = equiSteps/data_entry_interval; i < steps/data_entry_interval; i++)
    {
        save_SMC_polymer_xyz(out, POS,SPOS,i);
        // save_SMC_links(outLinks, Slinks, i);

        display_progress((steps - equiSteps) / data_entry_interval, i + 1 - equiSteps/data_entry_interval, msg_simulate);
        for (int j = 0; j < data_entry_interval; j++)
        {
            // ####### Validating if the simulation haven't broken #######
            // int index; string msg;
            // if(validate_sim(POS,SPOS,index,msg)){
                
            //     ofstream slinkData("dataSMC/slinkData.txt", ios::app);
            //     for(int i = 0; i < 2*NUM_SMC; i++){
            //         slinkData << Slinks[i] << " ";
            //     }
            //     slinkData << endl;
            //     slinkData.close();
            //     save_SMC_polymer_xyz(out,POS,SPOS,i+1);
            //     out.close();
            //     std::cerr << "Simulation Broke "<<" Particle : "<<msg<<index<<" Step: "<<i*data_entry_interval+j<< std::endl;
            //     std::exit(EXIT_FAILURE);
            // }
            //----------------------------------------------------------------

            // ####### Handling Detachment and Attachment #######
            if((i*data_entry_interval + j + 1)%DISSOCIATION_STEP_INTERVAL == 0){
                random_detach_and_attach(SPOS, POS, Slinks, sites_availability);
            }
            //-----------------------------------------------------------------

            // ####### Checking if Neighbour list needs to be updated #######
            if (check_if_NList_requires_update(POS0, POS)){
                neighbour_list(NL, NLsize, POS);
                copy_POS_data(POS, POS0);
            }
            //----------------------------------------------------------------

            eval_force_using_NL(F, POS, PE, NL, NLsize); // Evaluating force on RingPolymer (Excluding SMCs)
            eval_force_SMCs(SPOS,SF,Slinks);             // Evaluating forces on SMCs (Excluding Ring Polymer)
            inject_focre_SMC_poly(POS,F,SPOS,SF,Slinks); // Now Injecting forces on Ring due to SMCs and vice versa
            velocity_verlet_SMC(SPOS,SF);                // Updating positions on SMCs
            velocity_verlet_BD_using_NL(POS,F);          // Updating force on Ring Polymer
            check_proximity_and_take_action(POS,SPOS,Slinks, if_Slinks_blocked); // Checking if adjacent beads on ring polymer is closer to SMCs and if yes, updating links
        }
    }
    out.close();
    // outLinks.close();
}


int main()
{
    int steps = 10000000;
    // int steps = 1000000;
    int equiSteps = 1000000;
    int interval = 2000;
    // int interval = 1000;
    auto start = high_resolution_clock::now();
    smc_simulateBDwithNL(steps, equiSteps, interval);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Simulation time : " << duration.count() << " microseconds" << endl;

}
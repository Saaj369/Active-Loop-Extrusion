// Core logic for SMC-related operations.
// Implements functions declared in SMC.h.

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
#include "SMC_helpers.h"

using namespace std;
typedef std::array<double, 3> vec_t;
typedef std::array<std::array<double, 3>, NUM_PARTICLES> c2D_t;
typedef std::array<int, NUM_PARTICLES * NUM_PARTICLES> NList_t;
typedef std::array<std::array<int, 2>, NUM_PARTICLES> NLindex2d_t;
typedef std::array<std::array<double, 3>,2*NUM_SMC> s2D_t;
typedef std::array<int,2*NUM_SMC> links_t;


/*
This function initialises SMCs and attached them to the Ring Polymer.
By updating Slinks and sites_availability arrays
*/
void attach_SMCs(c2D_t &POS, s2D_t &SPOS, links_t &Slinks, std::array<int,NUM_PARTICLES> &sites_availability){
    // Inititlising the site_availability array to zero (All available)
    for(int i = 0; i < NUM_PARTICLES; i++){
        sites_availability[i] = 0; //Available : 0 and Not available : 1
    }

    for(int i = 0; i < 2*NUM_SMC; i = i+2){
        int index = ranndom_site_chooser_for_smcs(sites_availability);
        Slinks[i] = index;
        Slinks[i+1] = (index+1)%NUM_PARTICLES;
        SPOS[i] = POS[index];
        SPOS[i+1] = POS[(index+1)%NUM_PARTICLES];
        sites_availability[index] = 1;
        sites_availability[(index+1)%NUM_PARTICLES] = 1;
    }

}

//---------------------------------------------------------------------------------------------


/* This function evaluates forces between SMC monomers only*/
void eval_force_SMCs(s2D_t &SPOS, s2D_t &SF, links_t &Slinks){
    fill(&SF[0][0], &SF[0][0] + (NUM_SMC*2) * 3, 0);
    for(int i = 0; i < 2*NUM_SMC; i = i+2){
        vec_t r12 = difference(SPOS[i], SPOS[(i+1)]);
        double r = magnitude(r12);
        vec_t Fsp = springforce_smcs(r12,r); // changed and made a dedicated function
        for(int j = 0; j < 3; j++){
            SF[i][j] += Fsp[j];
            SF[i+1][j] -= Fsp[j];
        }
    }
}
//------------------------------------------------------------------------------------------------


/* This function injects forces due to SMCs on Ring Polymer and also viceversa */
void inject_focre_SMC_poly(c2D_t &POS, c2D_t &F,s2D_t &SPOS, s2D_t &SF, links_t &Slinks ){
    // Updating forces between SMC and polymer
    for(int i = 0; i < 2*NUM_SMC; i = i+2){
        vec_t r12_1 = difference(SPOS[i],POS[Slinks[i]]);
        vec_t r12_2 = difference(SPOS[i+1],POS[Slinks[i+1]]);
        double r1 = magnitude(r12_1); double r2 = magnitude(r12_2);
        vec_t Fdm1 = springforce_smc_polymer(r12_1,r1);
        vec_t Fdm2 = springforce_smc_polymer(r12_2,r2);
        for(int j = 0; j < 3; j++){
            SF[i][j] += Fdm1[j];
            SF[i+1][j] += Fdm2[j];
            F[Slinks[i]][j] -= Fdm1[j];
            F[Slinks[i+1]][j] -= Fdm2[j];
        }

    }

    // Injecting Reeling force
    for(int i = 0; i < 2*NUM_SMC; i = i + 2){
        vec_t r12_l = difference(POS[Slinks[i]], POS[(Slinks[i]-1+NUM_PARTICLES)%NUM_PARTICLES]);
        vec_t r12_r = difference(POS[Slinks[i+1]], POS[(Slinks[i+1]+1)%NUM_PARTICLES]);
        double rl = magnitude(r12_l); double rr = magnitude(r12_r);
        double fac_l = REE_FORCE/rl;
        double fac_r = REE_FORCE/rr;
        vec_t leftReelForce = multiply(r12_l,fac_l);
        vec_t rightReelForce = multiply(r12_r,fac_r);

        for(int j = 0; j < 3; j++){
            F[(Slinks[i]-1+NUM_PARTICLES)%NUM_PARTICLES][j] += leftReelForce[j];
            F[(Slinks[i+1]+1)%NUM_PARTICLES][j] += rightReelForce[j];
            SF[i][j] -= leftReelForce[j];
            SF[i+1][j] -= rightReelForce[j];
        }
    }
}
//------------------------------------------------------------------------------------------------

/* This function updates position of SMCs */
void velocity_verlet_SMC(s2D_t &SPOS,s2D_t &SF){
    random_device rd{};
    mt19937_64 G{rd()};
    normal_distribution<double> dis(0.0, 1.0);

    for(int i = 0; i < 2*NUM_SMC; i++){
        for(int j = 0; j < 3; j++){
            SPOS[i][j] = SPOS[i][j] + SF[i][j]*dt + sqrt(2*dt)*dis(G);
        }
    }
}
//-----------------------------------------------------------------------------------------------


/* A function to handle random detachment and then attachment */
void random_detach_and_attach(s2D_t &SPOS, c2D_t &POS, links_t &Slinks, std::array<int,NUM_PARTICLES> &sites_availability) {
    calculate_sites_availability(Slinks, sites_availability); // updating sites_availability

    // Detachment
    random_device rd{};
    mt19937 gen{rd()};
    uniform_int_distribution<> dis(0, NUM_SMC - 1);  // Corrected indexing
    int smc_to_detach = dis(gen) * 2;  // SMCs are indexed in pairs

    sites_availability[Slinks[smc_to_detach]] = 0; // Freeing left site
    sites_availability[Slinks[smc_to_detach + 1]] = 0; // Freeing right site

    // Attachment
    int site_to_attach = ranndom_site_chooser_for_smcs(sites_availability);
    if (site_to_attach != -1) {  // Ensure a valid site was found
        Slinks[smc_to_detach] = site_to_attach;
        Slinks[smc_to_detach + 1] = (site_to_attach + 1) % NUM_PARTICLES;
        SPOS[smc_to_detach] = POS[site_to_attach];
        SPOS[smc_to_detach + 1] = POS[(site_to_attach + 1) % NUM_PARTICLES];

        sites_availability[site_to_attach] = 1;
        sites_availability[(site_to_attach + 1) % NUM_PARTICLES] = 1;
    }
    else{
        // Detach Failed
        sites_availability[Slinks[smc_to_detach]] = 0; // Freeing left site
        sites_availability[Slinks[smc_to_detach + 1]] = 0; // Freeing right site
    }
}


//---------------------------------------------------------------------------------------------------


/*This function is responsible for pulling ring polymer by SMCs*/
void check_proximity_and_take_action(c2D_t &POS, s2D_t &SPOS, links_t &Slinks, links_t &if_Slinks_blocked){

    update_if_Slinks_blocked(if_Slinks_blocked,Slinks);
    //left SMC monomer
    for(int i = 0; i < 2*NUM_SMC; i = i+2){
        // If left site is unblocked
        if(if_Slinks_blocked[i] == 0){
            vec_t rleft_1 = difference(SPOS[i],POS[Slinks[i]]);
            vec_t rleft_2 = difference(SPOS[i],POS[(Slinks[i]-1+NUM_PARTICLES)%NUM_PARTICLES]);
            double rleft1 = magnitude(rleft_1); double rleft2 = magnitude(rleft_2);
            bool change_left = rleft1 > rleft2;
            if(change_left){
                Slinks[i] = (Slinks[i]-1+NUM_PARTICLES)%NUM_PARTICLES;
                check_and_update_blockness(i,Slinks,if_Slinks_blocked);
            }
        }
        // If Right is unblocked
        if(if_Slinks_blocked[i+1] == 0){
            vec_t rright_1 = difference(SPOS[i+1],POS[Slinks[i+1]]);
            vec_t rright_2 = difference(SPOS[i+1],POS[(Slinks[i+1]+1)%NUM_PARTICLES]);
            double rright1 = magnitude(rright_1); double rright2 = magnitude(rright_2);
            bool change_right = rright1 > rright2;
            
            if(change_right){
                Slinks[i+1] = (Slinks[i+1]+1)%NUM_PARTICLES;
                check_and_update_blockness(i+1,Slinks,if_Slinks_blocked);
            }
        }

    }    
}

//--------------------------------------------------------------------------------------------


/* This function saves position data in computer for SMCs*/
void save_SMC_polymer_xyz(std::ofstream& out,c2D_t &POS, s2D_t &SPOS, int step){
    // ofstream out(file_name, ios::app);
    out<<NUM_PARTICLES+2*NUM_SMC<<"\n"<<step<<"\n";
    for(int i = 0; i <NUM_PARTICLES; i++){
        out<<"Xp_"<<i<<" ";
        for(int j = 0; j < 3; j++){
            out<<setprecision(15)<<POS[i][j]<<" ";
        }
        // validate_sim(POS[i],i, step,"Xp_",flag);
        out<<"\n";
    }
    for(int i = 0; i <2*NUM_SMC; i = i+2){
        out<<"Sp_"<<i<<" ";
        for(int j = 0; j < 3; j++){
            out<<setprecision(15)<<SPOS[i][j]<<" ";
        }
        // validate_sim(SPOS[i],i, step,"Sp_",flag);
        out<<"\n";
        out<<"Sp_"<<i+1<<" ";
        for(int j = 0; j < 3; j++){
            out<<setprecision(15)<<SPOS[i+1][j]<<" ";
        }
        // validate_sim(SPOS[i+1],i+1, step,"Sp_",flag);
        out<<"\n";
    }
    // return flag;
}

//----------------------------------------------------------------------------------------

/* This function sotres SMC Slink data in a file*/
void save_SMC_links(std::ofstream& out, links_t &Slinks, int step){
    out<<step<<" ";
    for(int i = 0; i < 2*NUM_SMC; i++){
        out<<Slinks[i]<<" ";
    }
    out<<"\n";
}
//----------------------------------------------------------------------------------------

/* This function writes a configuration file for a running simulation */
void write_config(string configFile , int steps, int equisteps, int interval){
    ofstream config(configFile, ios::out);
    config << "************For Ring Polymer**************"<<endl;
    config << "NOP : "<<NUM_PARTICLES<<endl;
    config << "TEMP : " << TEMPERATURE << endl;
    config << "Ks(Ring Polymer) : " << SPRING_CONSTANT << endl;
    config << "l (Ring Polymer) : " << SPRING_L << endl;
    config << "Ring Diameter : " << RING_DIAMETER << endl;
    config << "Cutoff Radius (LJ) : " << CUTTOFF_RADIUS << endl;
    config << "************For SMCs**************"<<endl;
    config << "NUM_SMC : " << NUM_SMC << endl;
    config << "Ks (b/w SMC and RP) : " << K_SMC_POLY << endl;
    config << "Ks (SMCS) : " << K_SMCs << endl;
    config << "Reeling Force : " << REE_FORCE << endl;
    config << "Dissociation time interval : " << DISSOCIATION_TIME_INTERVAL << endl;
    config << "Dissociation step interval : " << DISSOCIATION_STEP_INTERVAL << endl;
    config << "************For Neighbour List**************"<<endl;
    config << "NL_cutoff_sq : " << NL_cutoff_sq << endl;
    config << "SKIN : " << SKIN << endl;
    config << "************Simulation**************"<<endl;
    config << "dt : " << dt << endl;
    config << "Total Steps : " << steps << endl;
    config << "Equilibrium Steps : " << equisteps <<endl;
    config << "Interval : " << interval << endl;
    config.close();

}

/* This is a debugging function to validate if the simulation is running fine. */
bool validate_sim(c2D_t &POS, s2D_t &SPOS, int &index, string &msg){
    auto CalcCOM = [](vec_t &COM, c2D_t &POS, s2D_t &SPOS) ->void{
        for(int i = 0; i < NUM_PARTICLES; i++){
            for (int j = 0; j < 3; j++){
                COM[j] += POS[i][j];
            }
        }
        for(int i = 0; i < 2*NUM_SMC; i++){
            for(int j = 0; j <3; j++){
                COM[j] += SPOS[i][j];
            }
        }
        for(int j = 0; j < 3; j++){
            COM[j] = COM[j]/(NUM_PARTICLES+2.0f*NUM_SMC);
        }
    };

    vec_t COM = {0,0,0};
    CalcCOM(COM,POS, SPOS);

    for(int i = 0; i < NUM_PARTICLES; i++){
        for(int j = 0; j < 3; j++){
            if(abs(POS[i][j] - COM[j]) > ideal_sim_max_distance){
                index = i; msg = "Xp_";
                return true;
            }
        }
    }

    for(int i = 0; i < 2*NUM_SMC; i++){
        for(int j = 0; j < 3; j++){
            if(abs(SPOS[i][j] - COM[j]) > ideal_sim_max_distance){
                index = i; msg = "Sp_";
                return true;
            }
        }
    }
    return false;
}
//--------------------------------------------------------------------------------------


/* This function was implemented only to do quick debugging, it will be removed in future */
void save_debug_data(links_t &Slinks, links_t &if_Slinks_blocked, s2D_t &SF, c2D_t &F, std::array<int,NUM_PARTICLES> &sites_availability, string fileSMC, string fileC, int step){
    ofstream fSMC(fileSMC, ios::app);
    fSMC<<step<<"\n";
    for(int i = 0; i <2*NUM_SMC; i++){
        fSMC<<i<<" "<<Slinks[i]<<" "<<if_Slinks_blocked[i]<<" ";
        for(int j = 0; j < 3; j++){
            fSMC<<setprecision(15)<<SF[i][j]<<" ";
        }
        fSMC<<"\n";
    }

    ofstream Cfile(fileC, ios::app);
    Cfile<<step<<"\n";
    for(int i = 0; i <NUM_PARTICLES; i++){
        Cfile<<i<<" "<<sites_availability[i]<<" ";
        for(int j = 0; j < 3; j++){
            Cfile<<setprecision(15)<<F[i][j]<<" ";
        }
        Cfile<<"\n";
    }
}
//-----------------------------------------------------------------------------------
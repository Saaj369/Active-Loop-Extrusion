// Extra helper functions supporting the polymer MD simulation.
// Used for specific tasks not covered in helpers.cpp.

# include <iostream>
# include <cmath>
# include <random>
# include <fstream>
# include <cstring>
# include <iomanip>
# include <array>
# include "parameters.h"
# include "helpers.h"
# include "helpers2.h"

using namespace std;
typedef std::array<std::array<double,3>, NUM_PARTICLES> c2D_t;
typedef std::array<int, NUM_PARTICLES*NUM_PARTICLES> NList_t;

void initialise(c2D_t &POS, c2D_t &VEL){
    random_device rd{};
    mt19937 gen{rd()};
	uniform_real_distribution<double> dis(-0.5, 0.5);
	double avgV[3] = {};
	for(int p=0; p<NUM_PARTICLES; p++)
	{
		POS[p][0] = 0.5*RING_DIAMETER*std::cos(2.0*PI*p/NUM_PARTICLES);
		POS[p][1] = 0.5*RING_DIAMETER*std::sin(2.0*PI*p/NUM_PARTICLES);
		POS[p][2] = 0;
		VEL[p][0] = dis(gen)*std::sqrt(12.0); 
		VEL[p][1] = dis(gen)*std::sqrt(12.0); 
		VEL[p][2] = dis(gen)*std::sqrt(12.0);

		avgV[0] += VEL[p][0]; avgV[1] += VEL[p][1]; avgV[2] += VEL[p][2];
	}
	avgV[0] = avgV[0]/NUM_PARTICLES; avgV[1] = avgV[1]/NUM_PARTICLES; avgV[2] = avgV[2]/NUM_PARTICLES;

	for(int p=0; p<NUM_PARTICLES; p++)
	{
		VEL[p][0] -= avgV[0]; VEL[p][1] -= avgV[1]; VEL[p][2] -= avgV[2];
	}
}

void initialiseBD(c2D_t &POS){
	for(int p=0; p<NUM_PARTICLES; p++)
	{
		POS[p][0] = 0.5*RING_DIAMETER*std::cos(2.0*PI*p/NUM_PARTICLES);
		POS[p][1] = 0.5*RING_DIAMETER*std::sin(2.0*PI*p/NUM_PARTICLES);
		POS[p][2] = 0;
	}
}

void eval_force(c2D_t &NF, c2D_t &POS, double &PE){
    fill(&NF[0][0], &NF[0][0] + NUM_PARTICLES * 3, 0);
    PE = 0;
    //Spring Force
    for(int i = 0; i < NUM_PARTICLES; i++){
        array<double,3> r12 = difference(POS[i], POS[(i+1)%NUM_PARTICLES]);
        double r = magnitude(r12);
        array<double,3> Fsp = spring_force(r12,r);
        for(int j = 0; j < 3; j++){
            NF[i][j] += Fsp[j];
            NF[(i+1)%NUM_PARTICLES][j] -= Fsp[j];
        }
        PE += 0.5*(SPRING_CONSTANT/TEMPERATURE)*pow((r - SPRING_L),2.0);
        
    }

    //LJ Force
    for(int i = 0; i < NUM_PARTICLES-1; i++){
        for(int j = i+1; j < NUM_PARTICLES; j++){
            array<double,3> rij = difference(POS[i],POS[j]);
            double r = magnitude(rij);

            if( r <= CUTTOFF_RADIUS){
                array<double,3> Flj = lj_force(rij,r);
                for(int k = 0; k < 3; k++){
                    NF[i][k] += Flj[k];
                    NF[j][k] -= Flj[k];
                }
                PE += (4.0/TEMPERATURE)*(pow(r, -12.0) - pow(r, -6.0)+0.25);
            }
        }
    }
}

void eval_force_using_NL(c2D_t &NF, c2D_t &POS, double &PE, NList_t &NL, array< array<int,2>, NUM_PARTICLES> &NLsize){
    fill(&NF[0][0], &NF[0][0] + NUM_PARTICLES * 3, 0);
    PE = 0;
    //Spring Force
    for(int i = 0; i < NUM_PARTICLES; i++){
        array<double,3> r12 = difference(POS[i], POS[(i+1)%NUM_PARTICLES]);
        double r = magnitude(r12);
        array<double,3> Fsp = spring_force(r12,r);
        for(int j = 0; j < 3; j++){
            NF[i][j] += Fsp[j];
            NF[(i+1)%NUM_PARTICLES][j] -= Fsp[j];
        }
        PE += 0.5*(SPRING_CONSTANT/TEMPERATURE)*pow((r - SPRING_L),2.0);
        
    }

    //LJ Force
    for(int i = 0; i < NUM_PARTICLES-1; i++){
        for(int index = NLsize[i][0]; index < NLsize[i][1]; index++){
            int j = NL[index];
            array<double,3> rij = difference(POS[i],POS[j]);
            double r = magnitude(rij);

            if( r <= CUTTOFF_RADIUS){
                array<double,3> Flj = lj_force(rij,r);
                for(int k = 0; k < 3; k++){
                    NF[i][k] += Flj[k];
                    NF[j][k] -= Flj[k];
                }
                PE += (4.0/TEMPERATURE)*(pow(r, -12.0) - pow(r, -6.0)+0.25);
            }
        }
    }
}

void velocity_verlet_step(c2D_t &POS, c2D_t &VEL, c2D_t &F, c2D_t &NF, double &PE){
    for(int i = 0; i < NUM_PARTICLES; i++){
        for(int j = 0; j < 3; j++){
            POS[i][j] = POS[i][j] + VEL[i][j]*dt + 0.5*F[i][j]*dt*dt;
        }
    }
    eval_force(NF,POS,PE);
    for(int i = 0; i < NUM_PARTICLES; i++){
        for(int j = 0; j < 3; j++){
            VEL[i][j] = VEL[i][j] + 0.5 * dt * (F[i][j]+NF[i][j]);
        }
    }
    swap(F,NF);
}

void velocity_verlet_step_using_NL(c2D_t &POS,c2D_t &POS0, c2D_t &VEL, c2D_t &F, c2D_t &NF, double &PE,NList_t &NL, array< array<int,2>, NUM_PARTICLES> &NLsize){
    for(int i = 0; i < NUM_PARTICLES; i++){
        for(int j = 0; j < 3; j++){
            POS[i][j] = POS[i][j] + VEL[i][j]*dt + 0.5*F[i][j]*dt*dt;
        }
    }
    if(check_if_NList_requires_update(POS0,POS)){
        neighbour_list(NL,NLsize,POS);
        copy_POS_data(POS,POS0);
    }
    eval_force_using_NL(NF,POS,PE,NL,NLsize);
    for(int i = 0; i < NUM_PARTICLES; i++){
        for(int j = 0; j < 3; j++){
            VEL[i][j] = VEL[i][j] + 0.5 * dt * (F[i][j]+NF[i][j]);
        }
    }
    swap(F,NF);
}

void velocity_verlet_BD(c2D_t &POS, c2D_t &F){
    random_device rd{};
    mt19937_64 G{rd()};
    normal_distribution<double> dis(0.0, 1.0);
    // double PE = 0; eval_force(F,POS,PE);
    for(int i = 0; i < NUM_PARTICLES; i++){
        for(int j = 0; j < 3; j++){
            POS[i][j] = POS[i][j] + F[i][j]*dt + sqrt(2*dt)*dis(G);
        }
    }

}

void velocity_verlet_BD_using_NL(c2D_t &POS,c2D_t &F){
    random_device rd{};
    mt19937_64 G{rd()};
    normal_distribution<double> dis(0.0, 1.0);

    // ofstream c_gauss("dataSMC/10_c_gauss.txt",ios :: app);
    // if(check_if_NList_requires_update(POS0,POS)){
    //     neighbour_list(NL,NLsize,POS);
    //     copy_POS_data(POS,POS0);
    // }
    // double PE = 0; eval_force_using_NL(F,POS,PE,NL,NLsize);
    for(int i = 0; i < NUM_PARTICLES; i++){
        for(int j = 0; j < 3; j++){
            // double g = dis(G);
            POS[i][j] = POS[i][j] + F[i][j]*dt + sqrt(2*dt)*dis(G);
            // c_gauss<<setprecision(15)<<g<<" ";
            // POS[i][j] = POS[i][j] + F[i][j]*dt;
        }
        // c_gauss<<"\n";
    }
    // c_gauss<<"-----------------0x0------------------\n";

}

void store_position_data(c2D_t &POS, string file_name, int index){
    ofstream out(file_name,ios :: app);
    out<<NUM_PARTICLES<<"\n"<<index<<"\n";
    for(int i = 0; i < NUM_PARTICLES; i++){
        out<<"Xp_"<<i<<" ";
        for(int j = 0; j < 3; j++){
            out<<setprecision(15)<<POS[i][j]<<" ";
        }
        out<<"\n";
    }
}

void store_KE_PE_momentum(double &PE, c2D_t &VEL,string file_name,int index){
    double KE = 0;
    double momentum[3] = {};
    for(int i = 0; i < NUM_PARTICLES; i++){
        for(int j = 0; j < 3; j++){
            momentum[j] += VEL[i][j];
            KE += VEL[i][j]*VEL[i][j];
        }
    }
    KE = KE*0.5;

    ofstream out(file_name,ios :: app);
    out<<index<<" "<<setprecision(15)<<PE<<" "<<setprecision(15)<<KE<<" "<<setprecision(15)<<momentum[0]<<" "<<setprecision(15)<<momentum[1]<<" "<<setprecision(15)<<momentum[2]<<"\n";
}
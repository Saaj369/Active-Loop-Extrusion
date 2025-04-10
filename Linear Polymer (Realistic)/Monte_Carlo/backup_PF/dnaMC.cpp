#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <array>
#include <algorithm>
#include <omp.h>
#include <unordered_map> // for activeSMCs
// #include "pcg/pcg_random.hpp"
#include <p2rng/p2rng.hpp>
#include "parameters.h"
#include "monomer.h"

using namespace std;

typedef std::unordered_map<int, std::pair<int, int>> activeSMCs_t;
typedef std::array<Monomer, NUM_PARTICLES> polystate_t;

/*
    state = 0 => free monomer
    state = 1 => left smc
    state = 2 => right smc
    state = 3 => CTCF (+)
    state = 4 => CTCF (-)
*/
/******************************************************Binding step ********************************************* */
void try_binding(pcg32 &gen, polystate_t &polystate, activeSMCs_t &activeSMCs, int &smc_id)
{
    uniform_real_distribution<> dis(0, 1);
    for (int i = 0; i < NUM_PARTICLES - 1; i++)
    {
        if (polystate[i].getState() == 0 && polystate[i + 1].getState() == 0)
        {
            float random_num = dis(gen);
            if (random_num <= B_prob)
            {
                // Occupying the site
                polystate[i].setState(1);
                polystate[i + 1].setState(2);
                polystate[i].setIndex(i + 1);
                polystate[i + 1].setIndex(i);
                // binding action
                activeSMCs[smc_id] = {i, i + 1};
                smc_id++;
            }
        }
    }
}

/******************************************************Unbinding (Or) Translocation *******************************/
void try_unbinding_or_translocation(pcg32 &gen2, polystate_t &polystate, activeSMCs_t &activeSMCs, int &smc_id)
{
    uniform_real_distribution<float> rdis2(0, 1);
    int index_R;
    float rand_num;
    int ids_to_erase[smc_id] = {};
    int point = 0;
    for (auto it = activeSMCs.begin(); it != activeSMCs.end();)
    {
        int id_ = it->first;
        int leftIndex = it->second.first;
        int rightIndex = it->second.second;
        rand_num = rdis2(gen2);
        // Checking for unbinding
        if (rand_num <= Un_pob)
        {
            polystate[leftIndex].setState(0);  // set to 'free'
            polystate[rightIndex].setState(0); // Set to 'free'

            // erasing SMCS from active SMCs table
            it = activeSMCs.erase(it);
        }
        else
        {
            // Unbinding failed, translocation
            int l_state = leftIndex > 1 ? polystate[leftIndex - 1].getState() : -1;                    // state at leftindex - 1
            int l__state = leftIndex > 2 ? polystate[leftIndex - 2].getState() : -1;                   // state at leftindex - 2
            int r_state = rightIndex < NUM_PARTICLES - 1 ? polystate[rightIndex + 1].getState() : -1;  // state at rightindex + 1
            int r__state = rightIndex < NUM_PARTICLES - 2 ? polystate[rightIndex + 2].getState() : -1; // state at rightindex + 2
            //--------------------------------------------------------------------------
            // Left-side logic
            rand_num = rdis2(gen2);
            if (l_state == 0 && rand_num < Free_jump_prob)
            { // Free state jump
                polystate[leftIndex - 1] = polystate[leftIndex];
                polystate[leftIndex].setState(0);
                activeSMCs[id_].first = leftIndex - 1;
            }
            else if (l_state == 4 && l__state == 0 && rand_num < Free_jump_prob)
            { // (-) CTCF and free state jump
                polystate[leftIndex - 2] = polystate[leftIndex];
                polystate[leftIndex].setState(0);
                activeSMCs[id_].first = leftIndex - 2;
            }

            // Right-side logic
            rand_num = rdis2(gen2);
            if (r_state == 0 && rand_num < Free_jump_prob)
            { // Free state jump
                polystate[rightIndex + 1] = polystate[rightIndex];
                polystate[rightIndex].setState(0);
                activeSMCs[id_].second = rightIndex + 1;
            }
            else if (r_state == 3 && r__state == 0 && rand_num < Free_jump_prob)
            { // (-) CTCF and free state jump
                polystate[rightIndex + 2] = polystate[rightIndex];
                polystate[rightIndex].setState(0);
                activeSMCs[id_].second = rightIndex + 2;
            }

            it++; // next iterator
        }
    }
}

/*****************************************************Monte Carlo step*********************************************/
void execute_MonteCarloStep(pcg32 &gen, pcg32 &bind_G, pcg32 &unbind_G, polystate_t &polystate, activeSMCs_t &activeSMCs, int &smc_id)
{
    uniform_real_distribution<float> rdis(0, 1);
    float choice = rdis(gen);

    if (choice <= 0.5)
    {
        // binding then unbinding, translocation
        try_binding(bind_G, polystate, activeSMCs, smc_id);
        try_unbinding_or_translocation(unbind_G, polystate, activeSMCs, smc_id);
    }
    else
    {
        try_unbinding_or_translocation(unbind_G, polystate, activeSMCs, smc_id);
        try_binding(bind_G, polystate, activeSMCs, smc_id);
    }
}

/****************************************************Displaying Progress in Terminal**************************** */
void display_progress(int totalSteps, int finishedSteps)
{
    cout << "\033[2J\033[H";
    cout.flush();
    // Display progress
    cout << finishedSteps << "/" << totalSteps << "------->" << 100.0f * finishedSteps / totalSteps << " %" << "\r";
    cout.flush();
}

/***************************************************Save Polystate and Active SMCs (for debugging)************** */
void save_system_state(activeSMCs_t &activeSMCs, std::ofstream &out, int step)
{
    // saving active SMCs
    out << step << "\n";
    out << activeSMCs.size() << "\n";
    for (auto it = activeSMCs.begin(); it != activeSMCs.end();)
    {
        int id_ = it->first;
        int leftIndex = it->second.first;
        int rightIndex = it->second.second;
        out << id_ << " " << leftIndex << " " << rightIndex << "\n";
        it++;
    }
}

//***************************************************Loading CTCF Data***************************************** */
void load_ctcf(string fileName, polystate_t &polystate)
{
    std::ifstream file(fileName);
    if (!file.is_open())
    {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // Define the 2D vector to store data
    std::vector<std::vector<int>> data;

    std::string line;

    // Skip the header line
    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        int coarseIndex, numPeaks;
        char orientation;
        char comma;

        // Read and parse each line
        ss >> coarseIndex >> comma >> numPeaks >> comma >> orientation;

        // Calculate numPeaks based on orientation
        numPeaks = (orientation == '-') ? -numPeaks : numPeaks;

        // Store in a vector (or array) as per requirement
        data.push_back({coarseIndex, numPeaks});
    }

    file.close();

    // prepare polymerState
    for (const auto &row : data)
    {
        if (row[1] > 0)
        {
            if (row[0] < NUM_PARTICLES)
            {
                polystate[row[0]].setState(3); // (+) oriented CTCF
            }
        }
        else
        {
            if (row[0] < NUM_PARTICLES)
            {
                polystate[row[0]].setState(4); // (-) oriented CTCF
            }
        }
    }
}

//******************************************* Checking Is TAD fully Extruded ***********************************
bool isFullyExtruded(activeSMCs_t &activeSMCs, int Ni, int Nf){
    for (auto it = activeSMCs.begin(); it != activeSMCs.end();)
    {
        int id_ = it->first;
        int leftIndex = it->second.first;
        int rightIndex = it->second.second;

        if(leftIndex == Ni + 1 & rightIndex == Nf - 1){
            return true;
        }
        it++;
    }
    return false;
}

void integrate_and_simulate(int steps, int saveFreq)
{

    // ------------------------------------- Container Arrays -------------------------------------------

    // for(int i = 0; i < NUM_PARTICLES; i++){
    //     polystate[i].setState(0);               // State = 0
    // }

    //--------------------------------------Files and Data locations-------------------------------------
    string stateFile = "dataDNA/smc_data_7FEB.txt";
    string CTCFfile = "CTCF_coarse_positions.csv";

    // -------------------------------------Hints for simulation ---------------------------------------
    string msg_init = "(Equilibrating)";
    string msg_simulate = "(Simulation Running)";

    // ------------------------------------Seeded Random Engines---------------------------------------
    pcg32 binding_GS{binding_RAND_SEED}; // binding
    pcg32 unbind_GS{unbind_RAND_SEED};   // unbinding and translocation
    pcg32 choice_GS{choice_RAND_SEED};   // choice of action

    // -----------------------------------Not seeded Random Engines------------------------------------
    pcg32 binding_G{}; // binding
    pcg32 unbind_G{};  // unbinding and translocation
    pcg32 choice_G{};  // choice of action

    omp_set_num_threads(5);

    // ------------------------------------------Loading CTCF----------------------------------------
    // load_ctcf(CTCFfile, polystate);
    // polystate[100].setState(3); // (+) oriented CTCF
    // polystate[300].setState(4); // (-) oriented CTCF

    ofstream out("dataDNA/extrusion_probabilities_40.txt");
    if (!out) {
        cerr << "Error: Could not open file for writing." << endl;
        exit(1);
    }

    #pragma omp parallel for
    for (int tad = 0; tad < 20; tad++)
    {
        polystate_t polystate;   // Polymer state (0(F), 1(LS), 2(RS), 3(+CTCF), 4(-CTCF))
        activeSMCs_t activeSMCs; // Active SMCs
        int smc_id = 0;          // unique id for smc

        int Ni = (NUM_PARTICLES / 2) - (TAD[tad] / 2);
        int Nf = Ni + TAD[tad];

        for (int i = 0; i < NUM_PARTICLES; i++)
        {
            polystate[i].setState(0); // State = 0
        }

        polystate[Ni].setState(3); // (+) oriented CTCF
        polystate[Nf].setState(4); // (-) oriented CTCF

        // string stateF = "dataDNA/smc_data_13FEB_47_" + to_string(TAD[tad]) + ".txt";
        // // ---------------------------------Opening POS File---------------------------------------------
        // ofstream out(stateF);
        // if (!out)
        // {
        //     cerr << "Error: Could not open stateFile for writing." << endl;
        //     exit(1);
        // }

        // ---------------------------------Simulating Monte Carlo-------------------------------

        bool prev_ext_check = false;
        bool curr_ext_check = false;
        int full_extrusion_time = 0;
        int total_time = 0;
        cout << "Equilibriating for => " << TAD[tad] << endl;
        for (int i = 0; i < 1000; i++)
        {
            execute_MonteCarloStep(choice_G, binding_G, unbind_G, polystate, activeSMCs, smc_id);

        }
        cout << "Simulating for => " << TAD[tad] << endl;
        for (int i = 0; i < steps; i++)
        {
            if (i % saveFreq == 0)
            {
                display_progress(steps, i);
                // save_system_state(activeSMCs, out, i);
            }
            execute_MonteCarloStep(choice_G, binding_G, unbind_G, polystate, activeSMCs, smc_id);
            curr_ext_check = isFullyExtruded(activeSMCs, Ni, Nf);
            if(curr_ext_check & prev_ext_check){
                full_extrusion_time++;
            }
            prev_ext_check = curr_ext_check;
            total_time++;
        }

        double full_extrusion_probability = full_extrusion_time / (total_time-1);
         #pragma omp critical
        {
            // Writing the result to the file (ensuring no race condition)
            out << "Thread " << omp_get_thread_num() << " " << TAD[tad]<< " " << full_extrusion_time << " " << total_time << " " << setprecision(7) <<full_extrusion_probability << endl;
        }
        // out.close();
    }

    out.close();

    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
}

int main()
{
    int steps = 1e6;
    int saveFreq = 1e4;

    integrate_and_simulate(steps, saveFreq);
    return 0;
}
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
# include <p2rng/p2rng.hpp>
#include "parameters.h"
#include "monomer.h"

using namespace std;
// typedef std::array<double, 3> vec_t;
// typedef std::array<int, NUM_PARTICLES * NUM_PARTICLES> NList_t;
// typedef std::array<std::array<int, 2>, NUM_PARTICLES> NLindex2d_t;
// typedef std::array<std::array<double, 3>, NUM_PARTICLES> c2D_t;
typedef double vec_t[3];                                                // Raw
typedef double c2D_t[NUM_PARTICLES][3];                                 // Raw
typedef double g2D_t[6*NUM_PARTICLES];                               // Raw
// typedef double g2D_t[2*NUM_PARTICLES][3];

typedef std::unordered_map<int, std::pair<int, int>> activeSMCs_t;
typedef std::array<Monomer, NUM_PARTICLES> polystate_t;




/*******************************************************Linear Chain Initialisation*************************************/
void initialiseBD(c2D_t &POS)
{
    for (int p = 0; p < NUM_PARTICLES; p++)
    {
        POS[p][0] = SPRING_L * p;
        POS[p][1] = 0;
        POS[p][2] = 0;
    }
}

/*******************************************************Load Saved or Equilibrium state*********************************/
void LoadEquiBD(c2D_t &POS, polystate_t &polystate, activeSMCs_t &activeSMCs, int &smc_id){
    // -----------------------------------------------Loading position data--------------------------------------------------
    std::string filename = "eConfig.xyz";
    std::ifstream infile(filename);
    
    if (!infile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
    }


    // Skip the first line (number of particles)
    std::string line;
    std::getline(infile, line);

    // Skip the second line (comment)
    std::getline(infile, line);

    // Read the particle data (positions)
    int index = 0;
    while (std::getline(infile, line) && index < NUM_PARTICLES) {
        std::istringstream ss(line);
        
        // Skip particle name and color (first 4 columns)
        std::string particle_name;
        int r, g, b;
        ss >> particle_name >> r >> g >> b;
        
        // Read the position coordinates (5th, 6th, and 7th columns)
        double x, y, z;
        ss >> x >> y >> z;

        // Store the position in the array
        POS[index][0] = x;
        POS[index][1] = y;
        POS[index][2] = z;

        ++index;
    }

    // Close the file
    infile.close();

    // ---------------------------------Loading SMC data---------------------------------------------------

    std::string filename_ = "eConfigSMC.txt";  // Specify your file path here
    std::ifstream infile_(filename_);
    
    if (!infile_.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        // return 1;
    }

    // Define a vector to store the tuples
    std::vector<std::array<int, 3>> data;

    while (std::getline(infile_, line)) {
        // Replace commas with spaces to simplify parsing
        std::replace(line.begin(), line.end(), ',', ' ');

        std::istringstream ss(line);
        int x, y, z;
        char paren_open, paren_close;

        // Loop over the space-separated values in the line
        while (ss >> paren_open >> x >> y >> z >> paren_close) {
            // Store the tuple (x, y, z) into the data vector
            data.push_back({x, y, z});
        }
    }

    // Close the file
    infile_.close();

    // Occupying the site
    for(int i = 0; i < data.size(); ++i)
    {   
        int l_index = data[i][0];
        int r_index = data[i][2];
        polystate[l_index].setState(1);
        polystate[r_index].setState(2);
        polystate[l_index].setIndex(r_index);
        polystate[r_index].setIndex(l_index);
        // binding action
        activeSMCs[smc_id] = {l_index, r_index};
        smc_id++;
        }
}

/*******************************************************Force Calculation(NL)*******************************************/
void eval_force_using_NL(c2D_t &NF, c2D_t &POS)
{
    // fill(&NF[0][0], &NF[0][0] + NUM_PARTICLES * 3, 0);
    std::memset(&NF[0][0], 0, NUM_PARTICLES * 3 * sizeof(double));

    const int NUM_THREADS = omp_get_max_threads();
    // std::vector<c2D_t> localForces(NUM_THREADS, {}); // Local force containers for each thread
    double localForces[NUM_THREADS][NUM_PARTICLES][3] = {};


    // Parallelized Spring Force
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        c2D_t &localNF = localForces[thread_id];

        #pragma omp for
        for (int i = 0; i < NUM_PARTICLES - 1; i++)
        {
            array<double, 3> r12 = {};
            double rSquared = 0.0;

            // Compute the displacement and squared distance
            #pragma omp simd reduction(+:rSquared)
            for (int m = 0; m < 3; m++)
            {
                r12[m] = POS[i][m] - POS[i + 1][m];
                rSquared += r12[m] * r12[m];
            }

            double r = sqrt(rSquared);
            double factor = -SPRING_CONSTANT * (1.0 - 1.0 / r);

            // Compute the spring force
            array<double, 3> Fsp = {};
            #pragma omp simd
            for (int j = 0; j < 3; j++)
            {
                Fsp[j] = factor * r12[j];
            }

            // Update local forces
            for (int j = 0; j < 3; j++)
            {
                localNF[i][j] += Fsp[j];
                localNF[i + 1][j] -= Fsp[j];
            }
        }
    }


    // Merge local forces into the global NF array
    #pragma omp parallel for
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int t = 0; t < NUM_THREADS; t++)
            {
                NF[i][j] += localForces[t][i][j];
            }
        }
    }
}

/*******************************************************SMC Force Injection**************************************** */
void smc_bondforce_injection(c2D_t &NF, c2D_t &POS, activeSMCs_t &activeSMCs)
{   
    // #pragma omp parallel for
    for (auto it = activeSMCs.begin(); it != activeSMCs.end(); ++it)
    {
        int id_ = it->first;
        int leftIndex = it->second.first;
        int rightIndex = it->second.second;

        array<double, 3> r12 = {};
        double r = 0;
        for (int m = 0; m < 3; m++)
        {
            r12[m] = POS[leftIndex][m] - POS[rightIndex][m];
            r += r12[m] * r12[m];
        }
        r = sqrt(r);

        array<double, 3> Fsp = {};
        double factor = -K_SMC_POLY * (1.0 - 1.0/r);
        for (int j = 0; j < 3; j++)
        {
            Fsp[j] = factor * r12[j];
        }
        for (int j = 0; j < 3; j++)
        {
            // #pragma omp atomic
            NF[leftIndex][j] += Fsp[j];
            // #pragma omp atomic
            NF[rightIndex][j] -= Fsp[j];
        }
    }
}

/*******************************************************Integration Step****************************************** */
void underdamped_velocity_verlet_BD_using_NL(g2D_t &GRAN, c2D_t &POS, c2D_t &VEL, c2D_t &F, activeSMCs_t &activeSMCs)
{
    // normal_distribution<double> dis(0.0, 1.0);
    // using uniform real distribution in range of -1 to 1
    // uniform_real_distribution<double> dis(-1.0,1.0);
    int rand;
    #pragma omp parallel for
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            VEL[i][j] = VEL[i][j] + dth * F[i][j] + B * GRAN[i*3 + j]; // dth = 0.5 * dt
            POS[i][j] = POS[i][j] +  C * VEL[i][j]; // C = (2 * dt / (2 + dt) )
        }
    }

    eval_force_using_NL(F, POS);
    smc_bondforce_injection(F, POS, activeSMCs);

    #pragma omp parallel for
    for(int i = 0; i < NUM_PARTICLES; i++){
        for(int j = 0; j < 3; j++){
            VEL[i][j] = A * VEL[i][j] + B * GRAN[(i + NUM_PARTICLES)*3 + j] + dth * F[i][j]; 
            // A = ((2 - dt) / (2 + dt))
            // B = sqrt(0.5 * dt)
        }
    }

}

/******************************************************Saving Position Data****************************************/
void save_polymer_xyz(std::ofstream &out, c2D_t &POS, int step, polystate_t &polystate)
{
    // ^^^^^^^^^^^^^^^^^^^^^^^^Optimised^^^^^^^^^^^^^^^^^^^^^^^^^

    // state strings for all particles
    std::vector<std::string> state_strings(NUM_PARTICLES);
    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        int state = polystate[i].getState();
        state_strings[i] = "X_" + std::to_string(i) + " " + color_code[state];
    }

    out << NUM_PARTICLES << "\n" << step << "\n";

    std::string buffer;
    buffer.reserve(40011);  // Reserve space in the string (40kb for 1000 particles)

    for (int i = 0; i < NUM_PARTICLES; i++)
    {
        buffer += state_strings[i] + " ";

        std::ostringstream position_stream;
        position_stream << std::fixed << std::setprecision(7);
        for (int j = 0; j < 3; j++)
        {
            position_stream << POS[i][j] << " ";
        }

        buffer += position_stream.str();

        buffer += "\n";
    }

    // Write all the accumulated data to the file in one go
    out << buffer;
}

void save_polymer_binary_optimized(std::ofstream &out, c2D_t &POS, int step, polystate_t &polystate)
{


    // Calculate total buffer size: 19 bytes per particle * 1000 particles
    char* buffer = new char[BUFFER_SIZE];  // Allocate memory for the buffer

    char* buffer_ptr = buffer;  // Pointer to the current position in the buffer

    // Write number of particles (int) at the start
    std::memcpy(buffer_ptr, &NUM_PARTICLES, sizeof(NUM_PARTICLES));
    buffer_ptr += sizeof(NUM_PARTICLES);  // Move pointer forward

    // Iterate over each particle and write its data
    for (int i = 0; i < NUM_PARTICLES; ++i)
    {
        // Write particle ID (int)
        int particle_id = i;
        std::memcpy(buffer_ptr, &particle_id, sizeof(particle_id));
        buffer_ptr += sizeof(particle_id);

        // Write 3D position (3 floats)
        for (int j = 0; j < 3; ++j)
        {
            float position = static_cast<float>(POS[i][j]);
            std::memcpy(buffer_ptr, &position, sizeof(position));
            buffer_ptr += sizeof(position);
        }

        // Write color (3 unsigned bytes for RGB)
        int state = polystate[i].getState();  // Get particle state (color index)
        int color_code_idx = state; // Example: Red (0), Blue (1), etc.
        uint8_t r, g, b;
        sscanf(color_code[color_code_idx].c_str(), "%hhu %hhu %hhu", &r, &g, &b);

        std::memcpy(buffer_ptr, &r, sizeof(r));
        buffer_ptr += sizeof(r);
        std::memcpy(buffer_ptr, &g, sizeof(g));
        buffer_ptr += sizeof(g);
        std::memcpy(buffer_ptr, &b, sizeof(b));
        buffer_ptr += sizeof(b);
    }

    // Write the entire buffer to the file in one operation
    out.write(buffer, BUFFER_SIZE);

    delete[] buffer;
}

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
            if (random_num <= B_prob) //&& smc_id < 1
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
            rand_num =  rdis2(gen2);
            if (l_state == 0 && rand_num < Free_jump_prob) { // Free state jump
                polystate[leftIndex - 1] = polystate[leftIndex];
                polystate[leftIndex].setState(0);
                activeSMCs[id_].first = leftIndex - 1;
            }
            else if (l_state == 4 && l__state == 0 && rand_num < Free_jump_prob) { // (-) CTCF and free state jump
                polystate[leftIndex - 2] = polystate[leftIndex];
                polystate[leftIndex].setState(0);
                activeSMCs[id_].first = leftIndex - 2;
            }

            // Right-side logic
            rand_num =  rdis2(gen2);
            if (r_state == 0 && rand_num < Free_jump_prob) { // Free state jump
                polystate[rightIndex + 1] = polystate[rightIndex];
                polystate[rightIndex].setState(0);
                activeSMCs[id_].second = rightIndex + 1;
            }
            else if (r_state == 3 && r__state == 0 && rand_num < Free_jump_prob) { // (-) CTCF and free state jump
                polystate[rightIndex + 2] = polystate[rightIndex];
                polystate[rightIndex].setState(0);
                activeSMCs[id_].second = rightIndex + 2;
            }

            it++; // next iterator
        }
    }
}

/*****************************************************Monte Carlo step*********************************************/
void execute_MonteCarloStep(pcg32 &gen, pcg32 &bind_G, pcg32 &unbind_G, polystate_t &polystate, activeSMCs_t &activeSMCs, int &smc_id){
    uniform_real_distribution<float> rdis(0,1);
    float choice = rdis(gen);

    if(choice <= 0.5){
        // binding then unbinding, translocation
        try_binding(bind_G, polystate, activeSMCs, smc_id);
        try_unbinding_or_translocation(unbind_G, polystate, activeSMCs, smc_id);
    }
    else{
        try_unbinding_or_translocation(unbind_G, polystate, activeSMCs, smc_id);
        try_binding(bind_G, polystate, activeSMCs, smc_id);
    }

}

/****************************************************Displaying Progress in Terminal**************************** */
void display_progress(int totalSteps, int finishedSteps, string msg)
{
    cout << "\033[2J\033[H";
    cout.flush();
    // Display progress
    cout << finishedSteps << "/" << totalSteps << "------->" << msg << 100.0f * finishedSteps / totalSteps << " %" << "\r";
    cout.flush();
}

/***************************************************Save Polystate and Active SMCs (for debugging)************** */
void save_system_state(polystate_t &polystate, activeSMCs_t &activeSMCs, string stateFile, int step){
    ofstream out(stateFile, ios::app);
    out << step << "\n";
    //saving state
    for(int i = 0; i < NUM_PARTICLES; i++){
        out << polystate[i].getState() << ", ";
    }
    out << "\n";

    //saving active SMCs
    for(auto it = activeSMCs.begin(); it != activeSMCs.end();){
        int id_ = it -> first;
        int leftIndex = it -> second.first;
        int rightIndex = it -> second.second;
        out << "( " << leftIndex << ", " << id_ << ", " << rightIndex << " ), ";
        it++;
    }
    out << "\n";

    out.close();
}


//***************************************************Loading CTCF Data***************************************** */
void load_ctcf(string fileName, polystate_t &polystate){
    std::ifstream file(fileName);
    if (!file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // Define the 2D vector to store data
    std::vector<std::vector<int>> data;
    
    std::string line;
    
    // Skip the header line
    std::getline(file, line);
    
    while (std::getline(file, line)) {
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

    //prepare polymerState
    for (const auto& row : data) {
        if(row[1] > 0){
            if(row[0]<NUM_PARTICLES){
                polystate[row[0]].setState(3); // (+) oriented CTCF
            }
        }
        else{
            if(row[0] < NUM_PARTICLES){
                polystate[row[0]].setState(4); // (-) oriented CTCF
            }
        }
    }
    
}

inline void gen_normal_random(pcg32 &G, g2D_t &GRAN, int step) {
    p2rng::generate_n(
        GRAN,
        6*NUM_PARTICLES,
        p2rng::bind(trng::normal_dist(0.0, 1.0), pcg32(3*step))
    );
    // std::normal_distribution<double> dis(0.0, 1.0);
    // #pragma omp parallel for
    // for (int i = 0; i < 6*NUM_PARTICLES; ++i) {
    //     GRAN[i] = dis(G);
    // }
}

void integrate_and_simulate(float simDuration, float esimDuration, float eMCDuration, float freq){

    // --------------------------------------- Calculating steps from time------------------------------
    long long steps = simDuration/dt;
    long long equiSteps = esimDuration/dt;
    long long eMCSteps = eMCDuration/dt;
    int data_entry_interval = freq/dt;

    // ------------------------------------- Container Arrays -------------------------------------------
    c2D_t POS = {};                             // Stores positions
    c2D_t POS0 = {};                            // Previous positions for NL
    c2D_t VEL = {};                            // Velocity
    g2D_t GRAN = {};                            // Gaussian Random Numbers
    c2D_t F = {};                               // Forces

    double PE = 0;                              // PE

    polystate_t polystate;                      // Polymer state (0(F), 1(LS), 2(RS), 3(+CTCF), 4(-CTCF))
    activeSMCs_t activeSMCs;                    // Active SMCs
    int smc_id = 0;                             // unique id for smc

    for(int i = 0; i < NUM_PARTICLES; i++){
        polystate[i].setState(0);               // State = 0
    }

    //--------------------------------------Files and Data locations-------------------------------------
    string POSfile = "dataDNA/DNA_21FEB_real_nlj.xyz";
    // string stateFile = "dataDNA/polymer_state_7DEC3.txt";
    string CTCFfile = "CTCF_coarse_positions.csv";
    string EPOSfile = "dataDNA/EDNA_21FEB_NLJ.xyz";
    string MCPOSfile = "dataDNA/MCPOS_21FEB_nlj.xyz";

    // -------------------------------------Hints for simulation ---------------------------------------
    string msg_init = "(Equilibrating)";
    string msg_simulate = "(Simulation Running)";

    // ------------------------------------Seeded Random Engines---------------------------------------
    pcg32 BD_G{BD_RAND_SEED};             // brownian dynamics verlet
    pcg32 binding_G{binding_RAND_SEED};   // binding
    pcg32 unbind_G{unbind_RAND_SEED};     // unbinding and translocation
    pcg32 choice_G{choice_RAND_SEED};     // choice of action

    // -----------------------------------Not seeded Random Engines------------------------------------
    pcg32 BD_GS{};                          // brownian dynamics verlet
    pcg32 binding_GS{};                     // binding
    pcg32 unbind_GS{};                      // unbinding and translocation
    pcg32 choice_GS{};                      // choice of action

    omp_set_num_threads(10);

    //-----------------------------------Preparing Simulation-----------------------------------------
    initialiseBD(POS);                                          // Straight Linear polymer chain
    // LoadEquiBD(POS, polystate, activeSMCs, smc_id);          // Load saved configuration
    eval_force_using_NL(F, POS);                    // Evaluate forces

    ofstream Eout(EPOSfile);
    if (!Eout) {
        cerr << "Error: Could not open POSfile for writing." << endl;
        exit(1);
    }

    // ---------------------------------Equilibriating the Polymer (MD)-------------------------------
    cout<< "Equilibriating" << endl;
    for(int i = 0; i < equiSteps; i++){
        if(i%data_entry_interval == 0){
            save_polymer_xyz(Eout, POS, i, polystate);
            // save_polymer_binary_optimized(DEout, POS, i, polystate);
            // display_progress(steps, i, msg_simulate);
            cout << i << "/" << equiSteps << endl;
        }
        gen_normal_random(BD_G, GRAN, i);
        underdamped_velocity_verlet_BD_using_NL(GRAN, POS, VEL, F, activeSMCs);

    }
    Eout.close();

    
    // ------------------------------------------Loading CTCF----------------------------------------
    load_ctcf(CTCFfile, polystate);
    // polystate[100].setState(3); // (+) oriented CTCF
    // polystate[300].setState(4); // (-) oriented CTCF

    ofstream MCout(MCPOSfile);
    if (!MCout) {
        cerr << "Error: Could not open POSfile for writing." << endl;
        exit(1);
    }

    //------------------------------- Equilibriating MC --------------------------------------------
    cout<< "MC Equilibriating" << endl;
    for(int i = 0; i < eMCSteps; i++){
        if(i%data_entry_interval == 0){
            save_polymer_xyz(MCout, POS, i, polystate);
            // save_polymer_binary_optimized(DEout, POS, i, polystate);
            // display_progress(steps, i, msg_simulate);
            cout << i << "/" << eMCSteps << endl;
        }

        if(i%MC_interval == 0){
            execute_MonteCarloStep(choice_G, binding_G, unbind_G, polystate, activeSMCs, smc_id);
        }
        gen_normal_random(BD_G, GRAN, i);
        underdamped_velocity_verlet_BD_using_NL(GRAN, POS, VEL, F, activeSMCs);

    }
    MCout.close();


    // ---------------------------------Opening POS File---------------------------------------------
    ofstream out(POSfile);
    if (!out) {
        cerr << "Error: Could not open POSfile for writing." << endl;
        exit(1);
    }


    // --------------------------------------------MC-MD Simulation ----------------------------------
    cout << "Simulating" << endl;
    for(int i = 0; i < steps; i++){
        if(i%data_entry_interval == 0){
            save_polymer_xyz(out, POS, i, polystate);
            // save_polymer_binary_optimized(Dout, POS, i, polystate);
            // display_progress(steps, i, msg_simulate);
            cout << i << "/" << steps << endl;
        }

        if(i%MC_interval == 0){
            execute_MonteCarloStep(choice_G, binding_G, unbind_G, polystate, activeSMCs, smc_id);
        }
        gen_normal_random(BD_G, GRAN, i);
        underdamped_velocity_verlet_BD_using_NL(GRAN, POS, VEL, F, activeSMCs);
    }

    out.close();
    cout << "Simulation finished" << endl;

    // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
}

int main(){
    float simDuration =  800000; //100000;       // 25000 Full simulation time (2e7)
    float esimDuration =  200000; //5000;      // Equilibration time (1e6)
    float eMCDuration =  1000*10; //10000;
    float freq = 100;   //10*dtau;                  // Data saving frequency

    integrate_and_simulate(simDuration, esimDuration, eMCDuration, freq);
    return 0;
}
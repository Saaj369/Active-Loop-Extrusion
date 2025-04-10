// This file uses parallelised code to calculate HiC matrix from the .xyz position file

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <sstream>
#include <iomanip>
#include <array>
#include <string>

struct Particle {
    double x, y, z;
};

void processFrame(std::ifstream& file, int numParticles, int frameNumber, const double cutoff, std::array<std::array<int, 1000>, 1000> &hicM) {
    // Read the positions of particles
    std::vector<Particle> particles(numParticles);
    std::string line;
    for (int i = 0; i < numParticles; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        std::string dummy;
        Particle p;
        iss >> dummy >> dummy >> dummy >> dummy >> p.x >> p.y >> p.z;
        particles[i] = p;
    }

    // Calculate HiC using OpenMP
    
    #pragma omp parallel for
    for(int i = 0; i < numParticles; i++){
        for(int j = i+1; j < numParticles; j++){
            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;
            // double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            double dist = dx * dx + dy * dy + dz * dz;
            if(dist < cutoff){
                hicM[i][j] += 1;
                hicM[j][i] += 1;
            }
            hicM[i][i] += 1;
        }
    }
}
void saveMatrix(const std::array<std::array<int, 1000>, 1000>& hicMatrix, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Unable to open file " << filename << " for writing.\n";
        return;
    }

    for (const auto& row : hicMatrix) {
        for (const auto& elem : row) {
            outFile << elem << " "; // Write each element separated by a space
        }
        outFile << "\n"; // End each row with a newline
    }

    outFile.close();
    std::cout << "Matrix saved successfully to " << filename << "\n";
}

int main() {
    omp_set_num_threads(10);
    std::string fileName = "DNA_18DEC_full.xyz";
    std::ifstream file(fileName);
    if (!file) {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        return 1;
    }

    int numParticles = 1000;
    std::string line;
    int frameNumber = 0;
    long long totalFrames = 0;
    long long totalLines = 0;
    double cutoff = 1.5;
    std::array<std::array<int, 1000>, 1000> hicMatrix = {};
    std::string hicFIle = "hic2_.txt";


    // Calculate total frames for progress estimation
    {
        std::ifstream tempFile(fileName);
        while (std::getline(tempFile, line)) {
            totalLines++;
        }
        totalFrames = totalLines / (1000 + 2); // 1000 particles + 2 header lines
    }

    file.clear();
    file.seekg(0);

    while (std::getline(file, line)) {
        if(frameNumber == totalFrames){
            break;
        }
        numParticles = std::stoi(line); // First line: Number of particles
        std::getline(file, line);      // Second line: Step number
        processFrame(file, numParticles, frameNumber, cutoff, hicMatrix);
        frameNumber++;

        // Display progress
        if (frameNumber % 100 == 0 || frameNumber == totalFrames) {
            std::cout << "Progress: " << std::fixed << std::setprecision(2)
                      << (frameNumber * 100.0 / totalFrames) << "%\n";
        }
    }

    file.close();
    saveMatrix(hicMatrix, hicFIle);
    return 0;
}

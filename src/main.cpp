#include <iostream>
#include <cmath>

#include "Laplacian.h"
#include "LatticeMesh.h"
#include "LaplacianWithTensor.h"

#define JACOBI_SOLVER "JacobiSolver.usda"
//#define CG_SOLVER "ConjugateGradient.usda"

template<class T, int latticeSize>
void generateLattice(std::array<std::array<T, latticeSize>, latticeSize> &lattice, std::array<std::array<uint8_t, latticeSize>, latticeSize> &nodeType) {
    for(int i = 0; i < latticeSize; i++) {
        for(int j = 0; j < latticeSize; j++) {
            lattice[i][j] = 0.0;
            if(i< 4 && j < 4) {
                nodeType[i][j] = 1;
            }
            else if( i == 4 && j <= 4) {
                nodeType[i][j] = 2;
            }
            else if( i <= 4 && j == 4) {
                nodeType[i][j] = 2;
            }
            else if( (i > 4 && j == 0) || (i == 0 && j > 4) ) {
                nodeType[i][j] = 2;
            }
            else if (i == latticeSize - 1 || j == latticeSize - 1) {
                nodeType[i][j] = 2;
            }
            else {
                nodeType[i][j] = 0;
            }
            /*
            if(i == 0 || j == 0 || i == latticeSize-1 || j == latticeSize-1) {
                nodeType[i][j] = 2;
            }
            else {
                nodeType[i][j] = 0;
            }*/
        }
    }
}

int main(int argc, char *argv[])
{
    const int latticeSize = 15;
    std::array<std::array<float, latticeSize>, latticeSize> lattice;
    std::array<std::array<uint8_t, latticeSize>, latticeSize> nodeType;
    
    generateLattice<float, latticeSize>(lattice, nodeType);
    
    //Laplacian lattice
    Laplacian2D<float, latticeSize> delOperator;
    
    //To visualise it.
    LatticeMesh<float, latticeSize> simulationMesh;
    simulationMesh.m_cellSize = { latticeSize-1, latticeSize-1};
    simulationMesh.m_gridDX = 0.25;
    simulationMesh.m_nFrames = 1;
    
    std::function<void(float [][latticeSize])> f = [&] (float inputLattice[][latticeSize])
    {
        simulationMesh.regFrame(inputLattice);
    };
    
    delOperator.initialiseLattice(lattice, nodeType);
    delOperator.registerWriteFunctionCallback(f);
    
    
#ifdef JACOBI_SOLVER
    simulationMesh.initialize(delOperator.m_lattice, JACOBI_SOLVER);
    simulationMesh.writeFrame(0);
    delOperator.JacobiIterationSolver();
#else
    simulationMesh.initialize(delOperator.m_lattice, CG_SOLVER);
    simulationMesh.writeFrame(0);
    delOperator.conjugateGradientsSolver();
#endif
    // Write the entire timeline to USD
    simulationMesh.writeUSD();

    return 0;
}


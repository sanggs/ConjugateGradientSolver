#include <iostream>
#include <cmath>

#include "Laplacian.h"
#include "LatticeMesh.h"

int main(int argc, char *argv[])
{
    const int latticeSize = 15;
    Laplacian2D<latticeSize> delOperator;
    
    //To visualise it.
    LatticeMesh<float, latticeSize> simulationMesh;
    simulationMesh.m_cellSize = { latticeSize-1, latticeSize-1};
    simulationMesh.m_gridDX = 0.25;
    simulationMesh.m_nFrames = 1;
    
    std::function<void(float [latticeSize][latticeSize])> f = [&] (float inputLattice[latticeSize][latticeSize]) {
        simulationMesh.regFrame(inputLattice);
    };
    
    delOperator.initialiseLattice();
    delOperator.registerWriteFunctionCallback(f);
    
    // Initialize the simulation example
    simulationMesh.initialize(delOperator.m_lattice);

    // Output the initial shape of the surface
    simulationMesh.writeFrame(0);

    delOperator.conjugateGradientsSolver();
    
    // Write the entire timeline to USD
    simulationMesh.writeUSD();

    return 0;
}


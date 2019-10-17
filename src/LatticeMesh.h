#include <iostream>

#include "AnimatedMesh.h"

template<class T, int LatticeSize>
struct LatticeMesh : public AnimatedMesh<T, 4>
{
    using Base = AnimatedMesh<T, 4>;
    using Base::m_meshElements;
    using Base::m_particleX;
    using Base::initializeUSD;
    using Base::initializeTopology;
    using Base::initializeParticles;
    using Base::writeFrame;
    
    std::array<int, 2> m_cellSize; // dimensions in grid cells
    float m_latticeValues[LatticeSize][LatticeSize];
    T m_gridDX;
    int m_nFrames;
    
    void setLatticeValues(float inputLattice[LatticeSize][LatticeSize]) {
        for(int i = 0; i < LatticeSize; i++) {
            for(int j = 0; j < LatticeSize; j++) {
                m_latticeValues[i][j] = inputLattice[i][j];
            }
        }
    }
    
    void initialize(float inputLattice[LatticeSize][LatticeSize])
    {
        initializeUSD("ConjugateGradient.usda");
        
        // Create a Cartesian lattice topology
        for(int cell_i = 0; cell_i < m_cellSize[0]; cell_i++) {
            for(int cell_j = 0; cell_j < m_cellSize[1]; cell_j++) {
                m_meshElements.emplace_back(
                    std::array<int, 4>{
                        gridToParticleID(cell_i  , cell_j  ),
                        gridToParticleID(cell_i+1, cell_j  ),
                        gridToParticleID(cell_i+1, cell_j+1),
                        gridToParticleID(cell_i  , cell_j+1)
                    }
                    );
            }
        }
        
        initializeTopology();
        
        //Set m_latticeValues
        setLatticeValues(inputLattice);
        
        // Also initialize the associated particles
        for(int node_i = 0; node_i <= m_cellSize[0]; node_i++) {
            for(int node_j = 0; node_j <= m_cellSize[1]; node_j++) {
                m_particleX.emplace_back( (T)node_i*m_gridDX, (T)node_j*m_gridDX, (T)m_latticeValues[node_i][node_j]);
            }
        }
        initializeParticles();
    }
    
    void setFrame() {
        for(int node_i = 1; node_i <= m_cellSize[0]; node_i++) {
            for(int node_j = 1; node_j <= m_cellSize[1]; node_j++) {
                m_particleX[gridToParticleID(node_i  ,node_j)] = GfVec3f((T)node_i*m_gridDX, (T)node_j*m_gridDX, (T)m_latticeValues[node_i][node_j]);
                std::cout << "(" << node_i << "," << node_j << "): " << (T)m_latticeValues[node_i][node_j] << std::endl;
            }
        }
    }
    
    void regFrame(float inputLattice[LatticeSize][LatticeSize]) {
        static int frameNo = 1;
        setLatticeValues(inputLattice);
        setFrame();
        writeFrame(frameNo);
        frameNo += 1;
    }
    
private:
    inline int gridToParticleID(const int i, const int j) { return i * (m_cellSize[1]+1) + j; }
};


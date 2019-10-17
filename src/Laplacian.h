#include <iostream>
#include <vector>

#include "CGVectorWrapper.h"
#include "ConjugateGradientSolver.h"

template <int gridWidth>
struct Laplacian2D {
    float m_lattice[gridWidth][gridWidth];
    float m_h, m_hSquared;
    float m_b = 2.0;
    std::function<void(float [gridWidth][gridWidth]) > m_callback;
    
    void initialiseLattice() {
        m_h = gridWidth-1;
        m_hSquared = m_h*m_h;
        for(int i = 0; i < gridWidth; i++) {
            for(int j = 0; j < gridWidth; j++)
                m_lattice[i][j] = 0.0;
        }
    }
    
    void registerWriteFunctionCallback(std::function<void(float [gridWidth][gridWidth])> cb) {
        if(cb) {
            m_callback = cb;
        }
    }
    
    void captureXIteration(std::vector<float> x) {
        for(int i = 0; i < gridWidth; i++) {
            for(int j = 0; j < gridWidth; j++) {
                m_lattice[i][j] = x[gridToParticleID(i, j)];
            }
        }
        if (m_callback) {
            m_callback(m_lattice);
        }
    }
    
    void applyStencilOperation(const std::vector<float> &x, std::vector<float> &y) {
        float centralWeight = 4/m_hSquared;
        float edgeWeights = 1.0/m_hSquared;
        for(int i = 1; i < gridWidth-1; i++) {
            for(int j = 1; j < gridWidth-1; j++) {

                int pCenter = gridToParticleID(i,j);
                int pPlusX = gridToParticleID(i+1,j);
                int pMinusX = gridToParticleID(i-1,j);
                int pPlusY = gridToParticleID(i,j+1);
                int pMinusY = gridToParticleID(i,j-1);
                
                y[pCenter] =  (centralWeight * x[pCenter]) -
                (edgeWeights * (x[pPlusX] + x[pMinusX] + x[pPlusY] + x[pMinusY]));
            }
        }
        setBoundaryParticles(y);
    }
    
    void conjugateGradientsSolver() {
        std::function<void(const std::vector<float> &x, std::vector<float> &y)> lambdaMultiply = [&] (const std::vector<float> &x, std::vector<float> &y) {
            applyStencilOperation(x, y);
        };
        
        std::function<void(std::vector<float> &x)> lambdaSetBoundary = [&] ( std::vector<float> &x) {
            setBoundaryParticlesToZero(x);
        };

        std::vector<float> rhs(gridWidth*gridWidth, 0.0);
        std::vector<float> q(gridWidth*gridWidth, 0.0);
        std::vector<float> x(gridWidth*gridWidth, 0.0);
        std::vector<float> s(gridWidth*gridWidth, 0.0);
        std::vector<float> r(gridWidth*gridWidth, 0.0);
        
        CGVectorWrapper<float> rhsWrapper(rhs);
        CGVectorWrapper<float> qWrapper(q);
        CGVectorWrapper<float> xWrapper(x);
        CGVectorWrapper<float> sWrapper(s);
        CGVectorWrapper<float> rWrapper(r);
        
        qWrapper.registerMultiplyCb(lambdaMultiply);
        qWrapper.registerProjectCb(lambdaSetBoundary);
        rWrapper.registerProjectCb(lambdaSetBoundary);
        
        //initialise RHS to b/(h*h)
        for(int i = 0; i < rhs.size(); i++ ) {
            rhs[i] = -m_b/m_hSquared;
        }
        
        for(int i = 0; i < gridWidth; i++) {
            for(int j = 0; j < gridWidth; j++)
                x[gridToParticleID(i,j)] = 0.0;
        }
        setBoundaryParticles(x);
        
        std::function<void(std::vector<float>)> cbFunction = [=](std::vector<float> xIter) {
            this->captureXIteration(xIter);
        };
        
        ConjugateGradientSolver<float>::Solve(xWrapper, rhsWrapper, qWrapper, sWrapper, rWrapper, 1e-4, 50, cbFunction);
    }

    void setBoundaryParticlesToZero(std::vector<float> &v) {
        for(int i = 0; i < gridWidth; i++) {
            v[gridToParticleID(0,i)] = 0.0;
            v[gridToParticleID(gridWidth-1,i)] = 0.0;
            v[gridToParticleID(i, 0)] = 0.0;
            v[gridToParticleID(i, gridWidth-1)] = 0.0;
        }
    }
    
    void setBoundaryParticles(std::vector<float> &v) {
        for(int i = 0; i < gridWidth; i++) {
            v[gridToParticleID(0,i)] = -m_b/m_hSquared;
            v[gridToParticleID(gridWidth-1,i)] = -m_b/m_hSquared;
            v[gridToParticleID(i, 0)] = -m_b/m_hSquared;
            v[gridToParticleID(i, gridWidth-1)] = -m_b/m_hSquared;
        }
    }
    
private:
    int gridToParticleID(const int x, const int y) {
        return x*(gridWidth)+y;
    }
    
};

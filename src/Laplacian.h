#include <iostream>
#include <vector>
#include <memory>

#include "CGVectorWrapper.h"
#include "ConjugateGradientSolver.h"
#include "JacobiSolver.h"

template <class T ,int gridWidth>
struct Laplacian2D {
    using LatticeType = T(&)[gridWidth][gridWidth];
    using VectorType = T (&)[gridWidth*gridWidth];
    
    //aligned_alloc
    void *aMemory = NULL;
    std::size_t s = 8;
    std::size_t cap = sizeof(T) * gridWidth * gridWidth;
    int x = posix_memalign(&aMemory , s, cap);
    T **A = reinterpret_cast<T **>(aMemory);
    
    //Member Variables
    LatticeType m_lattice = reinterpret_cast<LatticeType>(*A);
    uint8_t **m_nodeType;
    float m_h, m_hSquared;
    float m_b = 2.0;
    std::function<void (T [][gridWidth])> m_callback;
    
    ~Laplacian2D() {
        free(m_lattice);
    }
    
    void initialiseLattice(std::array<std::array<T, gridWidth>, gridWidth> lattice, std::array<std::array<uint8_t, gridWidth>, gridWidth> nodetype) {
//        m_lattice = new T*[gridWidth];
        m_nodeType = new uint8_t*[gridWidth];
        for(int i = 0; i < gridWidth; i++) {
//            m_lattice[i] = new T[gridWidth];
            m_nodeType[i] = new uint8_t[gridWidth];
        }
        
        m_h = gridWidth-1;
        m_hSquared = m_h*m_h;
        for(int i = 0; i < gridWidth; i++) {
            for(int j = 0; j < gridWidth; j++) {
                //m_lattice[i][j] = 0.0;
                //std::cout << (*(m_lattice + i) + j) << std::endl;
                m_lattice[i][j] = (T)lattice[i][j];
                m_nodeType[i][j] = nodetype[i][j];
            }
        }
    }
    
    void registerWriteFunctionCallback(std::function<void(T [][gridWidth])> cb) {
        if(cb) {
            m_callback = cb;
        }
    }

    void captureXIteration(VectorType x) {
        for(int i = 0; i < gridWidth; i++) {
            for(int j = 0; j < gridWidth; j++) {
                m_lattice[i][j] = x[gridToParticleID(i, j)];
            }
        }
        if (m_callback) {
            m_callback(m_lattice);
        }
    }

    T getLatticeValueAtPosition(int x, int y, VectorType &lat) {
//        if (m_nodeType[x][y] == 2 || m_nodeType[x][y] == 1) {
//            return (T)0;
//        }
//        else {
            return lat[gridToParticleID(x,y)];
//        }
    }

    bool isActive(int x, int y) {
        if (m_nodeType[x][y] == 0) {
            return true;
        }
        return false;
    }

    void applyStencilOperation(VectorType &x, VectorType &y) {
        float centralWeight = 4/m_hSquared;
        float edgeWeights = 1.0/m_hSquared;
        for(int i = 1; i < gridWidth-1; i++) {
            for(int j = 1; j < gridWidth-1; j++) {
                
                if (isActive(i, j) == false) {
                    //std::cout << "Skipping : " << i << " " << j << std::endl;
                    continue;
                }
                
                T valCenter = getLatticeValueAtPosition(i, j, x);
                T valPlusX = getLatticeValueAtPosition(i+1,j, x);
                T valMinusX = getLatticeValueAtPosition(i-1,j, x);
                T valPlusY = getLatticeValueAtPosition(i,j+1, x);
                T valMinusY = getLatticeValueAtPosition(i,j-1, x);

                y[gridToParticleID(i,j)] =  (centralWeight * valCenter) -
                (edgeWeights * (valPlusX + valMinusX + valPlusY + valMinusY));
            }
        }
        setBoundaryParticles(y);
    }

    void checkConvergence(VectorType x) {
        std::vector<float> err;
        std::cout << "CONVERGED SOLUTION. Value close to: " << -m_b/m_hSquared << std::endl;
        float centralWeight = 4/m_hSquared;
        float edgeWeights = 1.0/m_hSquared;
        for(int i = 1; i < gridWidth-1; i++) {
            for(int j = 1; j < gridWidth-1; j++) {

                int pCenter = gridToParticleID(i,j);
                int pPlusX = gridToParticleID(i+1,j);
                int pMinusX = gridToParticleID(i-1,j);
                int pPlusY = gridToParticleID(i,j+1);
                int pMinusY = gridToParticleID(i,j-1);

                auto y = (centralWeight * x[pCenter]) -
                (edgeWeights * (x[pPlusX] + x[pMinusX] + x[pPlusY] + x[pMinusY]));
                //std::cout << y << std::endl;
                err.push_back((-m_b/m_hSquared) - y);
            }
        }
        double errT = 0.0;
        for(auto it = err.begin(); it != err.end(); it++) {
            errT += (*it) * (*it);
        }
        std::cout << "MSE: " << std::sqrt(errT) << std::endl;
    }

    VectorType getAlignedMemory(int cap) {
        void *aMemory = NULL;
        std::size_t s = 8;
        cap = sizeof(T) * cap;
        if(posix_memalign(&aMemory , s, cap)) {
            throw std::logic_error("Failed to create aligned memory");
        }
        T **A = reinterpret_cast<T **>(aMemory);
        
        return reinterpret_cast<VectorType>(*A);
    }
    
    void conjugateGradientsSolver() {
        std::function<void(VectorType &x, VectorType &y)> lambdaMultiply = [&] (VectorType &x, VectorType &y) {
            applyStencilOperation(x, y);
        };

        std::function<void(VectorType &x)> lambdaSetBoundary = [&](VectorType &x) {
            setBoundaryParticlesToZero(x);
        };

        int len = gridWidth*gridWidth;
        
        VectorType rhs = getAlignedMemory(len);
        VectorType q = getAlignedMemory(len);
        VectorType x = getAlignedMemory(len);
        VectorType s = getAlignedMemory(len);
        VectorType r = getAlignedMemory(len);

        for(int i = 0; i < len; i++) {
            q[i] = 0;
            r[i] = 0;
            s[i] = 0;
            x[i] = 0;
        }

        CGVectorWrapper<float, gridWidth*gridWidth> rhsWrapper(rhs);
        CGVectorWrapper<float, gridWidth*gridWidth> qWrapper(q);
        CGVectorWrapper<float, gridWidth*gridWidth> xWrapper(x);
        CGVectorWrapper<float, gridWidth*gridWidth> sWrapper(s);
        CGVectorWrapper<float, gridWidth*gridWidth> rWrapper(r);

        qWrapper.registerMultiplyCb(lambdaMultiply);
        qWrapper.registerProjectCb(lambdaSetBoundary);
        rWrapper.registerProjectCb(lambdaSetBoundary);

        //initialise RHS to b/(h*h)
        for(int i = 0; i < gridWidth; i++ ) {
            for(int j = 0; j < gridWidth; j++) {
                if(m_nodeType[i][j] != 0)
                    rhs[gridToParticleID(i, j)] = -m_b/m_hSquared;
                else {
                    rhs[gridToParticleID(i, j)] = (T)0;
                }
            }
        }

        for(int i = 0; i < gridWidth; i++) {
            for(int j = 0; j < gridWidth; j++)
                x[gridToParticleID(i,j)] = 0.0;
        }
        setBoundaryParticles(x);

        std::function<void(CGVectorWrapper<float, gridWidth*gridWidth>)> cbFunction = [&](CGVectorWrapper<float, gridWidth*gridWidth> xIter) {
            this->captureXIteration(xIter.m_data);
        };

        ConjugateGradientSolver<float, gridWidth>::Solve(xWrapper, rhsWrapper, qWrapper, sWrapper, rWrapper, 1e-5, 100, cbFunction);

        checkConvergence(x);

        free(rhs);
        free(r);
        free(s);
        free(q);
        free(x);

    }

    void JacobiIterationSolver() {
        std::function<void(VectorType &x, VectorType &y)> lambdaMultiply = [&] (VectorType &x, VectorType &y) {
            applyStencilOperation(x, y);
        };

        std::function<void(VectorType &x)> lambdaSetBoundary = [&](VectorType &x) {
            setBoundaryParticlesToZero(x);
        };
        int len = gridWidth * gridWidth;
        VectorType rhs = getAlignedMemory(len);
        VectorType q = getAlignedMemory(len);
        VectorType x = getAlignedMemory(len);
        VectorType dInverse = getAlignedMemory(len);
        VectorType r = getAlignedMemory(len);

        for(int i = 0; i < len; i++) {
            q[i] = 0;
            r[i] = 0;
            dInverse[i] = m_hSquared/4.0;
            if(i == 0 || i == len-1) {
                dInverse[i] = m_hSquared;
            }
        }

        CGVectorWrapper<float, gridWidth*gridWidth> rhsWrapper(rhs);
        CGVectorWrapper<float, gridWidth*gridWidth> qWrapper(q);
        CGVectorWrapper<float, gridWidth*gridWidth> xWrapper(x);
        CGVectorWrapper<float, gridWidth*gridWidth> dWrapper(dInverse);
        CGVectorWrapper<float, gridWidth*gridWidth> rWrapper(r);

        qWrapper.registerMultiplyCb(lambdaMultiply);
        qWrapper.registerProjectCb(lambdaSetBoundary);
        rWrapper.registerProjectCb(lambdaSetBoundary);

        //initialise RHS to b/(h*h)
        for(int i = 0; i < gridWidth; i++ ) {
            for(int j = 0; j < gridWidth; j++) {
                //if(m_nodeType[i][j] != 0)
                rhs[gridToParticleID(i, j)] = -m_b/m_hSquared;
                //else {
                //    rhs[gridToParticleID(i, j)] = (T)0;
                //}
            }
        }

        for(int i = 0; i < gridWidth; i++) {
            for(int j = 0; j < gridWidth; j++)
                x[gridToParticleID(i,j)] = 0.0;
        }
        setBoundaryParticles(x);

        std::function<void(CGVectorWrapper<float, gridWidth*gridWidth>)> cbFunction = [&](CGVectorWrapper<float, gridWidth*gridWidth> xIter) {
            this->captureXIteration(xIter.m_data);
        };

        JacobiSolver<float>::Solve(xWrapper, rhsWrapper, qWrapper, dWrapper, rWrapper, 1e-5, 500, cbFunction);
        setBoundaryParticles(x);
        checkConvergence(x);

        free(rhs);
        free(q);
        free(x);
        free(dInverse);
        free(r);
        
    }

    void setBoundaryParticlesToZero(VectorType &v) {
        for(int i = 0; i < gridWidth; i++) {
            for(int j = 0; j < gridWidth; j++) {
                if(m_nodeType[i][j] == 2 || m_nodeType[i][j] == 1) {
                    v[gridToParticleID(i, j)] = (T)0;
                }
            }
        }
    }

    void setBoundaryParticles(VectorType &v) {
        for(int i = 0; i < gridWidth; i++) {
            for(int j = 0; j < gridWidth; j++) {
                if(m_nodeType[i][j] == 2 || m_nodeType[i][j] == 1) {
                    v[gridToParticleID(i,j)] = -m_b/m_hSquared;
                }
            }
        }
    }

private:
    int gridToParticleID(const int x, const int y) {
        return x*(gridWidth)+y;
    }
};

#include <torch/torch.h>
//#include <torch/csrc/autograd/generated/variable_factories.h> //for from_blob
#include <ATen/Tensor.h>
#include <torch/csrc/utils/tensor_new.h>
#include <iostream>

template <class T, int gridWidth>
class LaplacianWithTensor {
public:
    torch::Tensor tLattice;
    torch::Tensor tnodeType;
    float m_h, m_hSquared;
    float m_b;
    
    LaplacianWithTensor(std::array<std::array<T, gridWidth>, gridWidth> lattice,
                        std::array<std::array<uint8_t, gridWidth>, gridWidth> nodeType) {
//        tLattice = torch::from_blob(&lattice, {gridWidth, gridWidth});
        auto options = torch::TensorOptions().dtype(torch::kFloat32);
        tLattice = torch::zeros({gridWidth, gridWidth}, options);
        tnodeType = torch::zeros({gridWidth, gridWidth}, options);
//        tLattice.new_tensor(lattice);
//        tnodeType.new_tensor(nodeType);
        m_b = 2.0;
        m_h = gridWidth-1;
        m_hSquared = m_h * m_h;
    }
    
    void JacobiIterationSolver() {
        auto options = torch::TensorOptions().dtype(torch::kFloat32);
        int len = gridWidth * gridWidth;
        torch::Tensor rhs = torch::full(len, -m_b/m_hSquared, options);
        torch::Tensor q = torch::zeros(len, options);
        torch::Tensor x = torch::zeros(len, options);
        torch::Tensor dInverse = torch::zeros(len, options);
        torch::Tensor r = torch::zeros(len, options);
        
        
    }
    
};

#include <iostream>
#include <functional>

template<class T, int gridWidth>
struct ConjugateGradientSolver
{
    template<class LatticeType>
    static bool Solve(LatticeType& x, const LatticeType& b,
                      LatticeType& q, LatticeType& p, LatticeType& r,
                      const T tolerance,const int max_iterations,
                      std::function<void(LatticeType)> cb=NULL)
    {
        T rhoOld = std::numeric_limits<T>::max();
        T convergence_norm = 0;
        r = b;
        q.MultiplyWithStencil(x);
        r -= q;
        r.ProjectToZero();
        for(int iterations=0; iterations < max_iterations ; iterations++){
            convergence_norm = r.ConvergenceNorm();
            std::cout << "Residual norm = " << convergence_norm << " [after " << iterations << " CG iterations]" << std::endl;
            if(convergence_norm <= tolerance)
                return true;
            if(iterations == max_iterations)
                break;
            T rho = r.InnerProduct(r);
            if(iterations == 0)
                p = r;
            else
                p.Saxpy(rho / rhoOld, p, r);
            q.MultiplyWithStencil(p);
            q.ProjectToZero();
            T p_dot_q = q.InnerProduct(p);
            if(p_dot_q <= 0)
                std::cout << "CG: matrix appears indefinite or singular, p_dot_q/p_dot_p="
                << p_dot_q / p.InnerProduct(p) << std::endl;
            T alpha = rho / p_dot_q;
            x.Saxpy(alpha, p, x);
            r.Saxpy(-alpha, q, r);
            rhoOld = rho;
            
            std::cout << "X norm : " << x.ConvergenceNorm() << std::endl;
            if(cb) {
                cb(x);
            }
        }
        std::cout << "cg not converged after " << max_iterations << " iterations, residual norm = " << convergence_norm << std::endl;
        return false;
    }
};

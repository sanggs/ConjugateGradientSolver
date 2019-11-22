#include <iostream>
#include <functional>

template<class T>
struct JacobiSolver
{
    template<class LatticeType>
    static bool Solve(LatticeType& x, const LatticeType& b,
                      LatticeType& q, LatticeType& dInverse, LatticeType& r,
                      const T tolerance,const int max_iterations,
                      std::function<void(LatticeType)> cb=NULL)
    {
        T convergence_norm = 0;
        r = b;
        q.MultiplyWithStencil(x);
        r -= q;
        r.ProjectToZero();
        for(int iterations=0; iterations < max_iterations ; iterations++) {
            //Check convergence
            convergence_norm = r.ConvergenceNorm();
            std::cout << "Residual norm = " << convergence_norm << " [after " << iterations << " CG iterations]" << std::endl;
            if(convergence_norm <= tolerance)
                return true;
            if(iterations == max_iterations)
                break;
            //Update x
            r.Multiply(dInverse);
            r.MultiplyWithScalar((T)2/3);
            r.ProjectToZero();
            x += r;
            r = b;
            q.MultiplyWithStencil(x);
            r -= q;
            r.ProjectToZero();
            //Callback to captureX
            std::cout << "X norm : " << x.ConvergenceNorm() << std::endl;
            if(cb) {
                cb(x);
            }
        }
        std::cout << "cg not converged after " << max_iterations << " iterations, residual norm = " << convergence_norm << std::endl;
        return false;
    }
};

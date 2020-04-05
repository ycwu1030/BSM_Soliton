#ifndef QUARTIC_SOLVER_H
#define QUARTIC_SOLVER_H

#include <complex>
#include <vector>

typedef std::complex<double> cdouble;

class QuarticSolver
{
private:
    cdouble a,b,c,d,e;
    cdouble P,Q,D,u,v;
    cdouble y;

    cdouble omg1;
    cdouble omg2;
    
public:
    QuarticSolver();
    QuarticSolver(double a4, double a3, double a2, double a1, double a0);
    ~QuarticSolver(){};

    void SetUpEquation(double a4=1, double a3=3, double a2=4, double a1=2, double a0=9);
    void Solve();
    void Solve(double a4, double a3, double a2, double a1, double a0);

    cdouble SOLUTIONS[4];
    std::vector<double> GetRealSolution();
};


#endif //QUARTIC_SOLVER_H
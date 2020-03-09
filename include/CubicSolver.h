#ifndef CUBICSOLVER_H
#define CUBICSOLVER_H

#ifndef Pi
#define Pi 3.14159265358979312
#endif

// Here I only list the situation for real solution, I don't care the complex solution
enum CUBICSTATES
{
    ONEREAL,
    THREEEQUALREAL,
    ONETWO,
    THREEREAL
};

class CubicSolver
{
public:
    CubicSolver();
    CubicSolver(double a3, double a2, double a1, double a0); // a3 x^3 + a2 x^2 + a1 x + a0 
    ~CubicSolver(){};

    void SetUpEquation(double a3=2, double a2=9, double a1=13, double a0=6);
    void Solve();
    void Solve(double a3, double a2, double a1, double a0);

    CUBICSTATES STATE;
    double SOLUTIONS[3];

private:
    double A0,A1,A2,A3;

    double b;// a2/(3 a3)
    double c;// a1/(6 a3)
    double d;// a0/(2 a3)

    double alpha; // - b^3 + 3bc - d
    double beta; // b^2 - 2c
    double Delta; // alpha^2 - beta^3
    double R1; // cbrt(alpha+sqrt(Delta))
    double R2; // cbrt(alpha-sqrt(Delta))
    double theta; // arccos(alpha/sqrt(beta^3))


};

#endif
#include <cmath>

#include "Basic_Model.h"
#include "DWSolver.h"

class ToyKink : public BSM_Soliton::BaseModel {
public:
    ToyKink() : v(1), lam(1), BaseModel(1){};
    ~ToyKink(){};

    virtual double V(VD field_values, double scale = 1) override {
        return lam * pow(pow(field_values[0], 2) - v * v, 2) / 4.0;
    }
    virtual VD dV(VD f, double scale = 1) override { return {lam * f[0] * (pow(f[0], 2) - v * v)}; }
    virtual VVD d2V(VD f, double scale = 1) override { return {{lam * (3 * pow(f[0], 2) - v * v)}}; }
    virtual double V_min(double scale = 1) override { return 0; }
    virtual VD Quartic_Couplings(VD f) override { return {lam}; }

protected:
    virtual void Calculate_Local_Extrema() override {
        Add_Local_Extremum({0});
        Add_Local_Extremum({v});
        Add_Local_Extremum({-v});
    }

private:
    double v;
    double lam;
};

int main(int argc, char const *argv[]) {
    ToyKink model;
    BSM_Soliton::DomainWallSolver DW(&model);
    DW.Set_Allocation_Parameters(0.0);
    DW.Solve({-1}, {1});
    DW.Dump_Solution("ToyKink_test.dat");
    return 0;
}

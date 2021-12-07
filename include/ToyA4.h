#ifndef TOY_A4_H
#define TOY_A4_H

#include "Basic_Model.h"

class ToyA4 : public Basic_Model {
 protected:
  double _mu2;
  double _g1;
  double _g2;

  double _v1global;
  double _v2global;
  double _v1;
  double _v2;
  // double _v3;

 public:
  ToyA4();
  ~ToyA4(){};

  void Set_Potential_Parameters(double mu2, double g1, double g2);

  virtual double Vtotal(VD field_values, double scale = 1);
  virtual VD dVtotal(VD field_values, double scale = 1);
  virtual VVD d2Vtotal(VD field_values, double scale = 1);
  virtual double V0_global(double scale = 1);
  virtual VD QuarticCoupling(VD field_values);

  virtual void FindLocalMinima();
  virtual void PrintParameters();

  double GetV1() { return _v1; }
  double GetV2() { return _v2; }
};

#endif  // TOY_A4_H

#ifndef RUNGEKUTTA_H_
#define RUNGEKUTTA_H_
#include"./../../Utilities/vecdef.h"

class rk4{
  private:
    flow3D k1, k2, k3, k4, fval_old, fval_new;
    double3D P;
    flow3D temp, temp2, var1, var2, var3;

  public:
    rk4();
    void solve();
    void calcvar3D(flow3D v_n, flow3D v_o, flow3D fl);
    virtual void gaslaw(flow3D f){};
    virtual void volflux(flow3D l, flow3D m);

};

#endif

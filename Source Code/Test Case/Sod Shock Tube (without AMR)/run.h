#include "main.decl.h"
#include "./Flux/Interface/Intflux.decl.h"
#include "Utilities/vecdef.h"
#include <cmath>

/*readonly*/ CProxy_Cell cellProxy;
/*readonly*/ CProxy_Intflux interfaceProxy;
/*readonly*/ int global_t_steps;

class Main: public CBase_Main{
  private:
    int iter, max_iter;

  public:
    Main(CkArgMsg*);
    void done();
};

class Cell: public CBase_Cell{
  protected:
    flow3D flux_c;
    flow3D val;
    double3D P;
    int numdiv;
    double gam;
    int local_tsteps;
    double local_dt, local_L;
    flow3D k1, k2, k3, k4, fval_old, fval_new;
    flow3D temp, temp2, var1, var2, var3;

  public:
    void gaslaw(flow3D);
    void volflux(flow3D,flow3D);
    void solve();

};

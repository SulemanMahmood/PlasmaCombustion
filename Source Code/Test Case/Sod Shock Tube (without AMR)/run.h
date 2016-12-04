#include "main.decl.h"
#include "vecdef.h"
#include <cmath>

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CProxy_Cell cellProxy;
/*readonly*/ CProxy_Flux fluxProxy;
/*readonly*/ CProxy_Intflux interfaceProxy;
/*readonly*/ int t_steps;

class Main: public CBase_Main{

  public:
    Main(CkArgMsg*);
    void done();
};

class Msg: public CMessage_Msg{
	
}

class Cell: public CBase_Cell{
  protected:
    flow3D flux_c;
    flow3D val;
    double3D P;
    int numdiv, iter;
    double gam;
    int t_steps;
    double dt, local_L;
    flow3D k1, k2, k3, k4, fval_old, fval_new;
    flow3D temp, temp2, var1, var2, var3;

  public:
    void gaslaw(flow3D);
    void solve();

};

class Flux: public CBase_Flux{
	private:
		int numdiv, iter;
};

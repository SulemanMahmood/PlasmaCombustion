#include "vecdef.h"
#include<string.h>
#include<stdlib.h>
#include <cmath>
#include "pup_stl.h"
#include "main.decl.h"

#define local_dt	0.23	// Value must be defined properly
#define N_DIM 123   // Value must go here

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CProxy_Cell cellProxy;
/*readonly*/ CProxy_Flux fluxProxy;
/*readonly*/ CProxy_intFlux interfaceProxy;
/*readonly*/ int t_steps;

class Main: public CBase_Main{

  public:
    Main(CkArgMsg*);
    void done();
};

/*class Msg: public CMessage_Msg{
	//int dir_cp;
	flow *val;
}*/

class Cell: public CBase_Cell{
	Cell_SDAG_CODE
  protected:
    flow3D flux_c;
    flow3D val;
    double3D P;
    int numdiv, iter;
    double gam;
    int t_steps;
    double dt, local_L;
    flow3D k1, k2, k3, k4, fval_old, fval_new, a;
    flow3D temp, temp2, var1, var2, var3;

  public:
		Cell();
    void gaslaw(flow3D);
    void solve();
		void initialize();
		void calcvar3D(flow3D v_n, flow3D v_o, flow3D fl);
		void WriteOutput();

};

class Flux: public CBase_Flux{
	Flux_SDAG_CODE
	
	private:
		int numdiv, iter;
	public:
		Flux();
		void inviscidFlux(flow3D fl[], flow3D a);
		void fluxFacetoCell();
};

class intFlux: public CBase_intFlux{
	intFlux_SDAG_CODE
	private:

	public:
		intFlux();
		void wall();
		void inviscidFlux();
};

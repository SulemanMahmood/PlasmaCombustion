#ifndef RUN_H
#define RUN_H

#include <vector>
#include "utility.h"
#include "run.decl.h"

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CProxy_Cell cellProxy;
/*readonly*/ CProxy_Flux fluxProxy;
/*readonly*/ CProxy_Interface interfaceProxy;
/*readonly*/ int dimX;
/*readonly*/ int dimY;
/*readonly*/ int dimZ;
/*readonly*/ int t_steps;
/*readonly*/ double dt;
/*readonly*/ double dx;
/*readonly*/ int ndiv;
/*readonly*/ double gma;

class Main: public CBase_Main{
	public:
		Main(CkArgMsg*);
		void done();
};

class Cell: public CBase_Cell{
	private:
		Cell_SDAG_CODE
		flow3D val_new, val_old;
		double3D P;
		int iter;
	public:
		Cell();
		Cell(CkMigrateMessage* m){}
		void update();
		void initialize();
		void gaslaw();
		void calcvar3D(flow3D,flow3D,flow3D);
};

class Flux: public CBase_Flux{
	private:
		Flux_SDAG_CODE
		flow3D flux_c, cell_val;
		flow4D flux_f;
		double3D P;
		int iter;
	public:
		Flux();
		Flux(CkMigrateMessage* m){}
		void inviscidFlux();
		void fluxFacetoCell();
};

class Interface: public CBase_Interface{
	private:
		Interface_SDAG_CODE
		flow2D val_l, val_r, flux;
		double2D P_left, P_right;
	public:
		Interface();
		Interface(CkMigrateMessage* m){}
		void calc();
		void wall(flow2D,flow2D,double2D,double2D);
};

#endif

#ifndef RUN_H
#define RUN_H

#include "pup_stl.h"
#include "utility.h"
#include "run.decl.h"
#include <stdlib.h>

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
		void InterfaceIsUp();
};

class Cell: public CBase_Cell{
	public:
		Cell_SDAG_CODE
		flow3D val_new, val_old;
		double3D P;
		int iter;

		Cell();
		Cell(CkMigrateMessage* m){}
		void pup(PUP::er &p){
			CBase_Cell::pup(p);
			__sdag_pup(p);
			p|val_new;
			p|val_old;
			p|P;
			p|iter;
		}
		void update();
		void initialize();
		void gaslaw();
		void calcvar3D(flow3D&,flow3D,flow3D);
};

class Flux: public CBase_Flux{
	public:
		Flux_SDAG_CODE
		flow3D flux_c, cell_val;
		flow4D flux_f;
		double3D P;
		//int iter;

		Flux();
		Flux(CkMigrateMessage* m){}
		void pup(PUP::er &p){
			CBase_Flux::pup(p);
			__sdag_pup(p);
			p|flux_c;
			p|flux_f;
			p|P;
			p|cell_val;
		}
		void inviscidFlux();
		void fluxFacetoCell();
};

class Interface: public CBase_Interface{
	public:
		Interface_SDAG_CODE
		flow2D val_l, val_r, flux;
		double2D P_left, P_right;
		int iter;

		Interface();
		Interface(CkMigrateMessage* m){}
		void pup(PUP::er &p){
			CBase_Interface::pup(p);
			__sdag_pup(p);
			p|val_l;
			p|val_r;
			p|P_left;
			p|P_right;
			p|flux;
		}
		void calc();
		void wall(flow2D&,flow2D,double2D&,double2D);
};

#endif

#ifndef RUN_H
#define RUN_H

#include "pup_stl.h"
#include "utility.h"
#include "run.decl.h"
#include <stdlib.h>

#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

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

// readonly variables for chemical reactions
/*readonly*/ double Te; //in K
/*readonly*/ double Tg; // in K
/*readonly*/ double end_time_chem; // in s
/*readonly*/ double dt_chem; // in s
/*readonly*/ int iter_chem;
/*readonly*/ double R; // Gas constat in J/mol.K
/*readonly*/ double Av; // Avogadro's number
/*readonly*/ double n; // Number density of air
/*readonly*/ double eq; // Equivalence Ratio

/*readonly*/ double1D H_f; // Enthalpy of rxn;
/*readonly*/ double2D K; // rxn rate constant for each rxn
/*readonly*/ int2D rs; // reactant species for each rxn
/*readonly*/ double1D adv; // advancement for each reaction
/*readonly*/ int2D p_rxn; // participating rxn for each species
/*readonly*/ double2D r_p; // reactants or products
/*readonly*/ double2D add_info; // Third body efficiencies and pressure dependence
/*readonly*/ int2D tb_sp; // Third body
/*readonly*/ int size; // Number of species
/*readonly*/ int rxn_size; // Species rxns
/*readonly*/ double1D sp; // species concentration
/*readonly*/ int wf; // Frequency of writing output to file
/*readonly*/ string1D species;
/*readonly*/ double1D Cp;

class Main: public CBase_Main{
	public:
		Main(CkArgMsg*);
		void done();
		void InterfaceIsUp();
        void read_file(); // for chemical reactions
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
    
        // functions for calculating chemical reactions
        void solve_rxn();
        void calc_change(double1D&, double1D&);
        void write_file(int);
        void calc_temp(double1D&);
        void initialize_chem();
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

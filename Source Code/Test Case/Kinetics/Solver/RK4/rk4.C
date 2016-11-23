#include <cmath>
#include "../../Utilities/vecdef.h"
#include "../../Utilities/vardef.h"
#include "rk4.h"

reaction::solve_rxn(double1D sp, double2D K, double T_e, double T_g, double end_time, double dt){
  reaction(sp.size(),T_e,T_g);
  iter = int(end_time/dt);
  for (int i = 0; i < iter; i++){
    sp1 = sp;
    calc_change(k1,sp1,K,dt);
    calc_val(sp2,sp1,k1,dt);
    calc_change(k2,sp2,K,dt);
    calc_val(sp3,sp2,k2,dt);
    calc_change(k3,sp3,K,dt);
    calc_val(sp4,sp3,k3,dt);
    calc_change(k4,sp4,K,dt);
    k = k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;
    calc_val(sp,sp1,k,dt);
  }
}

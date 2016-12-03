#include "flux.h"

void Flux::fluxFacetoCell(){
  for (int i = 0; i < numdiv; i++){
    for (int j = 0; j < numdiv; j++){
      for (int k = 0; k < numdiv; k++){
        flux_c[i][j][k] = flux_f[0][i+1][j][k] - flux_f[0][i][j][k] + flux_f[1][j][k][i] - flux_f[1][j+1][i][k] + flux_f[2][k][i][j] - flux_f[2][k+1][i][j];
      }
    }
  }
}

void Flux::volflux(flow3D a_n, flow3D a){

  inviscidFlux(flux_i,a);
  comm();
  //flux_f = 0.0 - flux_i;
  fluxFacetoCell();

}

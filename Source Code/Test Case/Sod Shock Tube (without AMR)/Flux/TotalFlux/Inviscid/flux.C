#include "flux.h"

void Flux::fluxFacetoCell(){
  for (int i = 0; i < numdiv; i++){
    for (int j = 0; j < numdiv; j++){
      for (int k = 0; k < numdiv; k++){
        flux_c[i][j][k] = flux_f[0][i][j][k] - flux_f[0][i+1][j][k] - flux_f[1][j+1][k][i] + flux_f[1][j][i][k] - flux_f[2][k+1][i][j] + flux_f[2][k][i][j];
      }
    }
  }
}

void Flux::volflux(){

  inviscidFlux();
  comm();
  fluxFacetoCell();

}

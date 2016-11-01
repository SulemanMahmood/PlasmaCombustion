#ifndef FLUX_H_
#define FLUX_H_

class totalFlux{
  public:

    void volflux(flow3D flux_f, flow3D val_c, doube3D P){

      flow3D flux_i[3];

      inviscidFlux(flux_i,val_c,P);
      comm(flux_i);
      fluxFacetoCell(flux_f,flux_i);
    }

    void fluxFacetoCell(flow3D cell, flow3D face[]){
      for (int i = 0; i < numdiv; i++){
        for (int j = 0; j < numdiv; j++){
          for (int k = 0; k < numdiv; k++){
            cell[i][j][k] = 0.0 - face[0][i+1][j][k] + face[0][i][j][k] - face[1][j+1][k][i] + face[1][j][i][k] - face[2][k+1][i][j] + face[2][k][i][j];
          }
        }
      }
    }

}

#endif

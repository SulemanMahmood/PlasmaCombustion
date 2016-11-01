#ifndef FLUX_H_
#define FLUX_H_

class totalFlux{
  public:

    void volflux(flow3D flux_f, flow3D val_c, doube3D P, int ndiv){

      flow3D flux_i[3], flux_v[3], flux_tb[3], flux_rd[3];

      inviscidFlux(flux_i,val_c,P,ndiv,gmix);
      viscousFlux(flux_v,val_c,P,ndiv);
      turbulentFlux(flux_tb,val_c,P,ndiv);
      reacDiffFlux(flux_rd,val_c,P,ndiv);
      vecselfsub4D(flux_v,flux_i,ndiv);
      vecselfadd4D(flux_tb,flux_rd,ndiv);
      vecselfadd4D(flux_v,flux_tb,ndiv);
      interfaceFlux();
      fluxFacetoCell(flux_f,flux_v,ndiv);
    }

    void fluxFacetoCell(flow3D cell, flow3D face[], int n_div){
      for (int i = 0; i < n_div; i++){
        for (int j = 0; j < n_div; j++){
          for (int k = 0; k < n_div; k++){
            cell[i][j][k] = face[0][i+1][j][k] - face[0][i][j][k] + face[1][j+1][k][i] - face[1][j][i][k] + face[2][k+1][i][j] - face[2][k][i][j];
          }
        }
      }
    }

}

#endif

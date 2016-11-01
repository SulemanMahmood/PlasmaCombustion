#ifndef IDEALGAS_H_
#define IDEALGAS_H_

class idealgas{
  private:
    double Cvmix;

  public:
    void gaslaw(double3D p, flow3D a, int ndiv){
      for (int x = 0; x < ndiv; x++){
        for (int y = 0; y < ndiv; y++){
          for (int z = 0; z < ndiv; z++){
            Cvmix = 0;
            for (int i = 0; i < a.Y.size(); i++){
              Cvmix += a[x][y][z].Y[i] * Cv[i];
            }
            p[x][y][z] = a[x][y][z].r*(a[x][y][z].E - 0.5*(a[x][y][z].u*a[x][y][z].u + a[x][y][z].v*a[x][y][z].v + a[x][y][z].w*a[x][y][z].w))/Cv_mix;
          }
        }
      }
    }
}

#endif

#ifndef IDEALGAS_H_
#define IDEALGAS_H_

class IDEALGAS{
  public:
    void gaslaw(double3D p, flow3D a, int ndiv){
      for (int x = 0; x < ndiv; x++){
        for (int y = 0; y < ndiv; y++){
          for (int z = 0; z < ndiv; z++){
            p[x][y][z] = (gam-1)*a[x][y][z].r*(a[x][y][z].E - 0.5*(a[x][y][z].u*a[x][y][z].u + a[x][y][z].v*a[x][y][z].v + a[x][y][z].w*a[x][y][z].w));
          }
        }
      }
    }
}

#endif

#ifndef RK4_H_
#define RK4_H_

class RK4{
  private:
    flow3D k1, k2, k3, k4, fval_old, fval_new;
    double3D P;
    flow3D temp, temp2, var1, var2, var3;

  public:

    RK4(){
      initializeflow3D(k1,numdiv);
      initializeflow3D(k2,numdiv);
      initializeflow3D(k3,numdiv);
      initializeflow3D(k4,numdiv);
      initializeflow3D(fval_old,numdiv);
      initializeflow3D(fval_new,numdiv);
      initializeflow3D(temp,numdiv);
      initializeflow3D(temp2,numdiv);
      initializeflow3D(var1,numdiv);
      initializeflow3D(var2,numdiv);
      initializeflow3D(var3,numdiv);
      initializedouble3D(P,numdiv);
    }

    void solve(){
      for (int i = 0; i < local_tsteps; i++){
        volflux(k1,fval_old,P);
        calcvar3D(var1,fval_old,k1);
        vecdiv3D(temp2,var1,double(2),numdiv);
        vecadd3D(temp,fval_old,temp2,numdiv);
        gaslaw(P,temp,numdiv);

        volflux(k2,temp,P,numdiv);
        calcvar3D(var2,var1,k2);
        vecdiv3D(temp2,var2,double(2),numdiv);
        vecadd3D(temp,fval_old,temp2,numdiv);
        gaslaw(P,temp,numdiv);

        volflux(k3,temp,P,numdiv);
        calcvar3D(var3,var2,k3);
        vecadd3D(temp,fval_old,var3,numdiv);
        gaslaw(P,temp,numdiv);

        volflux(k4,temp,P,numdiv);

        vecselfmul3D(k2,double(2),numdiv);
        vecselfmul3D(k3,double(2),numdiv);
        vecselfadd3D(k1,k2,numdiv);
        vecselfadd3D(k3,k4,numdiv);
        vecselfadd3D(k1,k3,numdiv);
        vecselfdiv3D(k1,double(6),numdiv);
        calcvar3D(fval_new,fval_old,k1);

        gaslaw(P,fval_new,numdiv);
      }
    }

    void calcvar3D(flow3D v_n, flow3D v_o, flow3D fl){
      double temp;
      for (int i = 0; i < numdiv; i++){
        for(int j = 0; j < numdiv; j++){
          for(int k = 0; k < numdiv; k++){
            v_n[i][j][k] = v_o[i][j][k] + local_dt *  fl[i][j][k] / local_L;
          }
        }
      }
    }

};

#endif

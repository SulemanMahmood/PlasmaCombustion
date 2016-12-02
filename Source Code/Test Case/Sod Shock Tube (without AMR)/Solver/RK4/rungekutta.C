rk4::rk4(){
  initialize(fval_old,P);
}

void rk4::solve(){
  for (int i = 0; i < local_tsteps; i++){
    volflux(k1,fval_old);
    calcvar3D(var1,fval_old,k1);
    vecdiv3D(temp2,var1,double(2),numdiv);
    vecadd3D(temp,fval_old,temp2,numdiv);
    gaslaw(temp);

    volflux(k2,temp);
    calcvar3D(var2,var1,k2);
    vecdiv3D(temp2,var2,double(2),numdiv);
    vecadd3D(temp,fval_old,temp2,numdiv);
    gaslaw(temp);

    volflux(k3,temp);
    calcvar3D(var3,var2,k3);
    vecadd3D(temp,fval_old,var3,numdiv);
    gaslaw(temp);

    volflux(k4,temp);

    vecselfmul3D(k2,double(2),numdiv);
    vecselfmul3D(k3,double(2),numdiv);
    vecselfadd3D(k1,k2,numdiv);
    vecselfadd3D(k3,k4,numdiv);
    vecselfadd3D(k1,k3,numdiv);
    vecselfdiv3D(k1,double(6),numdiv);
    calcvar3D(fval_new,fval_old,k1);

    gaslaw(fval_new);
  }
}


void rk4::calcvar3D(flow3D v_n, flow3D v_o, flow3D fl){
  for (int i = 0; i < numdiv; i++){
    for(int j = 0; j < numdiv; j++){
      for(int k = 0; k < numdiv; k++){
        v_n[i][j][k] = v_o[i][j][k] + local_dt *  fl[i][j][k] / local_L;
      }
    }
  }
}

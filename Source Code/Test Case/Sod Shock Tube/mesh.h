#ifndef MESH_H_
#define MESH_H_

class Mesh: public CBase_Mesh{
  friend class IDEALGAS, RK4, HRSLAU2;

  Mesh_SDAG_CODE

  private:
    int r_level, local_tsteps, numdiv;

  public:
    Mesh(){
      r_level = 1;
      numdiv = pow(2,r_level)*min_div;
      local_tsteps = ref_iter*pow(2,r_level);
    }

    Mesh(CkMigrateMessage *);

    void changemesh(int newlevel){
      flow3D temp_fval = fval_old;
      double3D temp_P = P;
      int temp_level = r_level;
      r_level = newlevel;
      numdiv = pow(2,r_level)*min_div;
      local_tsteps = ref_iter*pow(2,r_level);
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
      if (r_level < temp_level)
        coarsenmesh(fval_old,temp_fval,P,temp_P);
      else
        refinemesh(fval_old,temp_fval,P,temp_P);
    }

    void refinemesh(flow3D v_new, flow3D v_old, double3D P_new, double3D P_old){
      int n_size = v_new.size();
      int old_size = v_old.size();
      int s = int(n_size / old_size);
      for (int i = 0; i < n_size; i++){
        for (int j = 0; j < n_size; j++){
          for (int k = 0; k < n_size; k++){
            i_o = int(i/s);
            j_o = int(j/s);
            k_o = int(k/s);
            v_new[i][j][k] = v_old[i_o][j_o][k_o];
            P_new[i][j][k] = P_old[i_o][j_o][k_o];
          }
        }
      }
    }

    void coarsenmesh(flow3D v_new, flow3D v_old, double3D P_new, double3D P_old){
      
    }

}

#endif

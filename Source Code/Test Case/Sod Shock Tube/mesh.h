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
      numdiv = pow(2,r_level)*mindiv;
      local_tsteps = ref_iter*pow(2,r_level);
    }

    Mesh(CkMigrateMessage *);

    void refinemesh(){

    }

    void coarsenmesh(){
      
    }

}

#endif

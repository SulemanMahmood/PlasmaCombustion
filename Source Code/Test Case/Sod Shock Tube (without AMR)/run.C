#include "run.h"

Main::Main(CkArgMsg* msg){
  delete msg;
  iter = 0;
  cellProxy = CProxy_Cell::ckNew(dimX,dimY,dimZ);
  interfaceProxyface = CProxy_Intflux::ckNew(3,dimX+1,dimY,dimZ);
  cellProxy.run();
}

Main::void done(){
  if ((++iter) == max_iter){
    CkExit();
  }
}

void Cell::calcvar3D(flow3D v_n, flow3D v_o, flow3D fl){
  for (int i = 0; i < numdiv; i++){
    for(int j = 0; j < numdiv; j++){
      for(int k = 0; k < numdiv; k++){
        v_n[i][j][k] = v_o[i][j][k] + local_dt *  fl[i][j][k] / local_L;
      }
    }
  }
}

#include"main.def.h"

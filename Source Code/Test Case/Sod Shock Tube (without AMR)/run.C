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

#include"main.def.h"

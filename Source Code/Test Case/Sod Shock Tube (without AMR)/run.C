#include "main.decl.h"
#include "flux.h"
#include "idealgas.h"
#include "input.h"
#include "interface.h"
#include "mesh.h"
#include "rk4.h"

\*readonly*\ CProxy_Mesh cellProxy;
\*readonly*\ CProxy_Intflux interfaceProxy;

class Main: public CBase_Main{
  private:
    int iter, max_iter;

  public:
    Main(CkArgMsg* msg){
      delete msg;
      iter = 0;
      cellProxy = CProxy_Mesh::ckNew(dimX,dimY,dimZ);
      interfaceProxyface = CProxy_Intflux::ckNew(3,dimX+1,dimY,dimZ);
      cellProxy.run();
    }

    void done(){
      if ((++iter) == max_iter){
        CkExit();
      }
    }
}

#include"main.def.h"

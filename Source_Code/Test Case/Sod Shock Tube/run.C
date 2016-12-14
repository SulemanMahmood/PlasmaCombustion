#include "main.decl.h"
#include "flux.h"
#include "idealgas.h"
#include "input.h"
#include "interface.h"
#include "mesh.h"
#include "rk4.h"

\*readonly*\ CProxy_Mesh cell;
\*readonly*\ CProxy_Intflux face;

class Main: public CBase_Main{
  private:
    int iter;

  public:
    Main(CkArgMsg* msg){
      delete msg;
      iter = 0;
      cell = CProxy_Mesh::ckNew(dimX,dimY,dimZ);
      face = CProxy_Intflux::ckNew(3,dimX+1,dimY,dimZ);
      cell.run();
    }

    void done(){
      if ((++iter) == max_iter){
        CkExit();
      }
    }
}

#include"main.def.h"

#include "main.decl.h"
#include "./Flux/Interface/Intflux.decl.h"
#include "Utilities/vecdef.h"

/*readonly*/ CProxy_Cell cellProxy;
/*readonly*/// CProxy_Intflux interfaceProxy;

class Main: public CBase_Main{
  private:
    int iter, max_iter;

  public:
    Main(CkArgMsg*);
    void done();
};

class Cell: public CBase_Cell{
  protected:
    flow3D flux_c;
    double3D P;
    int numdiv;
    double gam;

  public:
    virtual void gaslaw(flow3D)=0;
    virtual void volflux(flow3D,flow3D)=0;

};

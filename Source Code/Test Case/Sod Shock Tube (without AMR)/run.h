#include "main.decl.h"

\*readonly*\ CProxy_Cell cellProxy;
\*readonly*\ CProxy_Intflux interfaceProxy;

class Main: public CBase_Main{
  private:
    int iter, max_iter;

  public:
    Main(CkArgMsg*);

    void done();
}

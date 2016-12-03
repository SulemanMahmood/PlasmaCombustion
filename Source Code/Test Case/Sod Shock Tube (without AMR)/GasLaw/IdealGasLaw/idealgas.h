#ifndef IDEALGAS_H_
#define IDEALGAS_H_

#include "../../run.h"

class IdealGas : public Cell{
  protected:
   double gam;
  public:
    void gaslaw(flow3D);
}

#endif

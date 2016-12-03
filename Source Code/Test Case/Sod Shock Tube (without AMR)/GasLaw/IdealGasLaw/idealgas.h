#ifndef IDEALGAS_H_
#define IDEALGAS_H_

#include "../../run.h"

class IdealGas : public Cell{
  public:
    void gaslaw(flow3D);
};

#endif

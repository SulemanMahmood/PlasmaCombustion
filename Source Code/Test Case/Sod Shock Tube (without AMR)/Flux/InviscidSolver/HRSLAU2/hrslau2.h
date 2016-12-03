#ifndef HRSLAU2_H_
#define HRSLAU2_H_

#include "../../TotalFlux/Inviscid/flux.h"
#include "../../../Utilities/vecdef.h"
#include <cmath>

class HRSLAU2: public Flux{
  public:
    void inviscidFlux(flow3D);
};

#endif

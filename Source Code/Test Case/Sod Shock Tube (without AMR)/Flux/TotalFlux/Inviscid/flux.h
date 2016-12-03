#ifndef FLUX_H_
#define FLUX_H_

#include "../../../run.h"

class Flux : public Cell{
  protected:
    flow3D flux_f[3];

  public:
    void fluxFacetoCell();
    void volflux();
    virtual void comm();
    virtual void inviscidFlux(flow3D);
};

#endif

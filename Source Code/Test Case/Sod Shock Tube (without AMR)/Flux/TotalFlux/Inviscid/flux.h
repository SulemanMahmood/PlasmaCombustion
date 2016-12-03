#ifndef FLUX_H_
#define FLUX_H_

#include "../../../run.h"

class Flux : public Cell{

  protected:
    flow3D flux_f[3], flux_i[3];

  public:
    void fluxFacetoCell();
    void volflux(flow3D,flow3D);
    virtual void comm();
    virtual void inviscidFlux(flow3D[],flow3D);
};

#endif

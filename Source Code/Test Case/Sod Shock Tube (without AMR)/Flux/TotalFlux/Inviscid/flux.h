#ifndef FLUX_H_
#define FLUX_H_

#include "../../../run.h"

class Flux{

  protected:
    flow3D flux_f[3], flux_i[3];

  public:
    void fluxFacetoCell();
    void volflux(flow3D,flow3D);
    void comm();
    void inviscidFlux(flow3D[],flow3D);
};

#endif

#ifndef FLUX_H_
#define FLUX_H_

class Flux : public Cell{
  protected:
    flow3D flux_i[3];

  public:
    void fluxFacetoCell();
    void volflux();
    virtual void comm(){}
    virtual void inviscidFlux(){}
}

#endif

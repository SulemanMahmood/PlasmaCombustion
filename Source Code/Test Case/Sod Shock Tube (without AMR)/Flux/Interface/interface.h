#ifndef INTERFACE_H_
#define INTERFACE_H_

class Interface : public CBase_Interface{

  private:
    double *f1, *f2;

  public:
    void inter_flux(Msg *m);

}


#endif

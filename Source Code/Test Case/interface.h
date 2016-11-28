#ifndef INTERFACE_H_
#define INTERFACE_H_


#include "../../Utilities/structdef.h"
#include "../../Utilities/vecdef.h"

class Interface : public CBase_Interface{

  private:
    double *f;

  public:

    void inter_flux(Msg *m){

      
      if (thisIndex.w == 0 && thisIndex.x == 0){
        inlet();
      }
      else if (thisIndex.w == 0 && thisIndex.x == dimX){
        outlet();
      }
      else if (thisIndex.w == 1 && (thisIndex.x == 0 || thisIndex.x == dimY)){
        wall();
      }
      else if (thisIndex.w == 2 && (thisIndex.x == 0 || thisIndex.x == dimZ)){
        wall();
      }
      else{
        calc();
      }

    }

}


#endif

#include "interface.h"

void Interface::inter_flux(Msg *m){

  int w = thisIndex.w;
  int x = thisIndex.x;
  int y = thisIndex.y;
  int z = thisIndex.z;

  if (w == 0){
    if (x == 0){
      wall();
      cellProxy(x,y,z).comm(m);
    }
    else if (x == dimX){
      wall();
      cellProxy(x-1,y,z).comm(m);
    }
    else{
      delete m;
      calc();
      cellProxy(x,y,z).comm(m1);
      cellProxy(x-1,y,z).comm(m2);
    }
  }

  else if (w == 1){
    if (x == 0){
      wall();
      cellProxy(z,x,y).comm(m);
    }
    else if (x == dimY){
      wall();
      cellProxy(z,x-1,y).comm(m);
    }
    else{
      delete m;
      calc();
      cellProxy(z,x,y).comm(m1);
      cellProxy(z,x-1,y).comm(m2);
    }
  }

  else{
    if (x == 0){
      wall();
      cellProxy(y,z,x).comm(m);
    }
    else if (x == dimZ){
      wall();
      cellProxy(y,z,x-1).comm(m);
    }
    else{
      delete m;
      calc();
      cellProxy(y,z,x).comm(m1);
      cellProxy(y,z,x-1).comm(m2);
    }
  }

}

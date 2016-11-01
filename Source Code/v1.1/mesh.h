#ifndef MESH_H_
#define MESH_H_

class Mesh: public CBase_Mesh{
  private:
    int r_level;

  public:
    Mesh(){
      r_level = 1;
    }
}

#endif

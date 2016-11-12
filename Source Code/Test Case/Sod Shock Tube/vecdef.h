#ifndef VECDEF_H_
#define VECDEF_H_

typedef vector<vector<vector<int>>> int3D;
typedef vector<vector<vector<double>>> double3D;
typedef vector<vector<vector<flow>>> flow3D;

void initializeflow3D(flow3D a, int nd){
  a.resize(nd);
  for (int i = 0; i < nd; i++){
    a[i].resize(nd);
    for (int j = 0; j < nd; j++){
      a[i][j].resize(nd);
    }
  }
}

void initializedouble3D(double3D a, int nd){
  a.resize(nd);
  for (int i = 0; i < nd; i++){
    a[i].resize(nd);
    for (int j = 0; j < nd; j++){
      a[i][j].resize(nd);
    }
  }
}


template<typename T>
typedef vector<vector<vector<T>>> T3D
void vecadd3D(T3D res, T3D lhs, T3D rhs, int dim){
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      for(int k = 0; k < dim; k++){
        res[i][j][k] = lhs[i][j][k] + rhs[i][j][k];
      }
    }
  }
}

template<typename T>
typedef vector<vector<vector<T>>> T3D
void vecdiv3D(T3D res, T3D lhs, double rhs, int dim){
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      for(int k = 0; k < dim; k++){
        res[i][j][k] = lhs[i][j][k] / rhs[i][j][k];
      }
    }
  }
}

template<typename T>
typedef vector<vector<vector<T>>> T3D
void vecselfadd3D(T3D lhs, T3D rhs, int dim){
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      for(int k = 0; k < dim; k++){
        lhs[i][j][k] += rhs[i][j][k];
      }
    }
  }
}

template<typename T>
typedef vector<vector<vector<T>>> T3D
void vecselfsub3D(T3D lhs, T3D rhs, int dim){
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      for(int k = 0; k < dim; k++){
        lhs[i][j][k] -= rhs[i][j][k];
      }
    }
  }
}

template<typename T>
typedef vector<vector<vector<T>>> T3D
void vecselfmul3D(T3D lhs, double rhs, int dim){
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      for(int k = 0; k < dim; k++){
        lhs[i][j][k] *= rhs;
      }
    }
  }
}

template<typename T>
typedef vector<vector<vector<T>>> T3D
void vecselfdiv3D(T3D lhs, double rhs, int dim){
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      for(int k = 0; k < dim; k++){
        lhs[i][j][k] /= rhs;
      }
    }
  }
}

template<typename T>
typedef vector<vector<vector<T>>> T3D
void vecselfadd4D(T3D lhs[], T3D rhs[], int dim){
  for(int n = 0; n < 3; n++){
    for(int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++){
        for(int k = 0; k < dim; k++){
          lhs[n][i][j][k] += rhs[n][i][j][k];
        }
      }
    }
  }
}

template<typename T>
typedef vector<vector<vector<T>>> T3D
void vecselfsub4D(T3D lhs[], T3D rhs[], int dim){
  for(int n = 0; n < 3; n++){
    for(int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++){
        for(int k = 0; k < dim; k++){
          lhs[n][i][j][k] -= rhs[n][i][j][k];
        }
      }
    }
  }
}

#endif

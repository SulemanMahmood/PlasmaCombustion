void IdealGas::gaslaw(flow3D a){
  for (int x = 0; x < numdiv; x++){
    for (int y = 0; y < numdiv; y++){
      for (int z = 0; z < numdiv; z++){
        P[x][y][z] = (gam-1)*a[x][y][z].r*(a[x][y][z].E - 0.5*(a[x][y][z].u*a[x][y][z].u + a[x][y][z].v*a[x][y][z].v + a[x][y][z].w*a[x][y][z].w));
      }
    }
  }
}

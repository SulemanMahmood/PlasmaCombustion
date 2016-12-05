#include "run.h"

Main::Main(CkArgMsg* msg){
  delete msg;
  mainProxy = thisProxy;
  cellProxy = CProxy_Cell::ckNew(dimX,dimY,dimZ);
  interfaceProxyface = CProxy_Intflux::ckNew(3,dimX+1,dimY,dimZ);
  fluxProxy = CProxy_Flux::ckNew(dimX,dimY,dimZ);
  cellProxy.run();
	fluxProxy.run_flux();
	interfaceProxy.run_interface();
}

Main::void done(){
  CkExit();
}

void Cell::calcvar3D(flow3D v_n, flow3D v_o, flow3D fl){
  for (int i = 0; i < numdiv; i++){
    for(int j = 0; j < numdiv; j++){
      for(int k = 0; k < numdiv; k++){
        v_n[i][j][k] = v_o[i][j][k] - local_dt *  fl[i][j][k] / local_L;
      }
    }
  }
}

void Cell::gaslaw(flow3D a){
  for (int x = 0; x < numdiv; x++){
    for (int y = 0; y < numdiv; y++){
      for (int z = 0; z < numdiv; z++){
        P[x][y][z] = (gam-1)*a[x][y][z].r*(a[x][y][z].E - 0.5*(a[x][y][z].u*a[x][y][z].u + a[x][y][z].v*a[x][y][z].v + a[x][y][z].w*a[x][y][z].w));
      }
    }
  }
}

void Flux::inviscidFlux(flow3D fl[], flow3D a){
  double r_l, r_r, r_h, P_l, P_r, g_mix_r, g_mix_l, C_p;
  flow a_l, a_r;
  double c_l, c_r, c_h, chi, g, g_HR, nx, ny, nz;
  double u_r, u_l, v_r, v_l, w_l, w_r, V, V_l, V_r;
  double M, M_L, M_R, f_L, f_R;
  double d_p, p, M_c, V_n, V_np, V_nn, m_dot;

  for(int n = 0; n < 3; n++){
    for (int i = 1; i < numdiv; i++){
      for (int j = 0; j < numdiv; j++){
        for (int k = 0; k < numdiv; k++){

          a_r = a[i][j][k];
          r_r = a_r.r;
          P_r = P[i][j][k];
          u_r = a_r.u/r_r;
          v_r = a_r.v/r_r;
          w_r = a_r.w/r_r;
          if (n == 0){
            a_l = a[i-1][j][k];
            r_l = a_l.r;
            P_l = P[i-1][j][k];
            u_l = a_l.u/r_l;
            v_l = a_l.v/r_l;
            w_l = a_l.w/r_l;
            nx = 1.0;
            ny = 0.0;
            nz = 0.0;
          }
          else if (n == 1){
            a_l = a[k][i-1][j];
            r_l = a_l.r;
            P_l = P[k][i-1][j];
            u_l = a_l.u/r_l;
            v_l = a_l.v/r_l;
            w_l = a_l.w/r_l;
            nx = 0.0;
            ny = 1.0;
            nz = 0.0;
          }
          else {
            a_l = a[j][k][i-1];
            r_l = a_l.r;
            P_l = P[j][k][i-1];
            u_l = a_l.u/r_l;
            v_l = a_l.v/r_l;
            w_l = a_l.w/r_l;
            nx = 0.0;
            ny = 0.0;
            nz = 1.0;
          }
          r_h = (r_l + r_r)/2.0;
          c_l = sqrt(gam*P_l/r_l);
          c_r = sqrt(gam*P_r/r_r);
          c_h = (c_l + c_r)/2.0;
          V = sqrt((u_l*u_l + v_l*v_l + w_l*w_l + u_r*u_r + v_r*v_r + w_r*w_r)/2.0);
          V_l = u_l*nx + v_l*ny + w_l*nz;
          V_r = u_r*nx + v_r*ny + w_r*nz;
          M = V_l/c_h;
          M_L = V_l/c_l;
          M_R = V_r/c_r;
          if (M >= 1.0){
            f_L = 1.0;
            f_R = 0.0;
          }
          else if (M <= -1.0){
            f_L = 0.0;
            f_R = 1.0;
          }
          else{
            f_L = (M + 1.0)*(M + 1.0)*(2.0 - M)/4.0;
            f_R = (M - 1.0)*(M - 1.0)*(2.0 + M)/4.0;
          }
          //g_HR = fmin(0.2,g_w);
					g_HR = 0.2;
          d_p = P_l - P_r;
          p = (P_l + P_r + (f_L - f_R)*d_p)/2.0 + g_HR*(f_L + f_R - 1.0)*r_h*c_h*V;
          M_c = fmin(1.0,V/c_h);
          chi = (1.0 - M_c)*(1.0 - M_c);
          g = 0.0 - fmax(fmin(M_L,0.0),-1.0)*fmin(fmax(M_R,0.0),1.0);
          V_n = (r_l*fabs(V_l) + r_r*fabs(V_r))/(r_l + r_r);
          V_np = (1.0 - g)*V_n + g*fabs(V_l);
          V_nn = (1.0 - g)*V_n + g*fabs(V_r);
          m_dot = 0.5*(r_l*(V_l + V_np) + r_r*(V_r - V_nn) - chi/c_h*d_p);

          // Flux calculation
          if (m_dot >= 0){
            fl[n][i][j][k] = (m_dot/r_l)*a_l;
            fl[n][i][j][k].u += P_l*nx;
            fl[n][i][j][k].v += P_l*ny;
            fl[n][i][j][k].w += P_l*nz;
            fl[n][i][j][k].E += P_l/r_l*m_dot;
          }
          else{
            fl[n][i][j][k] = m_dot*a_r;
            fl[n][i][j][k].r /= r_r;
            fl[n][i][j][k].u += P_r*nx;
            fl[n][i][j][k].v += P_r*ny;
            fl[n][i][j][k].w += P_r*nz;
            fl[n][i][j][k].E += P_r/r_r*m_dot;
          }

        }
      }
    }
  }
}

void Flux::fluxFacetoCell(){
  for (int i = 0; i < numdiv; i++){
    for (int j = 0; j < numdiv; j++){
      for (int k = 0; k < numdiv; k++){
        flux_c[i][j][k] = flux_f[0][i+1][j][k] - flux_f[0][i][j][k] + flux_f[1][j][k][i] - flux_f[1][j+1][i][k] + flux_f[2][k][i][j] - flux_f[2][k+1][i][j];
      }
    }
  }
}

void intFlux::wall(){
	for (int i = 0; i < numdiv; i++){
		for (int j = 0; j < numdiv; j++){
			fl[i][j] = 0.0;
		}
	}
}

void intFlux::inviscidFlux(){
  double r_l, r_r, r_h, P_l, P_r, g_mix_r, g_mix_l, C_p;
  //flow a_l, a_r;
  double c_l, c_r, c_h, chi, g, g_HR, nx, ny, nz;
  double u_r, u_l, v_r, v_l, w_l, w_r, V, V_l, V_r;
  double M, M_L, M_R, f_L, f_R;
  double d_p, p, M_c, V_n, V_np, V_nn, m_dot;

  //for(int n = 0; n < 3; n++){
    for (int i = 0; i < numdiv; i++){
      for (int j = 0; j < numdiv; j++){
        //for (int k = 0; k < numdiv; k++){
          u_r = a_r.u/r_r;
          v_r = a_r.v/r_r;
          w_r = a_r.w/r_r;
          /*if (n == 0){
            a_l = a[i-1][j];
            r_l = a_l.r;
            P_l = P[i-1][j];*/
            u_l = a_l.u/r_l;
            v_l = a_l.v/r_l;
            w_l = a_l.w/r_l;
            nx = 1.0;
            ny = 0.0;
            nz = 0.0;
          /*}
          else if (n == 1){
            a_l = a[k][i-1][j];
            r_l = a_l.r;
            P_l = P[k][i-1][j];
            u_l = a_l.u/r_l;
            v_l = a_l.v/r_l;
            w_l = a_l.w/r_l;
            nx = 0.0;
            ny = 1.0;
            nz = 0.0;
          }
          else {
            a_l = a[j][k][i-1];
            r_l = a_l.r;
            P_l = P[j][k][i-1];
            u_l = a_l.u/r_l;
            v_l = a_l.v/r_l;
            w_l = a_l.w/r_l;
            nx = 0.0;
            ny = 0.0;
            nz = 1.0;
          }*/
          r_h = (r_l + r_r)/2.0;
          c_l = sqrt(gam*P_l/r_l);
          c_r = sqrt(gam*P_r/r_r);
          c_h = (c_l + c_r)/2.0;
          V = sqrt((u_l*u_l + v_l*v_l + w_l*w_l + u_r*u_r + v_r*v_r + w_r*w_r)/2.0);
          V_l = u_l*nx + v_l*ny + w_l*nz;
          V_r = u_r*nx + v_r*ny + w_r*nz;
          M = V_l/c_h;
          M_L = V_l/c_l;
          M_R = V_r/c_r;
          if (M >= 1.0){
            f_L = 1.0;
            f_R = 0.0;
          }
          else if (M <= -1.0){
            f_L = 0.0;
            f_R = 1.0;
          }
          else{
            f_L = (M + 1.0)*(M + 1.0)*(2.0 - M)/4.0;
            f_R = (M - 1.0)*(M - 1.0)*(2.0 + M)/4.0;
          }
					//g_w =
          //g_HR = fmin(0.2,g_w);
					g_HR = 0.2;
          d_p = P_l - P_r;
          p = (P_l + P_r + (f_L - f_R)*d_p)/2.0 + g_HR*(f_L + f_R - 1.0)*r_h*c_h*V;
          M_c = fmin(1.0,V/c_h);
          chi = (1.0 - M_c)*(1.0 - M_c);
          g = 0.0 - fmax(fmin(M_L,0.0),-1.0)*fmin(fmax(M_R,0.0),1.0);
          V_n = (r_l*fabs(V_l) + r_r*fabs(V_r))/(r_l + r_r);
          V_np = (1.0 - g)*V_n + g*fabs(V_l);
          V_nn = (1.0 - g)*V_n + g*fabs(V_r);
          m_dot = 0.5*(r_l*(V_l + V_np) + r_r*(V_r - V_nn) - chi/c_h*d_p);

          // Flux calculation
          if (m_dot >= 0){
            fl[i][j] = (m_dot/r_l)*a_l;
            fl[i][j].u += P_l*nx;
            fl[i][j].v += P_l*ny;
            fl[i][j].w += P_l*nz;
            fl[i][j].E += P_l/r_l*m_dot;
          }
          else{
            fl[i][j] = m_dot*a_r;
            fl[i][j].r /= r_r;
            fl[i][j].u += P_r*nx;
            fl[i][j].v += P_r*ny;
            fl[i][j].w += P_r*nz;
            fl[i][j].E += P_r/r_r*m_dot;
          }

        //}
      }
    }
  //}
}

#include"main.def.h"

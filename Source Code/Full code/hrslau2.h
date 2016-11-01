#ifndef HRSLAU2_H_
#define HRSLAU2_H_

class hrslau2{
  public:
    void inviscidFlux(flow3D fl[], flow3D a, double3D p, int ndim){
      double r_l, r_r, r_h, P_l, P_r, g_mix_r, g_mix_l, C_p;
      vector<flow> a_l, a_r;
      double c_l, c_r, c_h, chi, g, g_HR, nx, ny nz;
      double u_r, u_l, v_r, v_l, w_l, w_r, V, V_l, V_r;
      double M, M_L, M_R, f_L, f_R;
      double d_p, p, M_c, V_n, V_np, V_nn, m_dot;

      for(int n = 0; n < 3; n++){
        for (int i = 1; i < ndim; i++){
          for (int j = 0; j < ndim; j++){
            for (int k = 0; k < ndim; k++){

              a_r = a[i][j][k];
              r_r = a_r.r;
              P_r = p[i][j][k];
              C_p = 0.0;
              for (int q = 0; q < a_r.Y.size(); q++){
                C_p += a_r.Y[q]*Cp[q];
              }
              g_mix_r = (C_p/(C_p - R));
              u_r = a_r.u/r_r;
              v_r = a_r.v/r_r;
              w_r = a_r.w/r_r;
              if (n == 0){
                a_l = a[i-1][j][k];
                r_l = a_l.r;
                P_l = p[i-1][j][k];
                C_p = 0.0;
                for (int q = 0; q < a_l.Y.size(); q++){
                  C_p += a_l.Y[q]*Cp[q];
                }
                g_mix_l = (C_p/(C_p - R));
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
                P_l = p[k][i-1][j];
                C_p = 0.0;
                for (int q = 0; q < a_l.Y.size(); q++){
                  C_p += a_l.Y[q]*Cp[q];
                }
                g_mix_l = (C_p/(C_p - R));
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
                P_l = p[j][k][i-1];
                C_p = 0.0;
                for (int q = 0; q < a_l.Y.size(); q++){
                  C_p += a_l.Y[q]*Cp[q];
                }
                g_mix_l = (C_p/(C_p - R));
                u_l = a_l.u/r_l;
                v_l = a_l.v/r_l;
                w_l = a_l.w/r_l;
                nx = 0.0;
                ny = 0.0;
                nz = 1.0;
              }
              r_h = (r_l + r_r)/2.0;
              c_l = sqrt(g_mix_l*P_l/r_l);
              c_r = sqrt(g_mix_r*P_r/r_r);
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
              g_HR = fmin(0.2,g_w);
              d_p = P_l - P_r;
              p = (P_l + P_r + (f_L - f_R)*d_p)/2.0 + g_HR*(f_L + f_R - 1.0)*r_h*c_h*V;
              M_c = fmin(1.0,V/c_h);
              chi = (1.0 - M_c)*(1.0 - M_c);
              g = 0.0 - fmax(fmin(M_L,0.0),-1.0)*fmin(fmax(M_R,0.0),1.0);
              V_n = (r_l*fabs(V_l) + r_r*fabs(V_r))/(r_l + r_r);
              V_np = (1.0 - g)*V_n + g*fabs(V_l);
              V_nn = (1.0 - g)*V_n + g*fabs(V_r);
              m_dot = 0.5*(r_l*(V_l + V_np) + r_r*(V_r - V_nn) - chi/c_h*d_p);

              if (m_dot >= 0){
                fl[n][i][j][k] = (m_dot/r_l)*a_l;
                fl[n][i][j][k].u += P_l*nx;
                fl[n][i][j][k].v += P_l*ny;
                fl[n][i][j][k].w += P_l*nz;
                fl[n][i][j][k].E += P_l/r_l*m_dot;
              }
              else{
                fl[n][i1][i2][i3] = m_dot*a_r;
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
}

#endif

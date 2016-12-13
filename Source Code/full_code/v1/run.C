#include "utility.h"
#include "run.h"

////////////////////////////////////////////////////////////
// Main Chare Functions
////////////////////////////////////////////////////////////

Main::Main(CkArgMsg* m){
	delete m;
	dimX = 5;
	dimY = 2;
	dimZ = 2;
	ndiv = 10;
	gma = 1.4;
	dx = double(1)/(dimX*ndiv);
	double end_time = 0.2;
	double Co = 0.1;
	double dt = Co*dx;
	t_steps = int(end_time/dt);
    
    // initialize readonly variables for chemical reactions
    double Te = 5.0*11604.0; //in K
    double Tg = 300.0; // in K
    double end_time_chem = 1.0e-8; // in s
    double dt_chem = 1.0e-15; // in s
    int iter_chem = int(end_time/dt);
    double R = 8.314; // Gas constat in J/mol.K
    double Av = 6.022e23; // Avogadro's number
    double n = 2.5e19; // Number density of air
    double eq = 0.1; // Equivalence Ratio
    int wf = 10000; // Frequency of writing output to file

  // call read_file();
    read_file();
  
	mainProxy = thisProxy;
	cellProxy = CProxy_Cell::ckNew(dimX,dimY,dimZ);
	fluxProxy = CProxy_Flux::ckNew(dimX,dimY,dimZ);
	int dimW = 3;
	CkArrayOptions opts(CkArrayIndex4D(0,0,0,0),
											CkArrayIndex4D(dimW,dimX+1,dimX+1,dimX+1),
											CkArrayIndex4D(1,1,1,1));
	interfaceProxy = CProxy_Interface::ckNew(opts);
	interfaceProxy.solve_i();
}

void Main::InterfaceIsUp(){
	CkPrintf("interface is up \n");
	cellProxy.solve_c();
}

void Main::done(){
	CkExit();
}

void Main::read_file(){
    std::ifstream myfile;
    std::string line;
    int temp, numreact, numproduct, temp1;
    double temp2, l1, l2, l3;
    std::string name;
    string1D secspecies;
    double1D sec, data, secmp;
    //data.resize(5);
    std::string a1, a2, a3, a4, a5, a6, p1, p2 , p3, p4, p5, asname;
    //double as;
    myfile.open("methane_air_plasma.txt");
    std::cout << "Reading file \n";
    getline(myfile,line);
    getline(myfile,line);
    std::istringstream iss1(line);
    iss1 >> size;
    std::cout << "Number of species : " << size << "\n";
    species.resize(size);
    Cp.resize(size);
    sp.resize(size);
    p_rxn.resize(size);
    r_p.resize(size);
    H_f.resize(size);
    getline(myfile,line);
    getline(myfile,line);
    for (int i = 0; i < size; i++){
        getline(myfile,line);
        std::istringstream iss2(line);
        iss2 >> temp >> species[i] >> Cp[i] >> H_f[i];
    }
    std::cout << "Thermodynamics reading complete \n";
    getline(myfile,line);
    getline(myfile,line);
    getline(myfile,line);
    std::istringstream iss(line);
    iss >> rxn_size;
    std::cout << "Number of reactions : " << rxn_size << "\n";
    K.resize(rxn_size);
    rs.resize(rxn_size);
    adv.resize(rxn_size);
    add_info.resize(rxn_size);
    tb_sp.resize(rxn_size);
    H_f.resize(rxn_size);
    for (int i = 0; i < rxn_size; i++){
        std::cout << "Reading reaction " << i+1 << "\n";
        getline(myfile,line);
        std::istringstream iss(line);
        iss >> temp >> numreact >> numproduct;
        temp = 0;
        if (numreact == 1){
            iss >> a1 >> a2;
            for (int j = 0; j < size; j++){
                if (a1 == species[j]){
                    rs[i].push_back(j);
                    r_p[j].push_back(-1.0);
                    p_rxn[j].push_back(i);
                }
            }
        }
        else if (numreact == 2){
            iss >> a1 >> a2 >> a3 >> a4;
            if (a3 == "M"){
                temp = 1;
            }
            for (int j = 0; j < size; j++){
                if (a1 == species[j]){
                    rs[i].push_back(j);
                    r_p[j].push_back(-1.0);
                    p_rxn[j].push_back(i);
                }
                if (a3 == species[j]){
                    rs[i].push_back(j);
                    r_p[j].push_back(-1.0);
                    p_rxn[j].push_back(i);
                }
            }
        }
        else{
            iss >> a1 >> a2 >> a3 >> a4 >> a5 >> a6;
            if (a5 == "(M)"){
                temp = 2;
            }
            if (a5 == "M"){
                temp = 1;
            }
            for (int j = 0; j < size; j++){
                if (a1 == species[j]){
                    rs[i].push_back(j);
                    r_p[j].push_back(-1.0);
                    p_rxn[j].push_back(i);
                }
                if (a3 == species[j]){
                    rs[i].push_back(j);
                    r_p[j].push_back(-1.0);
                    p_rxn[j].push_back(i);
                }
                if (a5 == species[j]){
                    rs[i].push_back(j);
                    r_p[j].push_back(-1.0);
                    p_rxn[j].push_back(i);
                }
            }
        }
        if (numproduct == 1){
            iss >> p1;
            for (int j = 0; j < size; j++){
                if (p1 == species[j]){
                    r_p[j].push_back(1.0);
                    p_rxn[j].push_back(i);
                }
            }
        }
        else if (numproduct == 2){
            iss >> p1 >> p2 >> p3;
            for (int j = 0; j < size; j++){
                if (p1 == species[j]){
                    r_p[j].push_back(1.0);
                    p_rxn[j].push_back(i);
                }
                if (p3 == species[j]){
                    r_p[j].push_back(1.0);
                    p_rxn[j].push_back(i);
                }
            }
        }
        else{
            iss >> p1 >> p2 >> p3 >> p4 >> p5;
            for (int j = 0; j < size; j++){
                if (p1 == species[j]){
                    r_p[j].push_back(1.0);
                    p_rxn[j].push_back(i);
                }
                if (p3 == species[j]){
                    r_p[j].push_back(1.0);
                    p_rxn[j].push_back(i);
                }
                if (p5 == species[j]){
                    r_p[j].push_back(1.0);
                    p_rxn[j].push_back(i);
                }
            }
        }
        data.resize(5);
        secspecies.resize(0);
        secmp.resize(0);
        iss >> data[0] >> data[1] >> data[2] >> data[3] >> data[4];
        data[2] = double(0) - data[2];
        data[4] = double(0) - data[4];
        K[i].push_back(data[0]);
        K[i].push_back(data[1]);
        K[i].push_back(data[2]);
        K[i].push_back(data[3]);
        K[i].push_back(data[4]);
        if (temp == 2){
            tb_sp[i].resize(size);
            add_info[i].resize(size+4);
            iss >> l1 >> l2 >> l3;
            add_info[i][0] = 1.0;
            add_info[i][1] = l1;
            add_info[i][2] = l2;
            add_info[i][3] = l3;
            iss >> temp1;
            for (int j = 0; j < temp1; j++){
                add_info[i][j+4] = 1.0;
                iss >> name >> temp2;
                secspecies.push_back(name);
                secmp.push_back(temp2);
            }
            for (int j = 0; j < size; j++){
                for (unsigned int k = 0; k < secspecies.size(); k++){
                    if (secspecies[k] == species[j]){
                        add_info[i][j+4] = secmp[k];
                    }
                }
            }
        }
        if (temp == 1){
            tb_sp[i].resize(size);
            add_info[i].resize(size+4);
            add_info[i][0] = 0.0;
            add_info[i][1] = 0.0;
            add_info[i][2] = 0.0;
            add_info[i][3] = 0.0;
            iss >> temp1;
            for (int j = 0; j < temp1; j++){
                add_info[i][j+4] = 1.0;
                iss >> name >> temp2;
                secspecies.push_back(name);
                secmp.push_back(temp2);
            }
            for (int j = 0; j < size; j++){
                for (unsigned int k = 0; k < secspecies.size(); k++){
                    if (secspecies[k] == species[j]){
                        add_info[i][j+4] = secmp[k];
                    }
                }
            }
        }
        else{
            add_info[i].push_back(0.0);
            add_info[i].push_back(0.0);
            add_info[i].push_back(0.0);
            add_info[i].push_back(0.0);
        }
        data.clear();
        secspecies.clear();
        secmp.clear();
    }
    myfile.close();
    
    std::ofstream myfile1;
    myfile1.open("Output.csv");
    myfile1 << "t (s)" << "," << "T_g (K)";
    for (int i = 0; i < size; i++){
        myfile1 << "," << species[i];
    }
    myfile1 << "\n";
    myfile1.close();
} // end read_file()

//////////////////////////////////////////////////////////////////////
// Cell [3D] Chare Array Functions
//////////////////////////////////////////////////////////////////////

Cell::Cell(){
	val_new.resize(ndiv);
	val_old.resize(ndiv);
	P.resize(ndiv);
	for (int i = 0; i < ndiv; i++){
		val_new[i].resize(ndiv);
		val_old[i].resize(ndiv);
		P[i].resize(ndiv);
		for (int j = 0; j < ndiv; j++){
			val_new[i][j].resize(ndiv);
			val_old[i][j].resize(ndiv);
			P[i][j].resize(ndiv);
		}
	}
	initialize();
}

void Cell::gaslaw(){
  for (int i = 0; i < ndiv; i++){
    for (int j = 0; j < ndiv; j++){
      for (int k = 0; k < ndiv; k++){
        P[i][j][k] = (gma-1)*val_new[i][j][k].r*(val_new[i][j][k].E - 0.5*(val_new[i][j][k].u*val_new[i][j][k].u + val_new[i][j][k].v*val_new[i][j][k].v + val_new[i][j][k].w*val_new[i][j][k].w));
      }
    }
  }
}

void Cell::update(){
	for (int i = 0; i < ndiv; i++){
		for (int j = 0; j < ndiv; j++){
			for (int k = 0; k < ndiv; k++){
				val_old[i][j][k] = val_new[i][j][k];
			}
		}
	}
}

void Cell::calcvar3D(flow3D &v_n, flow3D v_o, flow3D fl){
  for (int i = 0; i < ndiv; i++){
    for(int j = 0; j < ndiv; j++){
      for(int k = 0; k < ndiv; k++){
        v_n[i][j][k] = v_o[i][j][k] - dt *  fl[i][j][k] / dx;
      }
    }
  }
}

void Cell::initialize(){
	double P_i, r_i;
	if (thisIndex.x < dimX/2){
		P_i = 1.0;
		r_i = 1.0;
	}
	else{
		P_i = 0.1;
		r_i = 0.125;
	}
	for (int i = 0; i < ndiv; i++){
		for (int j = 0; j < ndiv; j++){
			for (int k = 0; k < ndiv; k++){
				P[i][j][j] = P_i;
				val_old[i][j][k].r = r_i;
				val_old[i][j][k].u = 0.0;
				val_old[i][j][k].v = 0.0;
				val_old[i][j][k].w = 0.0;
				val_old[i][j][k].E = P_i/(r_i*(gma-1));
			}
		}
	}
}

Flux::Flux(){
	flux_c.resize(ndiv);
	cell_val.resize(ndiv);
	P.resize(ndiv);
	for (int i = 0; i < ndiv; i++){
		flux_c[i].resize(ndiv);
		cell_val[i].resize(ndiv);
		P[i].resize(ndiv);
		for (int j = 0; j < ndiv; j++){
			flux_c[i][j].resize(ndiv);
			cell_val[i][j].resize(ndiv);
			P[i][j].resize(ndiv);
		}
	}
	flux_f.resize(3);
	for (int i = 0; i < 3; i++){
		flux_f[i].resize(ndiv+1);
		for (int j = 0; j <= ndiv; j++){
			flux_f[i][j].resize(ndiv);
			for (int k = 0; k < ndiv; k++){
				flux_f[i][j][k].resize(ndiv);
			}
		}
	}
}

void Flux::inviscidFlux(){
	double r_l, r_r, r_h, P_l, P_r;
  flow a_l, a_r;
  double c_l, c_r, c_h, chi, g, nx, ny, nz;
  double u_r, u_l, v_r, v_l, w_l, w_r, V, V_l, V_r;
  double M, M_L, M_R, f_L, f_R;
  double d_p, p, M_c, V_n, V_np, V_nn, m_dot;

	for(int n = 0; n < 3; n++){
    for (int i = 1; i < ndiv; i++){
      for (int j = 0; j < ndiv; j++){
        for (int k = 0; k < ndiv; k++){

					a_r = cell_val[i][j][k];
					r_r = a_r.r;
          P_r = P[i][j][k];
          u_r = a_r.u/r_r;
          v_r = a_r.v/r_r;
          w_r = a_r.w/r_r;
          if (n == 0){
            a_l = cell_val[i-1][j][k];
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
            a_l = cell_val[k][i-1][j];
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
            a_l = cell_val[j][k][i-1];
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
          c_l = sqrt(gma*P_l/r_l);
          c_r = sqrt(gma*P_r/r_r);
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
          d_p = P_l - P_r;
          p = (P_l + P_r + (f_L - f_R)*d_p)/2.0 + 0.5*(f_L + f_R - 1.0)*r_h*c_h*V;
          M_c = fmin(1.0,V/c_h);
          chi = (1.0 - M_c)*(1.0 - M_c);
          g = 0.0 - fmax(fmin(M_L,0.0),-1.0)*fmin(fmax(M_R,0.0),1.0);
          V_n = (r_l*fabs(V_l) + r_r*fabs(V_r))/(r_l + r_r);
          V_np = (1.0 - g)*V_n + g*fabs(V_l);
          V_nn = (1.0 - g)*V_n + g*fabs(V_r);
          m_dot = 0.5*(r_l*(V_l + V_np) + r_r*(V_r - V_nn) - chi/c_h*d_p);

          // Flux calculation
          if (m_dot >= 0){
            flux_f[n][i][j][k] = (m_dot/r_l)*a_l;
            flux_f[n][i][j][k].u += P_l*nx;
            flux_f[n][i][j][k].v += P_l*ny;
            flux_f[n][i][j][k].w += P_l*nz;
            flux_f[n][i][j][k].E += P_l/r_l*m_dot;
          }
          else{
            flux_f[n][i][j][k] = m_dot*a_r;
            flux_f[n][i][j][k].r /= r_r;
            flux_f[n][i][j][k].u += P_r*nx;
            flux_f[n][i][j][k].v += P_r*ny;
            flux_f[n][i][j][k].w += P_r*nz;
            flux_f[n][i][j][k].E += P_r/r_r*m_dot;
          }

				}
			}
		}
	}
}

void Flux::fluxFacetoCell(){
  for (int i = 0; i < ndiv; i++){
    for (int j = 0; j < ndiv; j++){
      for (int k = 0; k < ndiv; k++){
        flux_c[i][j][k] = flux_f[0][i+1][j][k] - flux_f[0][i][j][k] + flux_f[1][j+1][k][i] - flux_f[1][j][i][k] + flux_f[2][k+1][i][j] - flux_f[2][k][i][j];
      }
    }
  }
}

Interface::Interface(){
//	CkPrintf("Interface being created with index (%d,%d,%d,%d) \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z  );
	flux.resize(ndiv);
	val_l.resize(ndiv);
	val_r.resize(ndiv);
	P_left.resize(ndiv);
	P_right.resize(ndiv);
	for (int i = 0; i < ndiv; i++){
		flux[i].resize(ndiv);
		val_l[i].resize(ndiv);
		val_r[i].resize(ndiv);
		P_left[i].resize(ndiv);
		P_right[i].resize(ndiv);
	}
}

void Interface::calc(){
	double r_l, r_r, r_h, P_l, P_r;
  flow a_l, a_r;
  double c_l, c_r, c_h, chi, g, nx, ny, nz;
  double u_r, u_l, v_r, v_l, w_l, w_r, V, V_l, V_r;
  double M, M_L, M_R, f_L, f_R;
  double d_p, p, M_c, V_n, V_np, V_nn, m_dot;

	for (int i = 0; i < ndiv; i++){
		for (int j = 0; j < ndiv; j++){

			a_r = val_r[i][j];
			r_r = a_r.r;
			P_r = P_right[i][j];
			u_r = a_r.u/r_r;
			v_r = a_r.v/r_r;
			w_r = a_r.w/r_r;
			a_l = val_l[i][j];
			r_l = a_l.r;
			P_l = P_left[i][j];
			u_l = a_l.u/r_l;
			v_l = a_l.v/r_l;
			w_l = a_l.w/r_l;
			if (thisIndex.w == 0){
				nx = 1.0;
				ny = 0.0;
				nz = 0.0;
			}
			else if (thisIndex.w == 1){
				nx = 0.0;
				ny = 1.0;
				nz = 0.0;
			}
			else {
				nx = 0.0;
				ny = 0.0;
				nz = 1.0;
			}
			r_h = (r_l + r_r)/2.0;
			c_l = sqrt(gma*P_l/r_l);
			c_r = sqrt(gma*P_r/r_r);
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
			d_p = P_l - P_r;
			p = (P_l + P_r + (f_L - f_R)*d_p)/2.0 + 0.5*(f_L + f_R - 1.0)*r_h*c_h*V;
			M_c = fmin(1.0,V/c_h);
			chi = (1.0 - M_c)*(1.0 - M_c);
			g = 0.0 - fmax(fmin(M_L,0.0),-1.0)*fmin(fmax(M_R,0.0),1.0);
			V_n = (r_l*fabs(V_l) + r_r*fabs(V_r))/(r_l + r_r);
			V_np = (1.0 - g)*V_n + g*fabs(V_l);
			V_nn = (1.0 - g)*V_n + g*fabs(V_r);
			m_dot = 0.5*(r_l*(V_l + V_np) + r_r*(V_r - V_nn) - chi/c_h*d_p);

			// Flux calculation
			if (m_dot >= 0){
				flux[i][j] = (m_dot/r_l)*a_l;
				flux[i][j].u += P_l*nx;
				flux[i][j].v += P_l*ny;
				flux[i][j].w += P_l*nz;
				flux[i][j].E += P_l/r_l*m_dot;
			}
			else{
				flux[i][j] = m_dot*a_r;
				flux[i][j].r /= r_r;
				flux[i][j].u += P_r*nx;
				flux[i][j].v += P_r*ny;
				flux[i][j].w += P_r*nz;
				flux[i][j].E += P_r/r_r*m_dot;
			}

		}
	}

}

void Interface::wall(flow2D &f_n, flow2D f_o, double2D &P_n, double2D P_o){
	for (int i = 0; i < ndiv; i++){
		for (int j = 0; j < ndiv; j++){
			f_n[i][j] = f_o[i][j];
			P_n[i][j] = P_o[i][j];
			if (thisIndex.w == 0){
				f_n[i][j].u = 0.0 - f_o[i][j].u;
			}
			else if (thisIndex.w == 1){
				f_n[i][j].v = 0.0 - f_o[i][j].v;
			}
			else{
				f_n[i][j].w = 0.0 - f_o[i][j].w;
			}
		}
	}
}

#include "run.def.h"

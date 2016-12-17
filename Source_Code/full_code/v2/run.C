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
	double end_time = 0.1;
	double Co = 0.1;
	double dt = Co*dx;
	t_steps = int(end_time/dt);
    CkPrintf("t_steps: %d\n", t_steps);

    // initialize readonly variables for chemical reactions
    w_f = 320;
    Te = 5.0*11604.0; //in K
    //Tg = 300.0; // in K
		conc_i = 41.5;
		conc_f = 41.5;
    end_time_chem = dt; // in s, -8 originally
    dt_chem = 1.0e-15; // in s, 1.0e-15 originally
    iter_chem = int(end_time_chem / dt_chem);
    //CkPrintf("iter_chem = %d\n", iter_chem);
    R = 8.314; // Gas constat in J/mol.K
    Av = 6.022e23; // Avogadro's number
    n = 2.5e19; // Number density of air
    eq = 0.1; // Equivalence Ratio
    wf = 10000; // Frequency of writing output to file

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
  myfile.open("rxn.txt");
	std::cout << "Reading file \n";
	getline(myfile,line);
	getline(myfile,line);
	std::istringstream iss1(line);
	iss1 >> size;
	std::cout << "Number of species : " << size << "\n";
	species.resize(size);
	Cp.resize(size);
	//sp.resize(size);
	p_rxn.resize(size);
	r_p.resize(size);
	H_f.resize(size);
	S_f.resize(size);
	getline(myfile,line);
	getline(myfile,line);
	for (int i = 0; i < size; i++){
		getline(myfile,line);
		std::istringstream iss2(line);
		iss2 >> temp >> species[i] >> Cp[i] >> H_f[i] >> S_f[i];
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
    double2D ps;//v2
	ps.resize(rxn_size);
	//adv.resize(rxn_size); //v2
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
					ps[i].push_back(j);
					r_p[j].push_back(1.0);
					p_rxn[j].push_back(i);
				}
			}
		}
		else if (numproduct == 2){
			iss >> p1 >> p2 >> p3;
			for (int j = 0; j < size; j++){
				if (p1 == species[j]){
					ps[i].push_back(j);
					r_p[j].push_back(1.0);
					p_rxn[j].push_back(i);
				}
				if (p3 == species[j]){
					ps[i].push_back(j);
					r_p[j].push_back(1.0);
					p_rxn[j].push_back(i);
				}
			}
		}
		else{
			iss >> p1 >> p2 >> p3 >> p4 >> p5;
			for (int j = 0; j < size; j++){
				if (p1 == species[j]){
					ps[i].push_back(j);
					r_p[j].push_back(1.0);
					p_rxn[j].push_back(i);
				}
				if (p3 == species[j]){
					ps[i].push_back(j);
					r_p[j].push_back(1.0);
					p_rxn[j].push_back(i);
				}
				if (p5 == species[j]){
					ps[i].push_back(j);
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
				//add_info[i][j+4] = 0.0;
				iss >> name >> temp2;
				secspecies.push_back(name);
				secmp.push_back(temp2);
			}
			for (int j = 0; j < size; j++){
				add_info[i][j+4] = 0.0;
				for (unsigned int k = 0; k < secspecies.size(); k++){
					if (secspecies[k] == species[j]){
						add_info[i][j+4] = secmp[k];
					}
				}
			}
			add_info[i][23] = 1.0;
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
				//add_info[i][j+4] = 1.0;
				iss >> name >> temp2;
				secspecies.push_back(name);
				secmp.push_back(temp2);
			}
			for (int j = 0; j < size; j++){
				add_info[i][j+4] = 0.0;
				for (unsigned int k = 0; k < secspecies.size(); k++){
					if (secspecies[k] == species[j]){
						add_info[i][j+4] = secmp[k];
					}
				}
			}
			add_info[i][23] = 1.0;
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
}

//////////////////////////////////////////////////////////////////////
// Cell [3D] Chare Array Functions
//////////////////////////////////////////////////////////////////////

Cell::Cell(){
    sp.resize(ndiv);//v2
	val_new.resize(ndiv);
	val_old.resize(ndiv);
	P.resize(ndiv);
	adv.resize(ndiv);
	Tg.resize(ndiv);
	for (int i = 0; i < ndiv; i++){
        sp[i].resize(ndiv);//v2
		val_new[i].resize(ndiv);
		val_old[i].resize(ndiv);
		P[i].resize(ndiv);
		adv[i].resize(ndiv);
		Tg[i].resize(ndiv);
		for (int j = 0; j < ndiv; j++){
            sp[i][j].resize(ndiv);//v2
			val_new[i][j].resize(ndiv);
			val_old[i][j].resize(ndiv);
			P[i][j].resize(ndiv);
			adv[i][j].resize(ndiv);
			Tg[i][j].resize(ndiv);
			for (int k = 0; k < ndiv; k++){
                sp[i][j][k].resize(size);//v2
				adv[i][j][k].resize(size);
			}
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
        v_n[i][j][k] = v_o[i][j][k] - dt*fl[i][j][k]/dx + S[i][j][k]*dt;
      }
    }
  }
}

void Cell::initialize(){
	for (int i = 0; i < ndiv; i++){
		for (int j = 0; j < ndiv; j++){
			for (int k = 0; k < ndiv; k++){
				P[i][j][j] = P_i;
				val_old[i][j][k].r = r_i;
				val_old[i][j][k].u = u_i;
				val_old[i][j][k].v = v_i;
				val_old[i][j][k].w = w_i;
				val_old[i][j][k].E = E_i;
				for (int l = 0; l < val_old[i][j][k].Y.size(); l++){
					val_old[i][j][k].Y[l] = 0.0;
				}
				val_old[i][j][k].Y[0] = 0.2095*conc_i;
				val_old[i][j][k].Y[1] = 0.7809*conc_i;
				val_old[i][j][k].Y[2] = 0.0093*conc_i;
				val_old[i][j][k].Y[3] = 0.0003*conc_i;
			}
		}
	}
}

void Cell::initialize_chem(){
    ds = std::vector<std::vector<std::vector<std::vector<double>>>>(ndiv, std::vector<std::vector<std::vector<double>>>(ndiv, std::vector<std::vector<double>>(size, std::vector<double>(ndiv, 0)))); // v2

	for (int i = 0; i < ndiv; i++){
    for(int j = 0; j < ndiv; j++){
      for(int k = 0; k < ndiv; k++){

				for(int l = 0; l < val_old[i][j][k].Y.size(); l++){
	        sp[i][j][k][l] = val_old[i][j][k].Y[l];
	      }
      }
    }
  }
}

void Cell::solve_rxn(){
    double4D sp1, sp_temp;
    double4D k, k1, k2, k3, k4;
		sp1.resize(ndiv);
		sp_temp.resize(ndiv);
		k1.resize(ndiv);
		k2.resize(ndiv);
		k3.resize(ndiv);
		k4.resize(ndiv);
		double dnH, n_total;
        double3D Cp_mix; // v2
		for (int i = 0; i < ndiv; i++){
			sp1[i].resize(ndiv);
			sp_temp[i].resize(ndiv);
			k1[i].resize(ndiv);
			k2[i].resize(ndiv);
			k3[i].resize(ndiv);
			k4[i].resize(ndiv);
	    for (int j = 0; j < ndiv; j++){
				sp1[i][j].resize(ndiv);
				sp_temp[i][j].resize(ndiv);
				k1[i][j].resize(ndiv);
				k2[i][j].resize(ndiv);
				k3[i][j].resize(ndiv);
				k4[i][j].resize(ndiv);
	      for (int k = 0; k < ndiv; k++){
					sp1[i][j][k].resize(size);
			    sp_temp[i][j][k].resize(size);
			    k1[i][j][k].resize(size);
			    k2[i][j][k].resize(size);
			    k3[i][j][k].resize(size);
			    k4[i][j][k].resize(size);
	      }
	    }
	  }
    initialize_chem();
    //write_file(-1);
    //CkPrintf("iter_chem = %d\n", iter_chem);
    for (int i = 0; i < iter_chem; i++){
			calc_change(k1,sp);
			for (int x = 0; x < ndiv; x++){
				for (int y = 0; y < ndiv; y++){
					for (int z = 0; z < ndiv; z++){
						sp_temp[x][y][z] = sp[x][y][z] + k1[x][y][z]*(dt_chem/double(2));
					}
				}
			}
			calc_change(k2,sp_temp);
			for (int x = 0; x < ndiv; x++){
				for (int y = 0; y < ndiv; y++){
					for (int z = 0; z < ndiv; z++){
						sp_temp[x][y][z] = sp[x][y][z] + k2[x][y][z]*(dt_chem/double(2));
					}
				}
			}
			calc_change(k3,sp_temp);
			for (int x = 0; x < ndiv; x++){
				for (int y = 0; y < ndiv; y++){
					for (int z = 0; z < ndiv; z++){
						sp_temp[x][y][z] = sp[x][y][z] + k3[x][y][z] * dt_chem;
					}
				}
			}
			calc_change(k4,sp_temp);
			for (int x = 0; x < ndiv; x++){
				for (int y = 0; y < ndiv; y++){
					for (int z = 0; z < ndiv; z++){
						k[x][y][z] = k1[x][y][z]/double(6) + k2[x][y][z]/double(3) + k3[x][y][z]/double(3) + k4[x][y][z]/double(6);
						ds[x][y][z] = ds[x][y][z] + k[x][y][z];
						sp[x][y][z] = sp[x][y][z] + k[x][y][z] * dt_chem;
						for (int i = 0; i < sp[0][0][0].size(); i++){
							S[x][y][z].Y[i] = ds[x][y][z][i];
						}
						dnH = 0.0;
				    n_total = 0.0;
				    //Cp_mix = 0.0;
				    for (int i = 0; i < size; i++){
				        dnH += (k[x][y][z][i]*H_f[i]);
				        n_total += sp[x][y][z][i];
				        Cp_mix[x][y][z] += (sp[x][y][z][i]*Cp[i]);
				    }
				    Cp_mix[x][y][z] /= n_total;
						S[x][y][z].E = dnH;
					}
				}
			}
    }
}

void Cell::calc_change(double4D& t_k, double4D& t_s){
    double k_o, k_f, k_inf, Pr, tb_mp, conc;
		for (int x = 0; x < ndiv; x++){
			for (int y = 0; y < ndiv; y++){
				for (int z = 0; z < ndiv; z++){
					for (int i = 0; i < rxn_size; i++){
			        k_o = K[i][0]*pow(Tg[x][y][z],K[i][1])*exp(K[i][2]/(R*Tg[x][y][z]))*pow(300.0/Te,K[i][3])*exp(K[i][4]/(R*Te));
			        tb_mp = 0.0;
			        for (unsigned int j = 0; j < tb_sp[i].size(); j++){
			            tb_mp += add_info[i][j+4]*t_s[x][y][z][tb_sp[i][j]];
			        }
			        if (add_info[i][0] == double(1)){
			            k_inf = add_info[i][1]*pow(Tg[x][y][z],add_info[i][2])*exp(add_info[i][3]/(R*Tg[x][y][z]));
			            Pr = k_o*tb_mp/k_inf;
			            k_f = k_inf * Pr/(1.0 + Pr);
			        }
			        else{
			            k_f = k_o;
			        }
			        conc = 1.0;
			        for (unsigned int j = 0; j < rs[i].size(); j++){
			            conc *= t_s[x][y][z][rs[i][j]];
			        }
			        if (tb_sp[i].size() == 0){
			            tb_mp = 1.0;
			        }
			        adv[x][y][z][i] = k_f * conc * tb_mp;
			    }
			    for (int i = 0; i < size; i++){
			        t_k[x][y][z][i] = 0.0;
			        for (unsigned int j = 0; j < p_rxn[i].size(); j++){
			            t_k[x][y][z][i] += (r_p[i][j]*adv[x][y][z][p_rxn[i][j]]);
			        }
			    }
				}
			}
		}
}


//////////////////////////////////////////////////////////////////////////
// Flux [3D] Chare Array Functions
//////////////////////////////////////////////////////////////////////////

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


////////////////////////////////////////////////////////////////////
// Interface [4D] Chare Array Functions
////////////////////////////////////////////////////////////////////

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

void Interface::inlet(flow2D &f_n, flow2D f_o, double2D &P_n, double2D P_o){
	for (int i = 0; i < ndiv; i++){
		for (int j = 0; j < ndiv; j++){
			double M = sqrt(f_o[i][j].u*f_o[i][j].u + f_o[i][j].v*f_o[i][j].v + f_o[i][j].w*f_o[i][j].w)*f_o[i][j].r/(gma*P_o[i][j]);
			if (M >= 1.0){
				P_n[i][j] = P_i;
				f_n[i][j].r = r_i;
				f_n[i][j].u = u_i;
				f_n[i][j].v = v_i;
				f_n[i][j].w = w_i;
				f_n[i][j].E = E_i;
			}
			else{
				P_n[i][j] = P_o[i][j];
				double V = sqrt((pow(Pt/P_n[i][j],(gma-1)/gma)-1)*2.0*gma*P_n[i][j]/((gma-1)*f_o[i][j].r));
				double M_n = sqrt(gma*P_n[i][j]/f_n[i][j].r);
				f_n[i][j].r = rt/pow((1+(gma-1)/2*M_n*M_n),1/(gma-1));
				f_n[i][j].u = V;
				f_n[i][j].v = 0.0;
				f_n[i][j].w = 0.0;
				f_n[i][j].E = P_n[i][j]/(f_n[i][j].r*(gma-1)) + 0.5*V*V;
			}
			for (int k = 4; k < f_o[i][j].Y.size(); k++){
				f_n[i][j].Y[k] = 0.0;
			}
			f_n[i][j].Y[0] = 0.2095*conc_i;
			f_n[i][j].Y[1] = 0.7809*conc_i;
			f_n[i][j].Y[2] = 0.0093*conc_i;
			f_n[i][j].Y[3] = 0.0003*conc_i;
		}
	}
}

void Interface::outlet(flow2D &f_n, flow2D f_o, double2D &P_n, double2D P_o){
	for (int i = 0; i < ndiv; i++){
		for (int j = 0; j < ndiv; j++){
			f_n[i][j] = f_o[i][j];
			P_n[i][j] = P_o[i][j];
			double M = sqrt(f_o[i][j].u*f_o[i][j].u + f_o[i][j].v*f_o[i][j].v + f_o[i][j].w*f_o[i][j].w)*f_o[i][j].r/(gma*P_o[i][j]);
			if (M >= 1.0){
				P_n[i][j] = P_o[i][j];
				f_n[i][j].E = P_n[i][j]/(f_n[i][j].r*(gma-1)) + 0.5*(f_n[i][j].u*f_n[i][j].u + f_n[i][j].v*f_n[i][j].v + f_n[i][j].w*f_n[i][j].w);
			}
		}
	}
}

void Interface::fuelinlet(flow2D &f_n, flow2D f_o, double2D &P_n, double2D P_o){
	for (int i = 0; i < ndiv; i++){
		for (int j = 0; j < ndiv; j++){
			double M = sqrt(f_o[i][j].u*f_o[i][j].u + f_o[i][j].v*f_o[i][j].v + f_o[i][j].w*f_o[i][j].w)*f_o[i][j].r/(gma*P_o[i][j]);
			if (M >= 1.0){
				P_n[i][j] = P_f;
				f_n[i][j].r = r_f;
				f_n[i][j].u = u_f;
				f_n[i][j].v = v_f;
				f_n[i][j].w = w_f;
				f_n[i][j].E = E_f;
			}
			else{
				P_n[i][j] = P_o[i][j];
				double V = sqrt((pow(Pt_f/P_n[i][j],(gma-1)/gma)-1)*2.0*gma*P_n[i][j]/((gma-1)*f_o[i][j].r));
				f_n[i][j].r = rt_f/pow((1+(gma-1)/2*M*M),1/(gma-1));
				f_n[i][j].u = V;
				f_n[i][j].v = 0.0;
				f_n[i][j].w = 0.0;
				f_n[i][j].E = P_n[i][j]/(f_n[i][j].r*(gma-1)) + 0.5*V*V;
			}
			for (int k = 0; k < f_o[i][j].Y.size(); k++){
				f_n[i][j].Y[k] = 0.0;
			}
			f_n[i][j].Y[4] = conc_f;
		}
	}
}

#include "run.def.h"

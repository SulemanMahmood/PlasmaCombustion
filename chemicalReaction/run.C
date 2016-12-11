#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "utility.h"

//std::string species_info = ;
double Te = 5.0*11604.0; //in K
double Tg = 2000.0; // in K
double end_time = 1.0e3; // in s
double dt = 1.0e-2; // in s
int iter = int(end_time/dt);
double R = 8.314; // Gas constat in J/mol.K
double Av = 6.022e23; // Avogadro's number
double n = 2.5e19; // Number density of air
double eq = 1; // Equivalence Ratio

double3D K;
int size;
string1D species;
double1D sp;
double1D E;
double1D Cp;
int wf = 1000; // Frequency of writing output to file

void read_file();
void solve_rxn();
void calc_change(double1D&, double1D&);
void write_file(int);
void calc_temp(double1D&);
void initialize();

int main(){
  read_file();
  solve_rxn();
  return 0;
}

void read_file(){
	int react1, react2, react3, prod1, prod2, prod3;
  std::ifstream myfile;
	std::string line;
	int temp, numreact, numproduct, t;
	string1D secspecies;
	double1D sec, data;
	data.resize(5);
	std::string a1, a2, a3, a4, a5, a6, p1, p2 , p3, p4, p5, asname;
	//double as;
  myfile.open("methane_air_plasma_without_3rdbody.txt");
	std::cout << "Reading file \n";
	getline(myfile,line);
	//while(!myfile.eof()){
		getline(myfile,line);
		std::istringstream iss1(line);
		iss1 >> size;
		std::cout << "Number of species : " << size << "\n";
		species.resize(size);
		Cp.resize(size);
		E.resize(size);
		sp.resize(size);
		K.resize(size);
		getline(myfile,line);
		getline(myfile,line);
		for (int i = 0; i < size; i++){
			getline(myfile,line);
			std::istringstream iss2(line);
			iss2 >> temp >> species[i] >> Cp[i] >> E[i];
		}
		std::cout << "Thermodynamics reading complete \n";
		getline(myfile,line);
		getline(myfile,line);
		getline(myfile,line);
		std::istringstream iss(line);
		iss >> t;
		std::cout << "Number of reactions : " << t << "\n";
		for (int i = 0; i < t; i++){
			getline(myfile,line);
			std::istringstream iss(line);
			iss >> temp >> numreact >> numproduct;
			temp = 0;
			if (numreact == 1){
				iss >> a1 >> a2;
				for (int j = 0; j < size; j++){
					if (a1 == species[j]){
						react1 = j;
					}
				}
			}
			else if (numreact == 2){
				iss >> a1 >> a2 >> a3 >> a4;
				for (int j = 0; j < size; j++){
					if (a1 == species[j]){
						react1 = j;
					}
					if (a3 == species[j]){
						react2 = j;
					}
				}
			}
			else{
				iss >> a1 >> a2 >> a3 >> a4 >> a5 >> a6;
				for (int j = 0; j < size; j++){
					if (a1 == species[j]){
						react1 = j;
					}
					if (a3 == species[j]){
						react2 = j;
					}
					if (a5 == species[j]){
						react3 = j;
					}
				}
			}


			if (numproduct == 1){
				iss >> p1;
				for (int j = 0; j < size; j++){
					if (p1 == species[j]){
						prod1 = j;
					}
				}
			}
			else if (numproduct == 2){
				iss >> p1 >> p2 >> p3;
				if (p1 == "M"){
					temp = 1;
					prod1 = size;
				}
				else if (p3 == "M"){
					temp = 1;
					prod2 = size;
				}
				for (int j = 0; j < size; j++){
					if (p1 == species[j]){
						prod1 = j;
					}
					if (p3 == species[j]){
						prod2 = j;
					}
				}
			}
			else{
				iss >> p1 >> p2 >> p3 >> p4 >> p5;
				if (p1 == "M"){
					temp = 1;
					prod1 = size;
				}
				else if (p3 == "M"){
					temp = 1;
					prod2 = size;
				}
				else if (p5 == "M"){
					temp = 1;
					prod3 = size;
				}
				for (int j = 0; j < size; j++){
					if (p1 == species[j]){
						prod1 = j;
					}
					if (p3 == species[j]){
						prod2 = j;
					}
					if (p5 == species[j]){
						prod3 = j;
					}
				}
			}
			std::cout << "Reading reactions : " << i << "\n";
			data.clear();
			data.resize(5);
			iss >> data[0] >> data[1] >> data[2] >> data[3] >> data[4];
			data[2] = double(0) - data[2];
			data[0] = double(0) - data[0];
			data[4] = double(0) - data[4];
			if (numreact == 1){
				data.push_back(react1);
				K[react1].push_back(data);
			}
			else if (numreact == 2){
				data.push_back(react1);
				data.push_back(react2);
				K[react1].push_back(data);
				K[react2].push_back(data);
			}
			else{
				data.push_back(react1);
				data.push_back(react2);
				data.push_back(react3);
				K[react1].push_back(data);
				K[react2].push_back(data);
				K[react3].push_back(data);
			}
			data[0] = double(0) - data[0];
			if (numproduct == 1){
				K[prod1].push_back(data);
			}
			else if (numproduct == 2){
				K[prod1].push_back(data);
				K[prod2].push_back(data);
			}
			else{
				K[prod1].push_back(data);
				K[prod2].push_back(data);
				K[prod3].push_back(data);
			}
		}
  myfile.close();
	std::ofstream myfile1;
	myfile1.open("Output.csv");
	//myfile1 << "Electron temperature : " << Te/double(11604) << " eV \n";
	myfile1 << "t (s)" << "," << "T_g (K)";
	for (int i = 0; i < size; i++){
		myfile1 << "," << species[i];
	}
	myfile1 << "\n";
	myfile1.close();

}

void initialize(){
	for (int i = 0; i < size; i++){
		sp[i] = 0.0;
	}
	//sp[21] = 1.0e13;
	sp[3] = 0.2095*n;
	sp[19] = 0.7809*n;
	sp[20] = 0.0093*n;
	sp[12] = 0.0003*n;
	sp[10] = eq*2.0*sp[3];
}

void solve_rxn(){
  double1D sp1, sp_temp;
  double1D k, k1, k2, k3, k4;
  sp1.resize(size);
  sp_temp.resize(size);
  k1.resize(size);
  k2.resize(size);
  k3.resize(size);
  k4.resize(size);
	initialize();
  write_file(-1);
  for (int i = 0; i < iter; i++){
    //sp1 = sp;
    calc_change(k1,sp);
    sp_temp = sp + k1*(dt/double(2));
    calc_change(k2,sp_temp);
		sp_temp = sp + k2*(dt/double(2));
    calc_change(k3,sp_temp);
		sp_temp = sp + k3*dt;
    calc_change(k4,sp_temp);
    k = k1/double(6) + k2/double(3) + k3/double(3) + k4/double(6);
    sp = sp + k*dt;
		calc_temp(k);
		std::cout << "Iteration : " << i << "\n";
    if (i%wf == (wf-1)){
      write_file(i);
    }
  }
}

void calc_change(double1D& t_k, double1D& t_s){
	double temp;
  for (unsigned int x = 0; x < K.size(); x++){
    t_k[x] = 0.0;
    for (unsigned int y = 0; y < K[x].size(); y++){
      temp = K[x][y][0]*pow(Tg,K[x][y][1])*exp(K[x][y][2]/(R*Tg))*pow(300.0/Te,K[x][y][3])*exp(K[x][y][4]/(R*Te));
			//std::cout << temp << "\t";
      for (unsigned int z = 5; z < K[x][y].size(); z++){
        temp *= t_s[int(K[x][y][z])];
      }
			//std::cout << temp << "\t";
      t_k[x] += temp;
    }
		//std::cout << "\n";
		//std::cout << t_k[x];
  }
}

void calc_temp(double1D& k){
	double dnE = 0.0;
	double n = 0.0;
	double nCp = 0.0;
	for (int i = 0; i < size; i++){
		dnE += (k[i]*E[i]);
		n += sp[i];
		nCp += (sp[i]*Cp[i]);
	}
	//std::cout << dnE*n/nCp*dt << "\n";
	Tg -= (dnE*n/nCp*dt*dt);
}

void write_file(int it){
  std::ofstream myfile;
  std::stringstream stream;
  myfile.open("Output.csv", std::ofstream::app);
	myfile << std::scientific << std::setprecision(5) << double((it+1)*dt) << "," << Tg;
  for (int i = 0; i < size; i++){
		myfile << "," << std::scientific << std::setprecision(5) << sp[i];
  }
  myfile << "\n";
  myfile.close();
}

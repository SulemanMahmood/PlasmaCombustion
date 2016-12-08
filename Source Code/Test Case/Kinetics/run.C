#include <cmath>
#include <sstream>
#include <fstream>
#include "utility.h"

std::string species_info = "methane_air_plasma_without_3rdbody.txt";
double Te = 100.0*11604.0 //in K
double Tg = 300.0 // in K
double end_time = 1.0e-6 // in s
double dt = 1.0e-9; // in s
int iter = int(end_time/dt);
double R = 8.314; // Gas constat in J/mol.K
double Av = 6.022e-23; // Avogadro's number

double3D K;
int size;
string1D species;
double1D sp;
double1D E;
double1D Cp;
int wf = 10; // Frequency of writing output to file

void read_file();
void solve_rxn();
void calc_change(double1D, double1D);
void write_file();
void calc_temp();

int main(){
  read_file();
  solve_rxn();
  return 0;
}

void read_file(){
  ifstream myfile;
	std::string line;
	std::istringstream iss;
	int temp, numreact, numproduct;
	string1D secspecies;
	double1D sec, data;
	data.resize(4);
	std::string a1, a2, a3, a4, a5, a6, p1, p2 , p3, p4, p5, asname;
	double as;
  myfile.open(species_info);
	getline(myfile,line);
	while(!myfile.eof()){
		getline(myfile,line);
		iss(line);
		iss >> size;
		getline(myfile,line);
		getline(myfile,line);
		for (int i = 0; i < size; i++){
			getline(myfile,line);
			iss(line);
			iss >> temp >> species[i] >> Cp[i] >> E[i];
		}
		getline(myfile,line);
		getline(myfile,line);
		getline(myfile,line);
		iss(line);
		iss >> temp;
		for (int i = 0; i < temp i++){
			getline(myfile,line);
			iss(line);
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
			if (numreact == 2){
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
			data.clear();
			data.resize(4);
			iss >> data[0] >> data[1] >> data[2] >> data[3];
			data[2] = 0.0 - data[2];
			data[0] = 0.0 - data[0];
			if (numreact == 1){
				data.push_back(react1);
				K[react1].push_back(double1D());
				K[react1][K[react1].size()-1].push_back(data);
			}
			else if (numreact == 2){
				data.push_back(react1);
				data.push_back(react2);
				K[react1].push_back(double1D());
				K[react1][K[react1].size()-1].push_back(data);
				K[react2].push_back(double1D());
				K[react2][K[react2].size()-1].push_back(data);
			}
			else{
				data.push_back(react1);
				data.push_back(react2);
				data.push_back(react3);
				K[react1].push_back(double1D());
				K[react1][K[react1].size()-1].push_back(data);
				K[react2].push_back(double1D());
				K[react2][K[react2].size()-1].push_back(data);
				K[react3].push_back(double1D());
				K[react3][K[react3].size()-1].push_back(data);
			}
			data[0] = 0.0 - data[0];
			if (numproduct == 1){
				K[prod1].push_back(double1D());
				K[prod1][K[prod1].size()-1].push_back(data);
			}
			else if (numproduct == 2){
				K[prod1].push_back(double1D());
				K[prod1][K[prod1].size()-1].push_back(data);
				K[prod2].push_back(double1D());
				K[prod2][K[prod2].size()-1].push_back(data);
			}
			else{
				K[prod1].push_back(double1D());
				K[prod1][K[prod1].size()-1].push_back(data);
				K[prod2].push_back(double1D());
				K[prod2][K[prod2].size()-1].push_back(data);
				K[prod3].push_back(double1D());
				K[prod3][K[prod3].size()-1].push_back(data);
			}
			/*if (temp == 1){
				iss >> numsec;
				secspecies.clear();
				sec.clear();
				for (int j = 0; j < numsec; j++){
					iss >> asname >> as;
					secspecies.push_back(asname);
					sec.push_back(as);
				}
			}*/
		}
	}
  myfile.close();
	std::ofstream myfile1;
	myfile1.open("Output.txt");
	myfile1 << "Electron temperature : " << Te/11604.0 << " eV \n";
	myfile1 << "Time (s)" << "\t" << "Gas Temperature (K)";
	for (int i = 0; i < size; i++){
		myfile1 << "\t" << species[i];
	}
	myfile1 << "\n";
	myfile1.close();
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
  write_file();
  for (int i = 0; i < iter; i++){
    //sp1 = sp;
    calc_change(k1,sp);
    sp_temp = sp + k1*(dt/2.0);
    calc_change(k2,sp_temp);
		sp_temp = sp + k2*(dt/2.0);
    calc_change(k3,sp_temp);
		sp_temp = sp + k3*dt;
    calc_change(k4,sp_temp);
    k = k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;
    sp += k*dt;
		calc_temp(k);
    if (i%wf == (wf-1)){
      write_file(i);
    }
  }
}

void calc_change(double1D t_k, double1D t_s){
	double temp;
  for (int x = 0; x < K.size(); x++){
    t_k[x] = 0.0;
    for (int y = 0; y < K[i].size(); y++){
      temp = K[x][y][0]*pow(300.0/Tg,K[x][y][1])*exp(K[x][y][2]/(R*Tg))*pow(300.0/Te,K[x][y][3]);
      for (int z = 4; z < K[x][y].size(); z++){
        temp *= t_s[K[x][y][z]];
      }
      t_k[x] += temp;
    }
  }
}

void calc_temp(double1D k){
	double dnE = 0.0;
	double n = 0.0;
	double nCp = 0.0;
	for (int i = 0; i < size; i++){
		dnE += (k[i]*E[i]);
		n += sp[i];
		nCp += (sp[i]*Cp[i]);
	}
	Tg += (dnE*n/nCp*dt*dt);
}

void write_file(int it = -1){
  std::ofstream myfile;
  std::stringstream stream;
  myfile.open("Output.txt", std::ofstream::app);
  stream << std::scientific << std::setprecision(1) << double((it+1)*dt);
  myfile << stream.str();
	stream << std::scientific << std::setprecision(1) << Tg;
  myfile << "\t" << stream.str();
  for (int i = 0; i < size; i++){
		stream << std::scientific << std::setprecision(1) << sp[i];
    myfile << "\t" << stream.str();
  }
  myfile << "\n";
  myfile.close();
}

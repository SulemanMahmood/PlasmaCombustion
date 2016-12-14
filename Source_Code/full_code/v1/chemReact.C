#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "utility.h"

//std::string species_info = ;
double Te = 5.0*11604.0; //in K
double Tg = 300.0; // in K
double end_time = 1.0e-8; // in s
double dt = 1.0e-15; // in s
int iter = int(end_time/dt);
double R = 8.314; // Gas constat in J/mol.K
double Av = 6.022e23; // Avogadro's number
double n = 2.5e19; // Number density of air
double eq = 0.1; // Equivalence Ratio

double1D H_f; // Enthalpy of rxn;
double2D K; // rxn rate constant for each rxn
int2D rs; // reactant species for each rxn
double1D adv; // advancement for each reaction
//double1D d_sp; // change in species concentration
int2D p_rxn; // participating rxn for each species
double2D r_p; // reactants or products
double2D add_info; // Third body efficiencies and pressure dependence
int2D tb_sp; // Third body
int size; // Number of species
int rxn_size; // Species rxns
string1D species; // species name
double1D sp; // species concentration
double1D Cp; // Specific heat
int wf = 10000; // Frequency of writing output to file

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
    sp[21] = 1.0e13;
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
    double k_o, k_f, k_inf, Pr, tb_mp, conc;
    for (int i = 0; i < rxn_size; i++){
        k_o = K[i][0]*pow(Tg,K[i][1])*exp(K[i][2]/(R*Tg))*pow(300.0/Te,K[i][3])*exp(K[i][4]/(R*Te));
        tb_mp = 0.0;
        for (unsigned int j = 0; j < tb_sp[i].size(); j++){
            tb_mp += add_info[i][j+4]*t_s[tb_sp[i][j]];
        }
        if (add_info[i][0] == double(1)){
            k_inf = add_info[i][1]*pow(Tg,add_info[i][2])*exp(add_info[i][3]/(R*Tg));
            Pr = k_o*tb_mp/k_inf;
            k_f = k_inf * Pr/(1.0 + Pr);
        }
        else{
            k_f = k_o;
        }
        conc = 1.0;
        for (unsigned int j = 0; j < rs[i].size(); j++){
            conc *= t_s[rs[i][j]];
        }
        if (tb_sp[i].size() == 0){
            tb_mp = 1.0;
        }
        adv[i] = k_f*conc*tb_mp;
    }
    for (int i = 0; i < size; i++){
        t_k[i] = 0.0;
        for (unsigned int j = 0; j < p_rxn[i].size(); j++){
            t_k[i] += (r_p[i][j]*adv[p_rxn[i][j]]);
        }
    }
}

void calc_temp(double1D& k){
    double dnH = 0.0;
    double n_total = 0.0;
    double Cp_mix = 0.0;
    for (int i = 0; i < size; i++){
        dnH += (k[i]*H_f[i]);
        n_total += sp[i];
        Cp_mix += (sp[i]*Cp[i]);
    }
    Cp_mix /= n_total;
    //std::cout << dnH/Cp_mix*dt*dt << "\n";
    Tg += (dnH/Cp_mix*dt*dt);
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

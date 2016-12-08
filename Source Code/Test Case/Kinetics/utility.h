#ifndef UTILITY_H_
#define UTILITY_H_

#include <vector>
#include <string>

typedef std::vector<double> double1D
typedef std::vector<std::vector<double>> double2D
typedef std::vector<std::vector<std::vector<double>>> double3D
typedef std::vector<std::string> string1D

double1D& operator*(const double1D &lhs, const double &rhs){
	double1D temp;
	temp.resize(lhs.size());
	for (int i = 0; i < lhs.size(); i++){
		temp[i] = lhs[i]*rhs;
	}
	return *temp;
}

double1D& operator/(const double1D &lhs, const double &rhs){
	double1D temp;
	temp.resize(lhs.size());
	for (int i = 0; i < lhs.size(); i++){
		temp[i] = lhs[i]/rhs;
	}
	return *temp;
}

#endif

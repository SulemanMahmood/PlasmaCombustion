#ifndef UTILITY_H
#define UTILITY_H

#include "pup_stl.h"
#include <vector>

struct flow{
  double r, u, v, w, E;
	double1D Y;

	void pup(PUP::er &p){
		p|r;
		p|u;
		p|v;
		p|w;
		p|E;
	}

	flow& operator=(const flow& k) {
		for (int i = 0; i < Y.size(); i++){
	    this->Y[i] = k.Y[i];
	  }
		this->r = k.r;
		this->u = k.u;
		this->v = k.v;
		this->w = k.w;
		this->E = k.E;
		return *this;
	}

	flow& operator+=(flow& rhs) {
		for (int i = 0; i < Y.size(); i++){
	    lhs->Y[i] += k.Y[i];
	  }
		this->r += rhs.r;
		this->u += rhs.u;
		this->v += rhs.v;
		this->w += rhs.w;
		this->E += rhs.E;
		return *this;
	}

	flow& operator-=(const flow& rhs) {
		this->r -= rhs.r;
		this->u -= rhs.u;
		this->v -= rhs.v;
		this->w -= rhs.w;
		this->E -= rhs.E;
		return *this;
	}

	flow& operator*=(const flow& rhs) {
		this->r *= rhs.r;
		this->u *= rhs.u;
		this->v *= rhs.v;
		this->w *= rhs.w;
		this->E *= rhs.E;
		return *this;
	}

	flow& operator/=(const flow& rhs) {
		this->r /= rhs.r;
		this->u /= rhs.u;
		this->v /= rhs.v;
		this->w /= rhs.w;
		this->E /= rhs.E;
		return *this;
	}

};

flow operator+(const flow& lhs, const flow& rhs) {
  flow temp;
  temp.r = lhs.r + rhs.r;
  temp.u = lhs.u + rhs.u;
  temp.v = lhs.v + rhs.v;
  temp.w = lhs.w + rhs.w;
  temp.E = lhs.E + rhs.E;
  return temp;
}

flow operator+(const flow& lhs, const double& k) {
  flow temp;
  temp.r = lhs.r + k;
  temp.u = lhs.u + k;
  temp.v = lhs.v + k;
  temp.w = lhs.w + k;
  temp.E = lhs.E + k;
  return temp;
}

flow operator+(const double& k, const flow& rhs) {
  flow temp;
  temp.r = rhs.r + k;
  temp.u = rhs.u + k;
  temp.v = rhs.v + k;
  temp.w = rhs.w + k;
  temp.E = rhs.E + k;
  return temp;
}

flow operator-(const flow& lhs, const flow& rhs) {
  flow temp;
  temp.r = lhs.r - rhs.r;
  temp.u = lhs.u - rhs.u;
  temp.v = lhs.v - rhs.v;
  temp.w = lhs.w - rhs.w;
  temp.E = lhs.E - rhs.E;
  return temp;
}

flow operator-(const flow& lhs, const double& k) {
  flow temp;
  temp.r = lhs.r - k;
  temp.u = lhs.u - k;
  temp.v = lhs.v - k;
  temp.w = lhs.w - k;
  temp.E = lhs.E - k;
  return temp;
}

flow operator-(const double& k, const flow& rhs) {
  flow temp;
  temp.r = k - rhs.r;
  temp.u = k - rhs.u;
  temp.v = k - rhs.v;
  temp.w = k - rhs.w;
  temp.E = k - rhs.E;
  return temp;
}

flow operator*(const flow& lhs, const flow& rhs) {
  flow temp;
  temp.r = lhs.r * rhs.r;
  temp.u = lhs.u * rhs.u;
  temp.v = lhs.v * rhs.v;
  temp.w = lhs.w * rhs.w;
  temp.E = lhs.E * rhs.E;
  return temp;
}

flow operator*(const flow& lhs, const double& k) {
  flow temp;
  temp.r = lhs.r * k;
  temp.u = lhs.u * k;
  temp.v = lhs.v * k;
  temp.w = lhs.w * k;
  temp.E = lhs.E * k;
  return temp;
}

flow operator*(const double& k, const flow& rhs) {
  flow temp;
  temp.r = k * rhs.r;
  temp.u = k * rhs.u;
  temp.v = k * rhs.v;
  temp.w = k * rhs.w;
  temp.E = k * rhs.E;
  return temp;
}

flow operator/(const flow& lhs, const flow& rhs) {
  flow temp;
  temp.r = lhs.r / rhs.r;
  temp.u = lhs.u / rhs.u;
  temp.v = lhs.v / rhs.v;
  temp.w = lhs.w / rhs.w;
  temp.E = lhs.E / rhs.E;
  return temp;
}

flow operator/(const flow& lhs, const double& k) {
  flow temp;
  temp.r = lhs.r / k;
  temp.u = lhs.u / k;
  temp.v = lhs.v / k;
  temp.w = lhs.w / k;
  temp.E = lhs.E / k;
  return temp;
}

flow operator/(const double& k, const flow& rhs) {
  flow temp;
  temp.r = k / rhs.r;
  temp.u = k / rhs.u;
  temp.v = k / rhs.v;
  temp.w = k / rhs.w;
  temp.E = k / rhs.E;
  return temp;
}

typedef std::vector<std::vector<double>> double2D;
typedef std::vector<std::vector<flow>> flow2D;
typedef std::vector<std::vector<std::vector<double>>> double3D;
typedef std::vector<std::vector<std::vector<flow>>> flow3D;
typedef std::vector<std::vector<std::vector<std::vector<flow>>>> flow4D;

template<typename T>
void copy3D(T &lhs, T rhs, int dim){
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      for(int k = 0; k < dim; k++){
        lhs[i][j][k] = rhs[i][j][k];
      }
    }
  }
}

template<typename T>
void copy2D(T &lhs, T rhs, int dim){
  for(int i = 0; i < dim; i++){
    for(int j = 0; j < dim; j++){
      lhs[i][j] = rhs[i][j];
    }
  }
}


// following are from the utility file originally for chemical reaction part.
typedef std::vector<std::vector<int> > int2D;
typedef std::vector<double> double1D;
typedef std::vector<std::string> string1D;

double1D operator*(const double1D &lhs, const double &rhs){
    double1D temp;
    temp.resize(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++){
        temp[i] = lhs[i]*rhs;
    }
    return temp;
}

double1D operator/(const double1D &lhs, const double &rhs){
    double1D temp;
    temp.resize(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++){
        temp[i] = lhs[i]/rhs;
    }
    return temp;
}

double1D operator+(const double1D &lhs, const double1D &rhs){
    double1D temp;
    temp.resize(lhs.size());
    for (unsigned int i = 0; i < lhs.size(); i++){
        temp[i] = lhs[i] + rhs[i];
    }
    return temp;
}

#endif

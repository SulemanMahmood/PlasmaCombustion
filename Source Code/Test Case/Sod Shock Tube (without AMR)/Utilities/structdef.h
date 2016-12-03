#ifndef STRUCTDEF_H_
#define STRUCTDEF_H_

#include <vector>

struct flow{
  double r, u, v, w, E;
  std::vector<double> Y;
};

/*flow& operator=(flow& lhs, const double& k) {
  for (int i = 0; i < lhs.Y.size(); i++){
    lhs.Y[i] = k;
  }
  lhs.r = k;
  lhs.u = k;
  lhs.v = k;
  lhs.w = k;
  lhs.E = k;
  return *lhs;
}*/

flow& operator+=(flow& lhs, const flow& rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    lhs.Y[i] += rhs.Y[i];
  }
  lhs->r += rhs->r;
  lhs.u += rhs.u;
  lhs.v += rhs.v;
  lhs.w += rhs.w;
  lhs.E += rhs.E;
  return *lhs;
}

flow& operator-=(flow& lhs, const flow& rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    lhs.Y[i] -= rhs.Y[i];
  }
  lhs.r -= rhs.r;
  lhs.u -= rhs.u;
  lhs.v -= rhs.v;
  lhs.w -= rhs.w;
  lhs.E -= rhs.E;
  return *lhs;
}

flow& operator*=(flow& lhs, const double rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    lhs.Y[i] *= rhs.Y[i];
  }
  lhs.r *= rhs;
  lhs.u *= rhs;
  lhs.v *= rhs;
  lhs.w *= rhs;
  lhs.E *= rhs;
  return *lhs;
}

flow& operator/=(flow& lhs, const double rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    lhs.Y[i] /= rhs.Y[i];
  }
  lhs.r /= rhs;
  lhs.u /= rhs;
  lhs.v /= rhs;
  lhs.w /= rhs;
  lhs.E /= rhs;
  return *lhs;
}

flow& operator+(const flow& lhs, const flow& rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    this.Y[i] = lhs.Y[i] + rhs.Y[i];
  }
  this.r = lhs.r + rhs.r;
  this.u = lhs.u + rhs.u;
  this.v = lhs.v + rhs.v;
  this.w = lhs.w + rhs.w;
  this.E = lhs.E + rhs.E;
  return *this;
}

flow& operator+(const flow& lhs, const double& k) {
  for (int i = 0; i < lhs.Y.size(); i++){
    this.Y[i] = lhs.Y[i] + k;
  }
  this.r = lhs.r + k;
  this.u = lhs.u + k;
  this.v = lhs.v + k;
  this.w = lhs.w + k;
  this.E = lhs.E + k;
  return *this;
}

flow& operator+(const double& k, const flow& rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    this.Y[i] = k + rhs.Y[i];
  }
  this.r = rhs.r + k;
  this.u = rhs.u + k;
  this.v = rhs.v + k;
  this.w = rhs.w + k;
  this.E = rhs.E + k;
  return *this;
}

flow& operator-(const flow& lhs, const flow& rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    this.Y[i] = lhs.Y[i] - rhs.Y[i];
  }
  this.r = lhs.r - rhs.r;
  this.u = lhs.u - rhs.u;
  this.v = lhs.v - rhs.v;
  this.w = lhs.w - rhs.w;
  this.E = lhs.E - rhs.E;
  return *this;
}

flow& operator-(const flow& lhs, const double& k) {
  for (int i = 0; i < lhs.Y.size(); i++){
    this.Y[i] = lhs.Y[i] - k;
  }
  this.r = lhs.r - k;
  this.u = lhs.u - k;
  this.v = lhs.v - k;
  this.w = lhs.w - k;
  this.E = lhs.E - k;
  return *this;
}

flow& operator-(const double& k, const flow& rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    this.Y[i] = k - rhs.Y[i];
  }
  this.r = k - rhs.r;
  this.u = k - rhs.u;
  this.v = k - rhs.v;
  this.w = k - rhs.w;
  this.E = k - rhs.E;
  return *this;
}

flow& operator*(const flow& lhs, const flow& rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    this.Y[i] = lhs.Y[i] * rhs.Y[i];
  }
  this.r = lhs.r * rhs.r;
  this.u = lhs.u * rhs.u;
  this.v = lhs.v * rhs.v;
  this.w = lhs.w * rhs.w;
  this.E = lhs.E * rhs.E;
  return *this;
}

flow& operator*(const flow& lhs, const double& k) {
  for (int i = 0; i < lhs.Y.size(); i++){
    this.Y[i] = lhs.Y[i] * k;
  }
  this.r = lhs.r * k;
  this.u = lhs.u * k;
  this.v = lhs.v * k;
  this.w = lhs.w * k;
  this.E = lhs.E * k;
  return *this;
}

flow& operator*(const double& k, const flow& rhs) {
  for (int i = 0; i < rhs.Y.size(); i++){
    this.Y[i] = k * rhs.Y[i];
  }
  this.r = k * rhs.r;
  this.u = k * rhs.u;
  this.v = k * rhs.v;
  this.w = k * rhs.w;
  this.E = k * rhs.E;
  return *this;
}

flow& operator/(const flow& lhs, const flow& rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    this.Y[i] = lhs.Y[i] / rhs.Y[i];
  }
  this.r = lhs.r / rhs.r;
  this.u = lhs.u / rhs.u;
  this.v = lhs.v / rhs.v;
  this.w = lhs.w / rhs.w;
  this.E = lhs.E / rhs.E;
  return *this;
}

flow& operator/(const flow& lhs, const double& k) {
  for (int i = 0; i < lhs.Y.size(); i++){
    this.Y[i] = lhs.Y[i] / k;
  }
  this.r = lhs.r / k;
  this.u = lhs.u / k;
  this.v = lhs.v / k;
  this.w = lhs.w / k;
  this.E = lhs.E / k;
  return *this;
}

flow& operator/(const double& k, const flow& rhs) {
  for (int i = 0; i < rhs.Y.size(); i++){
    this.Y[i] = k / rhs.Y[i];
  }
  this.r = k / rhs.r;
  this.u = k / rhs.u;
  this.v = k / rhs.v;
  this.w = k / rhs.w;
  this.E = k / rhs.E;
  return *this;
}

#endif

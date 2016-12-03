#ifndef STRUCTDEF_H_
#define STRUCTDEF_H_

#include <vector>

struct flow{
  double r, u, v, w, E;
  std::vector<double> Y;
};

flow operator=(const flow& k) {
  flow temp;
  for (int i = 0; i < temp.Y.size(); i++){
    temp.Y[i] = k.Y[i];
  }
  temp.r = k.r;
  temp.u = k.u;
  temp.v = k.v;
  temp.w = k.w;
  temp.E = k.E;
  return temp;
}

flow operator+=(flow& lhs, const flow& rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    lhs.Y[i] += rhs.Y[i];
  }
  lhs.r += rhs.r;
  lhs.u += rhs.u;
  lhs.v += rhs.v;
  lhs.w += rhs.w;
  lhs.E += rhs.E;
  return lhs;
}

flow operator-=(flow& lhs, const flow& rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    lhs.Y[i] -= rhs.Y[i];
  }
  lhs.r -= rhs.r;
  lhs.u -= rhs.u;
  lhs.v -= rhs.v;
  lhs.w -= rhs.w;
  lhs.E -= rhs.E;
  return lhs;
}

flow operator*=(flow& lhs, const double rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    lhs.Y[i] *= rhs;
  }
  lhs.r *= rhs;
  lhs.u *= rhs;
  lhs.v *= rhs;
  lhs.w *= rhs;
  lhs.E *= rhs;
  return lhs;
}

flow operator/=(flow& lhs, const double rhs) {
  for (int i = 0; i < lhs.Y.size(); i++){
    lhs.Y[i] /= rhs;
  }
  lhs.r /= rhs;
  lhs.u /= rhs;
  lhs.v /= rhs;
  lhs.w /= rhs;
  lhs.E /= rhs;
  return lhs;
}

flow operator+(const flow& lhs, const flow& rhs) {
  flow temp;
  for (int i = 0; i < lhs.Y.size(); i++){
    temp.Y[i] = lhs.Y[i] + rhs.Y[i];
  }
  temp.r = lhs.r + rhs.r;
  temp.u = lhs.u + rhs.u;
  temp.v = lhs.v + rhs.v;
  temp.w = lhs.w + rhs.w;
  temp.E = lhs.E + rhs.E;
  return temp;
}

flow operator+(const flow& lhs, const double& k) {
  flow temp;
  for (int i = 0; i < lhs.Y.size(); i++){
    temp.Y[i] = lhs.Y[i] + k;
  }
  temp.r = lhs.r + k;
  temp.u = lhs.u + k;
  temp.v = lhs.v + k;
  temp.w = lhs.w + k;
  temp.E = lhs.E + k;
  return temp;
}

flow operator+(const double& k, const flow& rhs) {
  flow temp;
  for (int i = 0; i < rhs.Y.size(); i++){
    temp.Y[i] = k + rhs.Y[i];
  }
  temp.r = rhs.r + k;
  temp.u = rhs.u + k;
  temp.v = rhs.v + k;
  temp.w = rhs.w + k;
  temp.E = rhs.E + k;
  return temp;
}

flow operator-(const flow& lhs, const flow& rhs) {
  flow temp;
  for (int i = 0; i < lhs.Y.size(); i++){
    temp.Y[i] = lhs.Y[i] - rhs.Y[i];
  }
  temp.r = lhs.r - rhs.r;
  temp.u = lhs.u - rhs.u;
  temp.v = lhs.v - rhs.v;
  temp.w = lhs.w - rhs.w;
  temp.E = lhs.E - rhs.E;
  return temp;
}

flow operator-(const flow& lhs, const double& k) {
  flow temp;
  for (int i = 0; i < lhs.Y.size(); i++){
    temp.Y[i] = lhs.Y[i] - k;
  }
  temp.r = lhs.r - k;
  temp.u = lhs.u - k;
  temp.v = lhs.v - k;
  temp.w = lhs.w - k;
  temp.E = lhs.E - k;
  return temp;
}

flow operator-(const double& k, const flow& rhs) {
  flow temp;
  for (int i = 0; i < rhs.Y.size(); i++){
    temp.Y[i] = k - rhs.Y[i];
  }
  temp.r = k - rhs.r;
  temp.u = k - rhs.u;
  temp.v = k - rhs.v;
  temp.w = k - rhs.w;
  temp.E = k - rhs.E;
  return temp;
}

flow operator*(const flow& lhs, const flow& rhs) {
  flow temp;
  for (int i = 0; i < lhs.Y.size(); i++){
    temp.Y[i] = lhs.Y[i] * rhs.Y[i];
  }
  temp.r = lhs.r * rhs.r;
  temp.u = lhs.u * rhs.u;
  temp.v = lhs.v * rhs.v;
  temp.w = lhs.w * rhs.w;
  temp.E = lhs.E * rhs.E;
  return temp;
}

flow operator*(const flow& lhs, const double& k) {
  flow temp;
  for (int i = 0; i < lhs.Y.size(); i++){
    temp.Y[i] = lhs.Y[i] * k;
  }
  temp.r = lhs.r * k;
  temp.u = lhs.u * k;
  temp.v = lhs.v * k;
  temp.w = lhs.w * k;
  temp.E = lhs.E * k;
  return temp;
}

flow operator*(const double& k, const flow& rhs) {
  flow temp;
  for (int i = 0; i < rhs.Y.size(); i++){
    temp.Y[i] = k * rhs.Y[i];
  }
  temp.r = k * rhs.r;
  temp.u = k * rhs.u;
  temp.v = k * rhs.v;
  temp.w = k * rhs.w;
  temp.E = k * rhs.E;
  return temp;
}

flow operator/(const flow& lhs, const flow& rhs) {
  flow temp;
  for (int i = 0; i < lhs.Y.size(); i++){
    temp.Y[i] = lhs.Y[i] / rhs.Y[i];
  }
  temp.r = lhs.r / rhs.r;
  temp.u = lhs.u / rhs.u;
  temp.v = lhs.v / rhs.v;
  temp.w = lhs.w / rhs.w;
  temp.E = lhs.E / rhs.E;
  return temp;
}

flow operator/(const flow& lhs, const double& k) {
  flow temp;
  for (int i = 0; i < lhs.Y.size(); i++){
    temp.Y[i] = lhs.Y[i] / k;
  }
  temp.r = lhs.r / k;
  temp.u = lhs.u / k;
  temp.v = lhs.v / k;
  temp.w = lhs.w / k;
  temp.E = lhs.E / k;
  return temp;
}

flow operator/(const double& k, const flow& rhs) {
  flow temp;
  for (int i = 0; i < rhs.Y.size(); i++){
    temp.Y[i] = k / rhs.Y[i];
  }
  temp.r = k / rhs.r;
  temp.u = k / rhs.u;
  temp.v = k / rhs.v;
  temp.w = k / rhs.w;
  temp.E = k / rhs.E;
  return temp;
}

#endif

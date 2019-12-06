#include <math.h>

double FISCOout (
  double f0,
  double m0,
  double p0,
  double r_H,
  double x,
  double f1,
  double m1,
  double p1,
  double f2,
  double m2,
  double p2)
{
  return(((0.128e3 * pow(-0.1e1 + x, 0.6e1) * (0.8e1 * x - 0.16e2) * m0 * m0 + (0.1024e4 * pow(-0.1e1 + x, 0.7e1) * (x + 0.1e1) * m1 + 0.128e3 * pow(-0.1e1 + x, 0.7e1) * (x * x - 0.1e1) * (x + 0.1e1) * m2) * m0 - 0.256e3 * pow(-0.1e1 + x, 0.8e1) * m1 * m1 * pow(x + 0.1e1, 0.2e1)) * f1 + (0.128e3 * pow(-0.1e1 + x, 0.6e1) * (0.4e1 * x * x - 0.4e1) * m0 * m0 + 0.128e3 * pow(-0.1e1 + x, 0.7e1) * (-0.1e1 * x * x + 0.1e1) * (x + 0.1e1) * m1 * m0) * f2) * f0 + (0.128e3 * pow(-0.1e1 + x, 0.6e1) * (-0.8e1 * x * x + 0.8e1) * m0 * m0 + 0.128e3 * pow(-0.1e1 + x, 0.7e1) * (0.2e1 * x * x - 0.2e1) * (x + 0.1e1) * m1 * m0) * f1 * f1);
}
#include <math.h>

double dFISCOdxout (
  double f0,
  double m0,
  double p0,
  double r_H,
  double x,
  double f1,
  double m1,
  double p1,
  double f2,
  double m2,
  double p2)
{
  return(((0.1280e4 * pow(-0.1e1 + x, 0.5e1) * (-0.1040000000e2 + 0.5600000000e1 * x) * m0 * m0 + (0.1280e4 * pow(-0.1e1 + x, 0.6e1) * (0.6400000000e1 * x + 0.4800000000e1) * m1 + 0.1280e4 * pow(-0.1e1 + x, 0.6e1) * (pow(x, 0.3e1) + 0.6000000000e0 * x * x - 0.1e1 * x - 0.6000000000e0) * m2) * m0 - 0.2560e4 * pow(-0.1e1 + x, 0.7e1) * m1 * m1 * (x + 0.6000000000e0) * (x + 0.1e1)) * f1 + (0.1280e4 * pow(-0.1e1 + x, 0.5e1) * (0.3200000000e1 * x * x - 0.8000000000e0 * x - 0.2400000000e1) * m0 * m0 + 0.1280e4 * pow(-0.1e1 + x, 0.6e1) * (-0.1e1 * pow(x, 0.3e1) - 0.6000000000e0 * x * x + x + 0.6000000000e0) * m1 * m0) * f2) * f0 + (0.1280e4 * pow(-0.1e1 + x, 0.5e1) * (-0.6400000000e1 * x * x + 0.1600000000e1 * x + 0.4800000000e1) * m0 * m0 + 0.1280e4 * pow(-0.1e1 + x, 0.6e1) * (0.2e1 * pow(x, 0.3e1) + 0.1200000000e1 * x * x - 0.2e1 * x - 0.1200000000e1) * m1 * m0) * f1 * f1);
}
double FLRout (
  double f0,
  double m0,
  double p0,
  double r_H,
  double x,
  double f1,
  double m1,
  double p1,
  double f2,
  double m2,
  double p2)
{
  return(((0.8e1 - 0.8e1 * x) * m0 + 0.2e1 * (-0.1e1 + x) * (x * x - 0.1e1) * m1) * f0 + 0.2e1 * (-0.1e1 + x) * (-0.1e1 * x * x + 0.1e1) * m0 * f1);
}
double dFLRdxout (
  double f0,
  double m0,
  double p0,
  double r_H,
  double x,
  double f1,
  double m1,
  double p1,
  double f2,
  double m2,
  double p2)
{
  return((-0.8e1 * m0 + (0.6e1 * x * x - 0.4e1 * x - 0.2e1) * m1) * f0 - 0.6e1 * m0 * (-0.1e1 + x) * f1 * (x + 0.3333333333e0));
}

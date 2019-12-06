double Gpsiout (
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
  return((0.4e1 * x + 0.4e1) / (x + 0.1e1) * p2 + 0.2e1 * p1 / f0 * f1 + 0.2e1 * p1 / m0 * m1 + 0.8e1 * p1 / (x + 0.1e1));
}
#include <math.h>

double GBout (
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
  return((((0.32e2 * pow(-0.1e1 + x, 0.2e1) * (-0.16e2 * x - 0.8e1) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m1 - 0.128e3 * pow(r_H, -0.2e1) / (x + 0.1e1) * pow(-0.1e1 + x, 0.3e1) * m2) * pow(m0, -0.2e1) + (0.32e2 * (0.3e1 * x + 0.7e1) * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * pow(r_H, -0.2e1) * m1 * m1 + 0.32e2 * pow(-0.1e1 + x, 0.4e1) * m2 * pow(r_H, -0.2e1) * m1) * pow(m0, -0.3e1) - 0.40e2 * pow(m1, 0.3e1) * pow(-0.1e1 + x, 0.4e1) * pow(r_H, -0.2e1) * pow(m0, -0.4e1)) * f1 + (0.32e2 * pow(-0.1e1 + x, 0.2e1) * (-0.4e1 * x * x + 0.4e1) * m1 * pow(x + 0.1e1, -0.2e1) * pow(r_H, -0.2e1) * pow(m0, -0.2e1) + 0.32e2 * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * (0.5000000000e0 * x * x - 0.5000000000e0) * m1 * m1 * pow(r_H, -0.2e1) * pow(m0, -0.3e1)) * f2) / f0 + (0.32e2 * pow(-0.1e1 + x, 0.2e1) * (0.2e1 * x * x - 0.2e1) * m1 * pow(x + 0.1e1, -0.2e1) * pow(r_H, -0.2e1) * pow(m0, -0.2e1) + 0.32e2 * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * (-0.2500000000e0 * x * x + 0.2500000000e0) * m1 * m1 * pow(r_H, -0.2e1) * pow(m0, -0.3e1)) * f1 * f1 * pow(f0, -0.2e1));
}
#include <math.h>

double GBcalcout (
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
  return((((0.32e2 * pow(-0.1e1 + x, 0.6e1) * (-0.16e2 * x - 0.8e1) * pow(x + 0.1e1, -0.2e1) * pow(r_H, -0.4e1) * m1 - 0.128e3 * pow(-0.1e1 + x, 0.7e1) * m2 / (x + 0.1e1) * pow(r_H, -0.4e1)) * pow(m0, -0.3e1) + (0.32e2 * pow(-0.1e1 + x, 0.7e1) / (x + 0.1e1) * (0.3e1 * x + 0.7e1) * pow(r_H, -0.4e1) * m1 * m1 + 0.32e2 * pow(-0.1e1 + x, 0.8e1) * m2 * pow(r_H, -0.4e1) * m1) * pow(m0, -0.4e1) - 0.40e2 * pow(-0.1e1 + x, 0.8e1) * pow(m1, 0.3e1) * pow(r_H, -0.4e1) * pow(m0, -0.5e1)) * f1 + (0.32e2 * pow(-0.1e1 + x, 0.6e1) * (-0.4e1 * x * x + 0.4e1) * m1 * pow(x + 0.1e1, -0.2e1) * pow(r_H, -0.4e1) * pow(m0, -0.3e1) + 0.32e2 * pow(-0.1e1 + x, 0.7e1) / (x + 0.1e1) * (0.5000000000e0 * x * x - 0.5000000000e0) * m1 * m1 * pow(r_H, -0.4e1) * pow(m0, -0.4e1)) * f2) / f0 + (0.32e2 * pow(-0.1e1 + x, 0.6e1) * (0.2e1 * x * x - 0.2e1) * m1 * pow(x + 0.1e1, -0.2e1) * pow(r_H, -0.4e1) * pow(m0, -0.3e1) + 0.32e2 * pow(-0.1e1 + x, 0.7e1) / (x + 0.1e1) * (-0.2500000000e0 * x * x + 0.2500000000e0) * m1 * m1 * pow(r_H, -0.4e1) * pow(m0, -0.4e1)) * f1 * f1 * pow(f0, -0.2e1));
}
#include <math.h>

double R2calcout (
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
  return((0.256e3 * pow(-0.1e1 + x, 0.8e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m1 * m1 + 0.64e2 * pow(-0.1e1 + x, 0.8e1) * (0.4e1 * x + 0.4e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m2 * m1 + 0.64e2 * pow(-0.1e1 + x, 0.8e1) * pow(r_H, -0.4e1) * m2 * m2) * pow(m0, -0.4e1) + (-0.192e3 * pow(-0.1e1 + x, 0.8e1) / (x + 0.1e1) * pow(r_H, -0.4e1) * pow(m1, 0.3e1) - 0.96e2 * pow(-0.1e1 + x, 0.8e1) * m2 * pow(r_H, -0.4e1) * m1 * m1) * pow(m0, -0.5e1) + 0.36e2 * pow(m0, -0.6e1) * pow(-0.1e1 + x, 0.8e1) * pow(r_H, -0.4e1) * pow(m1, 0.4e1) + (((0.256e3 * pow(-0.1e1 + x, 0.8e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m1 + 0.128e3 * pow(-0.1e1 + x, 0.8e1) / (x + 0.1e1) * pow(r_H, -0.4e1) * m2) * pow(m0, -0.3e1) + (-0.32e2 * pow(-0.1e1 + x, 0.8e1) / (x + 0.1e1) * pow(r_H, -0.4e1) * m1 * m1 + 0.32e2 * pow(-0.1e1 + x, 0.8e1) * m2 * pow(r_H, -0.4e1) * m1) * pow(m0, -0.4e1) - 0.24e2 * pow(-0.1e1 + x, 0.8e1) * pow(m1, 0.3e1) * pow(r_H, -0.4e1) * pow(m0, -0.5e1)) * f1 + ((0.64e2 * pow(-0.1e1 + x, 0.8e1) * (0.2e1 * x + 0.2e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m1 + 0.32e2 * pow(-0.1e1 + x, 0.8e1) * (0.2e1 * x + 0.2e1) / (x + 0.1e1) * pow(r_H, -0.4e1) * m2) * pow(m0, -0.3e1) - 0.24e2 * pow(-0.1e1 + x, 0.8e1) * (0.2e1 * x + 0.2e1) / (x + 0.1e1) * m1 * m1 * pow(r_H, -0.4e1) * pow(m0, -0.4e1)) * f2) / f0 + ((0.64e2 * pow(-0.1e1 + x, 0.8e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * pow(m0, -0.2e1) + (0.64e2 * pow(-0.1e1 + x, 0.8e1) * (-0.5000000000e0 * x - 0.5000000000e0) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m1 - 0.32e2 * pow(-0.1e1 + x, 0.8e1) * m2 * pow(r_H, -0.4e1)) * pow(m0, -0.3e1) + 0.28e2 * pow(-0.1e1 + x, 0.8e1) * m1 * m1 * pow(r_H, -0.4e1) * pow(m0, -0.4e1)) * f1 * f1 + (0.32e2 * pow(-0.1e1 + x, 0.8e1) * (0.2e1 * x + 0.2e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * pow(m0, -0.2e1) + 0.8e1 * pow(-0.1e1 + x, 0.8e1) * (0.2e1 * x + 0.2e1) * m1 / (x + 0.1e1) * pow(r_H, -0.4e1) * pow(m0, -0.3e1)) * f2 * f1 + 0.4e1 * pow(-0.1e1 + x, 0.8e1) * pow(m0, -0.2e1) * pow(0.2e1 * x + 0.2e1, 0.2e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * f2 * f2) * pow(f0, -0.2e1) + ((-0.32e2 * pow(-0.1e1 + x, 0.8e1) / (x + 0.1e1) * pow(r_H, -0.4e1) * pow(m0, -0.2e1) - 0.8e1 * pow(-0.1e1 + x, 0.8e1) * m1 * pow(r_H, -0.4e1) * pow(m0, -0.3e1)) * pow(f1, 0.3e1) - 0.8e1 * pow(-0.1e1 + x, 0.8e1) * pow(m0, -0.2e1) / (x + 0.1e1) * (0.2e1 * x + 0.2e1) * pow(r_H, -0.4e1) * f2 * f1 * f1) * pow(f0, -0.3e1) + 0.4e1 * pow(m0, -0.2e1) * pow(f0, -0.4e1) * pow(-0.1e1 + x, 0.8e1) * pow(r_H, -0.4e1) * pow(f1, 0.4e1));
}
#include <math.h>

double Rab2calcout (
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
  return((0.24e2 * (0.4e1 * x * x - 0.5333333333e1 * x + 0.5333333333e1) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m1 * m1 + 0.24e2 * (0.4e1 * pow(x, 0.3e1) - 0.2666666667e1 * x * x - 0.4e1 * x + 0.2666666667e1) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m2 * m1 + 0.24e2 * (pow(x, 0.4e1) - 0.2e1 * x * x + 0.1e1) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m2 * m2) * pow(m0, -0.4e1) + (-0.40e2 * pow(-0.1e1 + x, 0.7e1) * (0.2e1 * x - 0.8000000000e0) / (x + 0.1e1) * pow(r_H, -0.4e1) * pow(m1, 0.3e1) - 0.40e2 * pow(-0.1e1 + x, 0.7e1) * (x * x - 0.1e1) / (x + 0.1e1) * pow(r_H, -0.4e1) * m2 * m1 * m1) * pow(m0, -0.5e1) + 0.18e2 * pow(m0, -0.6e1) * pow(-0.1e1 + x, 0.8e1) * pow(r_H, -0.4e1) * pow(m1, 0.4e1) + (((0.4e1 * (0.16e2 * x * x + 0.32e2) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m1 + 0.4e1 * (0.8e1 * pow(x, 0.3e1) - 0.8e1 * x) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m2) * pow(m0, -0.3e1) - 0.16e2 * m1 * m1 * pow(-0.1e1 + x, 0.7e1) / (x + 0.1e1) * (0.2e1 * x + 0.3e1) * pow(r_H, -0.4e1) * pow(m0, -0.4e1) + 0.4e1 * pow(-0.1e1 + x, 0.8e1) * pow(m1, 0.3e1) * pow(r_H, -0.4e1) * pow(m0, -0.5e1)) * f1 + ((0.4e1 * (0.8e1 * pow(x, 0.3e1) - 0.8e1 * x) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m1 + 0.4e1 * (0.4e1 * pow(x, 0.4e1) - 0.8e1 * x * x + 0.4e1) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m2) * pow(m0, -0.3e1) - 0.16e2 * pow(-0.1e1 + x, 0.7e1) * (x * x - 0.1e1) * m1 * m1 / (x + 0.1e1) * pow(r_H, -0.4e1) * pow(m0, -0.4e1)) * f2) / f0 + ((-0.8e1 * (-0.4e1 * x * x - 0.8e1) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * pow(m0, -0.2e1) + (-0.8e1 * pow(-0.1e1 + x, 0.7e1) / (x + 0.1e1) * (0.2e1 * x + 0.4e1) * pow(r_H, -0.4e1) * m1 - 0.8e1 * pow(-0.1e1 + x, 0.7e1) * (x * x - 0.1e1) / (x + 0.1e1) * pow(r_H, -0.4e1) * m2) * pow(m0, -0.3e1) + 0.12e2 * pow(-0.1e1 + x, 0.8e1) * m1 * m1 * pow(r_H, -0.4e1) * pow(m0, -0.4e1)) * f1 * f1 - 0.8e1 * (-0.4e1 * pow(x, 0.3e1) + 0.4e1 * x) * pow(m0, -0.2e1) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * f2 * f1 - 0.8e1 * (-0.1e1 * pow(x, 0.4e1) + 0.2e1 * x * x - 0.1e1) * pow(m0, -0.2e1) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * f2 * f2) * pow(f0, -0.2e1) + (-0.16e2 * pow(m0, -0.2e1) * pow(-0.1e1 + x, 0.7e1) / (x + 0.1e1) * x * pow(r_H, -0.4e1) * pow(f1, 0.3e1) - 0.8e1 * pow(m0, -0.2e1) * pow(-0.1e1 + x, 0.7e1) / (x + 0.1e1) * (x * x - 0.1e1) * pow(r_H, -0.4e1) * f2 * f1 * f1) * pow(f0, -0.3e1) + 0.2e1 * pow(m0, -0.2e1) * pow(f0, -0.4e1) * pow(-0.1e1 + x, 0.8e1) * pow(r_H, -0.4e1) * pow(f1, 0.4e1));
}
#include <math.h>

double Rabcd2calcout (
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
  return((0.32e2 * (0.4e1 * x * x + 0.8e1) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m1 * m1 + 0.32e2 * (0.4e1 * pow(x, 0.3e1) - 0.4e1 * x) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m2 * m1 + 0.32e2 * (pow(x, 0.4e1) - 0.2e1 * x * x + 0.1e1) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * m2 * m2) * pow(m0, -0.4e1) + (-0.64e2 * pow(-0.1e1 + x, 0.7e1) / (x + 0.1e1) * (0.2e1 * x + 0.1e1) * pow(r_H, -0.4e1) * pow(m1, 0.3e1) - 0.64e2 * pow(-0.1e1 + x, 0.7e1) * (x * x - 0.1e1) / (x + 0.1e1) * pow(r_H, -0.4e1) * m2 * m1 * m1) * pow(m0, -0.5e1) + 0.36e2 * pow(m0, -0.6e1) * pow(-0.1e1 + x, 0.8e1) * pow(r_H, -0.4e1) * pow(m1, 0.4e1) + ((0.12e2 * (0.5333333333e1 * x * x + 0.1066666667e2 * x + 0.16e2) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * pow(m0, -0.2e1) - 0.16e2 * pow(-0.1e1 + x, 0.7e1) * m1 / (x + 0.1e1) * (0.2e1 * x + 0.6e1) * pow(r_H, -0.4e1) * pow(m0, -0.3e1) + 0.12e2 * pow(-0.1e1 + x, 0.8e1) * m1 * m1 * pow(r_H, -0.4e1) * pow(m0, -0.4e1)) * f1 * f1 + (0.12e2 * (0.5333333333e1 * pow(x, 0.3e1) + 0.5333333333e1 * x * x - 0.5333333333e1 * x - 0.5333333333e1) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * pow(m0, -0.2e1) - 0.16e2 * pow(-0.1e1 + x, 0.7e1) * m1 / (x + 0.1e1) * (x * x - 0.1e1) * pow(r_H, -0.4e1) * pow(m0, -0.3e1)) * f2 * f1 + 0.12e2 * pow(m0, -0.2e1) * (0.1333333333e1 * pow(x, 0.4e1) - 0.2666666667e1 * x * x + 0.1333333333e1) * pow(-0.1e1 + x, 0.6e1) * pow(r_H, -0.4e1) * pow(x + 0.1e1, -0.2e1) * f2 * f2) * pow(f0, -0.2e1) + ((-0.32e2 * pow(-0.1e1 + x, 0.7e1) * pow(r_H, -0.4e1) * pow(m0, -0.2e1) + 0.8e1 * pow(-0.1e1 + x, 0.8e1) * m1 * pow(r_H, -0.4e1) * pow(m0, -0.3e1)) * pow(f1, 0.3e1) + 0.8e1 * pow(m0, -0.2e1) * pow(-0.1e1 + x, 0.7e1) * (-0.2e1 * x + 0.2e1) * pow(r_H, -0.4e1) * f2 * f1 * f1) * pow(f0, -0.3e1) + 0.4e1 * pow(m0, -0.2e1) * pow(f0, -0.4e1) * pow(-0.1e1 + x, 0.8e1) * pow(r_H, -0.4e1) * pow(f1, 0.4e1));
}
double dGpsidf2out (
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
  return(0.0e0);
}
double dGpsidf1out (
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
  return(0.2e1 * p1 / f0);
}
#include <math.h>

double dGpsidf0out (
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
  return(-0.2e1 * p1 * pow(f0, -0.2e1) * f1);
}
#include <math.h>

double dGBdf2out (
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
  return((-0.128e3 * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * pow(r_H, -0.2e1) * m1 * pow(m0, -0.2e1) - 0.16e2 * pow(-0.1e1 + x, 0.3e1) * (-0.1e1 * x * x + 0.1e1) / (x + 0.1e1) * pow(r_H, -0.2e1) * m1 * m1 * pow(m0, -0.3e1)) / f0);
}
#include <math.h>

double dGBdf1out (
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
  return(((0.32e2 * pow(-0.1e1 + x, 0.2e1) * (-0.16e2 * x - 0.8e1) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m1 + 0.32e2 * pow(-0.1e1 + x, 0.2e1) * (-0.4e1 * x * x + 0.4e1) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m2) * pow(m0, -0.2e1) + (0.32e2 * (0.3e1 * x + 0.7e1) * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * pow(r_H, -0.2e1) * m1 * m1 + 0.32e2 * pow(-0.1e1 + x, 0.4e1) * m2 * pow(r_H, -0.2e1) * m1) * pow(m0, -0.3e1) - 0.40e2 * pow(m1, 0.3e1) * pow(-0.1e1 + x, 0.4e1) * pow(r_H, -0.2e1) * pow(m0, -0.4e1)) / f0 + (0.32e2 * pow(-0.1e1 + x, 0.2e1) * (0.4e1 * x * x - 0.4e1) * m1 * pow(x + 0.1e1, -0.2e1) * pow(r_H, -0.2e1) * pow(m0, -0.2e1) + 0.32e2 * pow(-0.1e1 + x, 0.3e1) * m1 * m1 / (x + 0.1e1) * (-0.5000000000e0 * x * x + 0.5000000000e0) * pow(r_H, -0.2e1) * pow(m0, -0.3e1)) * f1 * pow(f0, -0.2e1));
}
#include <math.h>

double dGBdf0out (
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
  return((((-0.32e2 * pow(-0.1e1 + x, 0.2e1) * (-0.16e2 * x - 0.8e1) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m1 + 0.128e3 * pow(r_H, -0.2e1) / (x + 0.1e1) * pow(-0.1e1 + x, 0.3e1) * m2) * pow(m0, -0.2e1) + (-0.32e2 * (0.3e1 * x + 0.7e1) * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * pow(r_H, -0.2e1) * m1 * m1 - 0.32e2 * pow(-0.1e1 + x, 0.4e1) * m2 * pow(r_H, -0.2e1) * m1) * pow(m0, -0.3e1) + 0.40e2 * pow(m1, 0.3e1) * pow(-0.1e1 + x, 0.4e1) * pow(r_H, -0.2e1) * pow(m0, -0.4e1)) * f1 + (-0.32e2 * pow(-0.1e1 + x, 0.2e1) * (-0.4e1 * x * x + 0.4e1) * m1 * pow(x + 0.1e1, -0.2e1) * pow(r_H, -0.2e1) * pow(m0, -0.2e1) - 0.32e2 * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * (0.5000000000e0 * x * x - 0.5000000000e0) * m1 * m1 * pow(r_H, -0.2e1) * pow(m0, -0.3e1)) * f2) * pow(f0, -0.2e1) + (-0.32e2 * pow(-0.1e1 + x, 0.2e1) * (0.4e1 * x * x - 0.4e1) * m1 * pow(x + 0.1e1, -0.2e1) * pow(r_H, -0.2e1) * pow(m0, -0.2e1) - 0.32e2 * pow(-0.1e1 + x, 0.3e1) * m1 * m1 / (x + 0.1e1) * (-0.5000000000e0 * x * x + 0.5000000000e0) * pow(r_H, -0.2e1) * pow(m0, -0.3e1)) * f1 * f1 * pow(f0, -0.3e1));
}
double dGpsidm2out (
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
  return(0.0e0);
}
double dGpsidm1out (
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
  return(0.2e1 * p1 / m0);
}
#include <math.h>

double dGpsidm0out (
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
  return(-0.2e1 * p1 * pow(m0, -0.2e1) * m1);
}
#include <math.h>

double dGBdm2out (
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
  return((-0.128e3 * pow(r_H, -0.2e1) / (x + 0.1e1) * pow(-0.1e1 + x, 0.3e1) * pow(m0, -0.2e1) - 0.32e2 * pow(-0.1e1 + x, 0.3e1) * (-0.1e1 * x * x + 0.1e1) * pow(r_H, -0.2e1) / (x + 0.1e1) * m1 * pow(m0, -0.3e1)) * f1 / f0);
}
#include <math.h>

double dGBdm1out (
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
  return(((0.32e2 * pow(-0.1e1 + x, 0.2e1) * (-0.16e2 * x - 0.8e1) * pow(x + 0.1e1, -0.2e1) * pow(r_H, -0.2e1) * pow(m0, -0.2e1) + (0.32e2 * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * (0.6e1 * x + 0.14e2) * pow(r_H, -0.2e1) * m1 + 0.32e2 * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * (x * x - 0.1e1) * pow(r_H, -0.2e1) * m2) * pow(m0, -0.3e1) - 0.120e3 * m1 * m1 * pow(-0.1e1 + x, 0.4e1) * pow(r_H, -0.2e1) * pow(m0, -0.4e1)) * f1 + (0.32e2 * pow(-0.1e1 + x, 0.2e1) * (-0.4e1 * x * x + 0.4e1) * pow(x + 0.1e1, -0.2e1) * pow(r_H, -0.2e1) * pow(m0, -0.2e1) + 0.32e2 * pow(-0.1e1 + x, 0.4e1) * pow(r_H, -0.2e1) * m1 * pow(m0, -0.3e1)) * f2) / f0 + (0.32e2 * pow(-0.1e1 + x, 0.2e1) * (0.2e1 * x * x - 0.2e1) * pow(x + 0.1e1, -0.2e1) * pow(r_H, -0.2e1) * pow(m0, -0.2e1) - 0.16e2 * pow(-0.1e1 + x, 0.4e1) * pow(r_H, -0.2e1) * m1 * pow(m0, -0.3e1)) * f1 * f1 * pow(f0, -0.2e1));
}
#include <math.h>

double dGBdm0out (
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
  return((((-0.96e2 * pow(-0.1e1 + x, 0.2e1) * (-0.1066666667e2 * x - 0.5333333333e1) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m1 + 0.256e3 * pow(r_H, -0.2e1) / (x + 0.1e1) * pow(-0.1e1 + x, 0.3e1) * m2) * pow(m0, -0.3e1) + (-0.96e2 * (0.3e1 * x + 0.7e1) * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * pow(r_H, -0.2e1) * m1 * m1 - 0.96e2 * pow(-0.1e1 + x, 0.4e1) * m2 * pow(r_H, -0.2e1) * m1) * pow(m0, -0.4e1) + 0.160e3 * pow(m1, 0.3e1) * pow(-0.1e1 + x, 0.4e1) * pow(r_H, -0.2e1) * pow(m0, -0.5e1)) * f1 + (-0.96e2 * pow(-0.1e1 + x, 0.2e1) * (-0.2666666667e1 * x * x + 0.2666666667e1) * m1 * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * pow(m0, -0.3e1) - 0.96e2 * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * (0.5000000000e0 * x * x - 0.5000000000e0) * m1 * m1 * pow(r_H, -0.2e1) * pow(m0, -0.4e1)) * f2) / f0 + (-0.96e2 * pow(-0.1e1 + x, 0.2e1) * (0.1333333333e1 * x * x - 0.1333333333e1) * m1 * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * pow(m0, -0.3e1) - 0.96e2 * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * (-0.2500000000e0 * x * x + 0.2500000000e0) * m1 * m1 * pow(r_H, -0.2e1) * pow(m0, -0.4e1)) * f1 * f1 * pow(f0, -0.2e1));
}
double dGpsidp2out (
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
  return(0.4e1);
}
double dGpsidp1out (
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
  return(0.8e1 / (x + 0.1e1) + (0.2e1 * x + 0.2e1) / (x + 0.1e1) * m1 / m0 + 0.2e1 / f0 * f1);
}
double dGpsidp0out (
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
  return(0.0e0);
}
double dGBdp2out (
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
  return(0.0e0);
}
double dGBdp1out (
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
  return(0.0e0);
}
double dGBdp0out (
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
  return(0.0e0);
}

double Gtout (
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
  return(0.8e1 / (x + 0.1e1) * m1 + (0.4e1 * x + 0.4e1) / (x + 0.1e1) * m2 - 0.3e1 * m1 * m1 / m0);
}
double Ttout (
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
  return(-0.2e1 * p1 * p1 * m0);
}
#include <math.h>

double Ktout (
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
  return(((0.256e3 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * pow(r_H, -0.2e1) * p2 + 0.256e3 * pow(-0.1e1 + x, 0.2e1) * exp(p0) * (x * x - 0.1e1) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p1 * p1 + 0.256e3 * pow(-0.1e1 + x, 0.2e1) * exp(p0) * (0.4e1 * x + 0.2e1) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p1) * m1 + 0.256e3 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * p1 * pow(r_H, -0.2e1) * m2) / m0 + ((-0.32e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * (x * x - 0.1e1) * pow(r_H, -0.2e1) * p2 - 0.32e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * (x * x - 0.1e1) * pow(r_H, -0.2e1) * p1 * p1 - 0.32e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) * (0.6e1 * x + 0.14e2) / (x + 0.1e1) * pow(r_H, -0.2e1) * p1) * m1 * m1 - 0.64e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * (x * x - 0.1e1) * p1 * pow(r_H, -0.2e1) * m1 * m2) * pow(m0, -0.2e1) - 0.32e2 * pow(-0.1e1 + x, 0.2e1) * exp(p0) * (-0.2500000000e1 * pow(x, 0.4e1) + 0.5e1 * x * x - 0.2500000000e1) * p1 * pow(m0, -0.3e1) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * pow(m1, 0.3e1));
}
double dGtdf2out (
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
double dGtdf1out (
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
double dGtdf0out (
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
double dTtdf2out (
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
double dTtdf1out (
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
double dTtdf0out (
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
double dKtdf2out (
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
double dKtdf1out (
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
double dKtdf0out (
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
double dGtdm2out (
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
double dGtdm1out (
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
  return(0.8e1 / (x + 0.1e1) + (-0.6e1 * x - 0.6e1) * m1 / (x + 0.1e1) / m0);
}
#include <math.h>

double dGtdm0out (
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
  return(0.3e1 * m1 * m1 * pow(m0, -0.2e1));
}
double dTtdm2out (
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
double dTtdm1out (
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
double dTtdm0out (
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
  return(-0.2e1 * p1 * p1);
}
#include <math.h>

double dKtdm2out (
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
  return(0.256e3 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) / m0 * p1 * pow(r_H, -0.2e1) + 0.64e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) * p1 * (-0.1e1 * x * x + 0.1e1) / (x + 0.1e1) * pow(r_H, -0.2e1) * m1 * pow(m0, -0.2e1));
}
#include <math.h>

double dKtdm1out (
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
  return((-0.64e2 * (-0.4e1 * x * x + 0.4e1) * pow(-0.1e1 + x, 0.2e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p1 * p1 - 0.64e2 * (-0.16e2 * x - 0.8e1) * pow(-0.1e1 + x, 0.2e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p1 - 0.64e2 * (-0.4e1 * x * x + 0.4e1) * pow(-0.1e1 + x, 0.2e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p2) / m0 + ((-0.64e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * (x * x - 0.1e1) * pow(r_H, -0.2e1) * p1 * p1 - 0.64e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) * (0.6e1 * x + 0.14e2) / (x + 0.1e1) * pow(r_H, -0.2e1) * p1 - 0.64e2 * pow(-0.1e1 + x, 0.4e1) * p2 * exp(p0) * pow(r_H, -0.2e1)) * m1 - 0.64e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * (x * x - 0.1e1) * p1 * pow(r_H, -0.2e1) * m2) * pow(m0, -0.2e1) + 0.240e3 * m1 * m1 * p1 * pow(-0.1e1 + x, 0.4e1) * exp(p0) * pow(r_H, -0.2e1) * pow(m0, -0.3e1));
}
#include <math.h>

double dKtdm0out (
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
  return(((0.64e2 * (-0.4e1 * x * x + 0.4e1) * pow(-0.1e1 + x, 0.2e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p1 * p1 + 0.64e2 * (-0.16e2 * x - 0.8e1) * pow(-0.1e1 + x, 0.2e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p1 + 0.64e2 * (-0.4e1 * x * x + 0.4e1) * pow(-0.1e1 + x, 0.2e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p2) * m1 + 0.64e2 * (-0.4e1 * x * x + 0.4e1) * pow(-0.1e1 + x, 0.2e1) * exp(p0) * p1 * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m2) * pow(m0, -0.2e1) + ((0.64e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * (x * x - 0.1e1) * pow(r_H, -0.2e1) * p1 * p1 + 0.64e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) * (0.6e1 * x + 0.14e2) / (x + 0.1e1) * pow(r_H, -0.2e1) * p1 + 0.64e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * (x * x - 0.1e1) * pow(r_H, -0.2e1) * p2) * m1 * m1 + 0.64e2 * exp(p0) * pow(-0.1e1 + x, 0.3e1) * (0.2e1 * x * x - 0.2e1) * p1 / (x + 0.1e1) * pow(r_H, -0.2e1) * m2 * m1) * pow(m0, -0.3e1) - 0.240e3 * pow(-0.1e1 + x, 0.4e1) * exp(p0) * pow(m1, 0.3e1) * p1 * pow(r_H, -0.2e1) * pow(m0, -0.4e1));
}
double dGtdp2out (
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
double dGtdp1out (
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
double dGtdp0out (
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
double dTtdp2out (
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
double dTtdp1out (
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
  return(-0.4e1 * p1 * m0);
}
double dTtdp0out (
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
#include <math.h>

double dKtdp2out (
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
  return(0.256e3 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) / m0 * pow(r_H, -0.2e1) * m1 + 0.32e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) * (-0.1e1 * x * x + 0.1e1) / (x + 0.1e1) * pow(r_H, -0.2e1) * m1 * m1 * pow(m0, -0.2e1));
}
#include <math.h>

double dKtdp1out (
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
  return(((-0.64e2 * pow(-0.1e1 + x, 0.2e1) * (-0.8e1 * x * x + 0.8e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p1 - 0.64e2 * (-0.16e2 * x - 0.8e1) * pow(-0.1e1 + x, 0.2e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1)) * m1 - 0.64e2 * (-0.4e1 * x * x + 0.4e1) * pow(-0.1e1 + x, 0.2e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m2) / m0 + ((-0.64e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * (x * x - 0.1e1) * p1 * pow(r_H, -0.2e1) - 0.64e2 * pow(-0.1e1 + x, 0.3e1) / (x + 0.1e1) * (0.3e1 * x + 0.7e1) * exp(p0) * pow(r_H, -0.2e1)) * m1 * m1 - 0.64e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * (x * x - 0.1e1) * pow(r_H, -0.2e1) * m1 * m2) * pow(m0, -0.2e1) + 0.80e2 * pow(-0.1e1 + x, 0.4e1) * pow(m1, 0.3e1) * exp(p0) * pow(r_H, -0.2e1) * pow(m0, -0.3e1));
}
#include <math.h>

double dKtdp0out (
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
  return(((-0.32e2 * pow(-0.1e1 + x, 0.2e1) * (-0.8e1 * x * x + 0.8e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p1 * p1 - 0.32e2 * pow(-0.1e1 + x, 0.2e1) * exp(p0) * (-0.32e2 * x - 0.16e2) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p1 - 0.32e2 * pow(-0.1e1 + x, 0.2e1) * (-0.8e1 * x * x + 0.8e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p2) * m1 - 0.32e2 * pow(-0.1e1 + x, 0.2e1) * (-0.8e1 * x * x + 0.8e1) * exp(p0) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * p1 * m2) / m0 + ((-0.32e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * (x * x - 0.1e1) * pow(r_H, -0.2e1) * p2 - 0.32e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) / (x + 0.1e1) * (x * x - 0.1e1) * pow(r_H, -0.2e1) * p1 * p1 - 0.32e2 * pow(-0.1e1 + x, 0.3e1) * exp(p0) * (0.6e1 * x + 0.14e2) / (x + 0.1e1) * pow(r_H, -0.2e1) * p1) * m1 * m1 - 0.32e2 * exp(p0) * pow(-0.1e1 + x, 0.3e1) * (0.2e1 * x * x - 0.2e1) * p1 / (x + 0.1e1) * pow(r_H, -0.2e1) * m2 * m1) * pow(m0, -0.2e1) + 0.80e2 * pow(-0.1e1 + x, 0.4e1) * exp(p0) * pow(m1, 0.3e1) * p1 * pow(r_H, -0.2e1) * pow(m0, -0.3e1));
}

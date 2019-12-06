double Grout (
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
  return(-0.8e1 * m0 / (x * x - 0.1e1) * m1 + m1 * m1 + (-0.8e1 / (x * x - 0.1e1) * m0 * m0 + 0.2e1 * m1 * m0) * f1 / f0);
}
double Trout (
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
  return(0.2e1 * p1 * p1 * m0 * m0);
}
#include <math.h>

double Krout (
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
  return((-0.48e2 * pow(-0.1e1 + x, 0.2e1) * (pow(x, 0.4e1) - 0.2e1 * x * x + 0.1e1) * p1 * pow(r_H, -0.2e1) / m0 * pow(x + 0.1e1, -0.2e1) * m1 * m1 - 0.48e2 * (-0.8e1 * x * x + 0.8e1) * pow(-0.1e1 + x, 0.2e1) * p1 * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m1 - 0.512e3 * m0 * pow(-0.1e1 + x, 0.2e1) * p1 * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1)) * f1 / f0);
}
double dGrdf2out (
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
double dGrdf1out (
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
  return((-0.8e1 / (x * x - 0.1e1) * m0 * m0 + 0.2e1 * m1 * m0) / f0);
}
#include <math.h>

double dGrdf0out (
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
  return((0.8e1 / (x * x - 0.1e1) * m0 * m0 - 0.2e1 * m1 * m0) * f1 * pow(f0, -0.2e1));
}
double dTrdf2out (
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
double dTrdf1out (
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
double dTrdf0out (
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
double dKrdf2out (
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

double dKrdf1out (
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
  return((-0.48e2 * pow(-0.1e1 + x, 0.2e1) * (pow(x, 0.4e1) - 0.2e1 * x * x + 0.1e1) * p1 * pow(r_H, -0.2e1) / m0 * pow(x + 0.1e1, -0.2e1) * m1 * m1 - 0.48e2 * (-0.8e1 * x * x + 0.8e1) * pow(-0.1e1 + x, 0.2e1) * p1 * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m1 - 0.512e3 * m0 * pow(-0.1e1 + x, 0.2e1) * p1 * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1)) / f0);
}
#include <math.h>

double dKrdf0out (
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
  return((0.512e3 * m0 * pow(-0.1e1 + x, 0.2e1) * p1 * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) + 0.48e2 * (-0.8e1 * x * x + 0.8e1) * pow(-0.1e1 + x, 0.2e1) * p1 * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m1 + 0.48e2 * pow(-0.1e1 + x, 0.2e1) * (pow(x, 0.4e1) - 0.2e1 * x * x + 0.1e1) * p1 * pow(r_H, -0.2e1) / m0 * pow(x + 0.1e1, -0.2e1) * m1 * m1) * f1 * pow(f0, -0.2e1));
}
double dGrdm2out (
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
double dGrdm1out (
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
  return(-0.8e1 * m0 / (x * x - 0.1e1) + (0.2e1 * x * x - 0.2e1) / (x * x - 0.1e1) * m1 + (0.2e1 * x * x - 0.2e1) * m0 / (x * x - 0.1e1) * f1 / f0);
}
double dGrdm0out (
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
  return(-0.8e1 * m1 / (x * x - 0.1e1) + (-0.16e2 * m0 / (x * x - 0.1e1) + (0.2e1 * x * x - 0.2e1) / (x * x - 0.1e1) * m1) * f1 / f0);
}
double dTrdm2out (
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
double dTrdm1out (
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
double dTrdm0out (
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
  return(0.4e1 * p1 * p1 * m0);
}
double dKrdm2out (
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

double dKrdm1out (
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
  return((0.384e3 * pow(r_H, -0.2e1) / (x + 0.1e1) * pow(-0.1e1 + x, 0.3e1) * p1 - 0.96e2 * pow(-0.1e1 + x, 0.3e1) * (x * x - 0.1e1) / (x + 0.1e1) * p1 * pow(r_H, -0.2e1) * m1 / m0) * f1 / f0);
}
#include <math.h>

double dKrdm0out (
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
  return((0.48e2 * pow(-0.1e1 + x, 0.2e1) * (pow(x, 0.4e1) - 0.2e1 * x * x + 0.1e1) * p1 * pow(r_H, -0.2e1) * pow(m0, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m1 * m1 - 0.512e3 * pow(-0.1e1 + x, 0.2e1) * p1 * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1)) * f1 / f0);
}
double dGrdp2out (
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
double dGrdp1out (
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
double dGrdp0out (
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
double dTrdp2out (
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
double dTrdp1out (
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
  return(0.4e1 * m0 * m0 * p1);
}
double dTrdp0out (
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
double dKrdp2out (
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

double dKrdp1out (
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
  return((-0.48e2 * pow(-0.1e1 + x, 0.2e1) * (pow(x, 0.4e1) - 0.2e1 * x * x + 0.1e1) * pow(r_H, -0.2e1) / m0 * pow(x + 0.1e1, -0.2e1) * m1 * m1 - 0.48e2 * (-0.8e1 * x * x + 0.8e1) * pow(-0.1e1 + x, 0.2e1) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1) * m1 - 0.512e3 * m0 * pow(-0.1e1 + x, 0.2e1) * pow(r_H, -0.2e1) * pow(x + 0.1e1, -0.2e1)) * f1 / f0);
}
double dKrdp0out (
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

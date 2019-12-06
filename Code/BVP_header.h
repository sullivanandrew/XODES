#ifndef HEADER_BVP
#define HEADER_BVP

/*
 *File List and Description
 *BVP.c - main
 *BVP_header.c - header
 *BVP_control.c - main loop called by BVP.c
 *BVP_IC.c - initial conditions
 *BVP_sys.c - system evaluation at current step
 *BVP_gridsize.c - calculate resized grid
 *BVP_out.c - output
 *
 *Linear System solver
 *BVP_LSsolver.c - linear solver
 *BVP_BICO.c - BICOnjugate gradient stabilized method: ~4N^2 flops
 *BVP_GMRES.c - Generalized Minimal RESidual method: ~2N^2 flops
 *BVP_LUdcmp.c - LU DeCoMPosition: ~1/2N^3 flops
 *
 *Executable Files
 *BVP.exe - main .exe
 *BVPbatch.sh - bash file
 *BVPcompile.sh - compile bash file
 */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//Main iteration loop
int BVP_control(int r, int nstep, int p, double tol, double alpha, double beta, double kappa, double r_H);
//Initial conditions
int BVP_IC(void *params, double x[], double Psi[]);
//Evaluation of system at current iteration
int BVP_sys(void *params, const int r, double x[], double Psi[], double Jac[], double bvec[], double Dvec[], double RGB[], double R2[], double Rab2[], double Rabcd2[]);
//Output
int BVP_out(void *params, double x[], double Psi[], double Jac[], double dx[], double b[], double D[], int it, double RGB[], double R2[], double Rab2[], double Rabcd2[]);
//Calculates and interpolates new grid size
int BVP_gridsize(void *params, int r, double x[], double Psi[], double dxo[], double dxn[], double **xp, double **pp);
//Interpolates system at new grid size
int BVP_interp(void *params, int r, double x[], double Psi[], int nnew, double xtemp[], double Psitemp[], int ntemp);
int BVP_physics(void *params, double x[], double Psi[], double bvec[], double Dvec[], double RGB[], double R2[], double Rab2[], double Rabcd2[]);

//Calculates stencil
int steninit(int s[], int *oneside, int *twoside, int k, int r, int N);
//Calculates Newton polynomial coefficients
int BVP_NPcalc(void *params, int r, int rn, int I, int s[], double x[], double y[], int pn, double a_a[], double ax_a[], double axx_a[], double *ydk, double *dydk, double *d2ydk, double xk);

//Linear Solver
int BVP_LSsolver(void *params, double x[], double Psi[], double Jac[], double dx[], double b[], double err[]);
//LU DeCoMPosition
int BVP_LUdcmp(void *params, double A[], double b[], double x[], double *Det);
//BICOnjugate gradient stabilized method
int BVP_BICO(void *params, double A[], double b[], double x[], double err[]);
//Generalized Minimal RESidual method
int BVP_GMRES(void *params, double A[], double b[], double x[], double err[]);

double norm(void *params, double v[]);
double Jac_norm(void *params, double v[]);

int stenBC(int s[], int *oneside, int k, int r, int N);
int BVP_NPcalcBC(void *params, int r, int rn, int Ip, int I, int s[], double x[], double y[], int pn, double a_a[], double ax_a[], double axx_a[], double *ydk, double *dydk, double *d2ydk);

//Initializing struct
struct param_type
{
	double alphaparam;
	double betaparam;
	double gammaparam;
	double kappaparam;
	int nparam;
	int mparam;
	int pparam;
	int Nparam;
	int rparam;
	double r_Hparam;
	double tolparam;
	double LStolparam;
	int imaxparam;
};

double Gtout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double Ttout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double Ktout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double dGtdf2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGtdf1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGtdf0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGtdm2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGtdm1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGtdm0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGtdp2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGtdp1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGtdp0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double dTtdf2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTtdf1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTtdf0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTtdm2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTtdm1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTtdm0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTtdp2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTtdp1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTtdp0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double dKtdf2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKtdf1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKtdf0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKtdm2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKtdm1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKtdm0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKtdp2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKtdp1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKtdp0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double Grout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double Trout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double Krout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double dGrdf2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGrdf1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGrdf0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGrdm2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGrdm1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGrdm0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGrdp2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGrdp1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGrdp0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double dTrdf2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTrdf1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTrdf0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTrdm2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTrdm1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTrdm0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTrdp2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTrdp1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dTrdp0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double dKrdf2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKrdf1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKrdf0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKrdm2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKrdm1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKrdm0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKrdp2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKrdp1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dKrdp0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double Gpsiout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double GBout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double GBcalcout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double R2calcout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double Rab2calcout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double Rabcd2calcout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double dGpsidf2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGpsidf1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGpsidf0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGpsidm2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGpsidm1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGpsidm0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGpsidp2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGpsidp1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGpsidp0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double dGBdf2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGBdf1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGBdf0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGBdm2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGBdm1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGBdm0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGBdp2out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGBdp1out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dGBdp0out (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double FISCOout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dFISCOdxout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);

double FLRout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);
double dFLRdxout (double f0, double m0, double p0, double r_H, double x, double f1, double m1, double p1, double f2, double m2, double p2);


#endif

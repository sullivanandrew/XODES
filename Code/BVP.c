/*
 *Andrew Sullivan
 *Montana State University
 *Newton-Rapson method for solving PDE BVP in 1D as coupled system
 *Procedure:
 *Notes:
 *precision - double
 */

#include "BVP_header.h"

int main(int argc, char *argv[])
{
	int i, j, k, l, status;
	//Size
	int r;
	int n;
	int p = 3;
	double alpha;
	double tol;
	double r_H;

	//External Input
	alpha = atof(argv[1])*1.0E-5;
	tol = pow(10.0,atof(argv[2]));
	n = atoi(argv[3]);
	r = atoi(argv[4]);
	r_H = atof(argv[5]);

	double kappa = 1.0;
	double beta = 1.0;
	//Print statements
	printf("\nRunning BVP with a coupling value alpha = %.5f, and tol = %8.1e.\n",alpha,tol);
	printf("r = %i. Number of points n = %i.\n",r,n);

	//Executable
	status = BVP_control(r, n, p, tol, alpha, beta, kappa, r_H);


	return 0;
}

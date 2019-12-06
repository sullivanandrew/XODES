#include "BVP_header.h"

int BVP_out(void *params, double x[], double Psi[], double Jac[], double dx[], double b[], double D[], int it, double RGB[], double R2[], double Rab2[], double Rabcd2[])
{
	struct param_type *my_ptr;
	my_ptr = params;
	const int n = my_ptr->nparam;
	const int p = my_ptr->pparam;
	const int N = my_ptr->Nparam;
	const double r_H = my_ptr->r_Hparam;
	const double alpha = my_ptr->alphaparam;
	const double beta = my_ptr->betaparam;
	const double gamma = my_ptr->gammaparam;
	const double kappa = my_ptr->kappaparam;
	const double tol = my_ptr->tolparam;
	int i,j,k,l, status;
	int pn;

	char poutname[256];
	FILE *poutfile;
	if(it == -3)
	{
		sprintf(poutname,"./Data/BVPout_sols/Sol_Funcs/sol_al%05i_tol%03i_it%02i_rH%03i.dat",(int) round(alpha*1.0E5), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2));

	}
	else
	{
		sprintf(poutname,"./Data/BVPout_sols/Sol_Funcs/sol_al%05i_tol%03i_it%02i_rH%03i.dat",(int) round(alpha*1.0E5), (int) round(-log10(tol)), it, (int) round(r_H*1.0E2));
	}
	poutfile = fopen(poutname,"w");
	if(poutfile == NULL)
	{//Open output file
		perror("Error opening file.");
	}//Open output file
	for (j = 0; j < n; ++j)
	{
		fprintf(poutfile,"%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e\n",x[j],Psi[j*p + 0], Psi[j*p + 1], Psi[j*p + 2], b[j*p +0], b[j*p +1], b[j*p +2], D[j*p +0], D[j*p +1], D[j*p +2], RGB[j], R2[j], Rab2[j], Rabcd2[j]);
	}
	fclose(poutfile);


	return 0;
}

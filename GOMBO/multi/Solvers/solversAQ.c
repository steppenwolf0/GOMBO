#include "problemDef.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../Amazigh/mathMatrix.h"
#include "../Amazigh/mathVector.h"
#include "../Gaussian/gaussian.h"
#ifdef OSisWindows
#include "./RV-GOMEA/RV-GOMEAWIN.h"
#else
#include "./RV-GOMEA/RV-GOMEA.h"
#endif
int maxEval ;

Individual* RVGOMEAUCB(double** Ainv, double** xtrain,double* param,double* ytrain,
	int n_param,int n_x, int Dimension, double kappa, int kernelType, int evals)
{
	maxEval = evals;
	type = kernelType;
	problem_Solving = 0;
	initializeproblemUCB( Ainv,  xtrain,  param,  ytrain, n_param, n_x, Dimension, kappa, kernelType);
	
	//-b -w -s -r -f 1 -1e2 1e2 0.35 10 25 0.9 1 0 1e-10 100 0.0  10000
	const char** argv[19];
	//argv[17]= { "","-s", "-r", "-f", "1", "-1e2", "1e2", "0.35", "10", "25", "0.9", "1", "0", "1e-10", "100", "0.0", "10" };
	argv[0] = " ";
	argv[1] = "-b";
	argv[2] = "-b";
	//argv[3] = "-s"; generational statistics
	argv[3] = "-b";
	argv[4] = "-b";
	argv[5] = "-f";
	argv[6] = "1";
	argv[7] = "0";
	argv[8] = "1";
	argv[9] = "0.35";
	argv[10] = "200"; //play with this one
	argv[11] = "5";
	argv[12] = "0.9";
	argv[13] = "1";
	argv[14] = "0";
	argv[15] = "-1e8"; //without -r this does not matter
	argv[16] = "100";
	argv[17] = "0.0";
	argv[18] = "30";

	//noError = noError && sscanf(argv[*index + 0], "%lf", &lower_user_range);
	//noError = noError && sscanf(argv[*index + 1], "%lf", &upper_user_range);
	//noError = noError && sscanf(argv[*index + 2], "%lf", &tau);
	//noError = noError && sscanf(argv[*index + 3], "%d", &base_population_size);
	//noError = noError && sscanf(argv[*index + 4], "%d", &maximum_number_of_populations);
	//noError = noError && sscanf(argv[*index + 5], "%lf", &distribution_multiplier_decrease);
	//noError = noError && sscanf(argv[*index + 6], "%lf", &st_dev_ratio_threshold);
	//noError = noError && sscanf(argv[*index + 7], "%lf", &maximum_number_of_evaluations);
	//noError = noError && sscanf(argv[*index + 8], "%lf", &vtr);
	//noError = noError && sscanf(argv[*index + 9], "%d", &maximum_no_improvement_stretch);
	//noError = noError && sscanf(argv[*index + 10], "%lf", &fitness_variance_tolerance);
	//noError = noError && sscanf(argv[*index + 11], "%lf", &maximum_number_of_seconds);
	
	int argc = 19;
	interpretCommandLine(argc, argv);
	
	Individual* ouput = run();

	closeProblem();
	
	return ouput;
}

Individual* RVGOMEAEI(double** Ainv, double** xtrain, double* param, double* ytrain, int n_param, int n_x, int Dimension, double bestSolution, int kernelType)
{
	type = kernelType;
	problem_Solving = 2;
	initializeproblemEI(Ainv, xtrain, param, ytrain, n_param, n_x, Dimension, bestSolution);

	
	//-b -w -s -r -f 1 -1e2 1e2 0.35 10 25 0.9 1 0 1e-10 100 0.0  10000
	const char** argv[19];
	//argv[17]= { "","-s", "-r", "-f", "1", "-1e2", "1e2", "0.35", "10", "25", "0.9", "1", "0", "1e-10", "100", "0.0", "10" };
	argv[0] = " ";
	argv[1] = "-b";
	argv[2] = "-b";
	//argv[3] = "-s"; generational statistics
	argv[3] = "-b";
	argv[4] = "-r";
	argv[5] = "-f";
	argv[6] = "1";
	argv[7] = "-1e2";
	argv[8] = "1e2";
	argv[9] = "0.35";
	argv[10] = "10";
	argv[11] = "25";
	argv[12] = "0.9";
	argv[13] = "1";
	argv[14] = "0";
	argv[15] = "-1e4";
	argv[16] = "100";
	argv[17] = "0.0";
	argv[18] = "10";

	

	/*printf("\n");
	printf("# Tau                     = %e\n", tau);
	printf("# Population size/normal  = %d\n", base_population_size);
	printf("# FOS element size        = %d\n", FOS_element_size);
	printf("# Max num of populations  = %d\n", maximum_number_of_populations);
	printf("# Dis. mult. decreaser    = %e\n", distribution_multiplier_decrease);
	printf("# St. dev. rat. threshold = %e\n", st_dev_ratio_threshold);
	printf("# Maximum numb. of eval.  = %lf\n", maximum_number_of_evaluations);
	printf("# Value to reach (vtr)    = %e\n", vtr);
	printf("# Max. no improv. stretch = %d\n", maximum_no_improvement_stretch);
	printf("# Fitness var. tolerance  = %e\n", fitness_variance_tolerance);
	printf("# Random seed             = %ld\n", random_seed);
	printf("#\n");*/

	int argc = 19;
	interpretCommandLine(argc, argv);

	Individual* ouput = run();

	closeProblem();

	return ouput;
}

/*
//Individual* RVGOMEA3A(double** AinvKexp, double** AinvKM52, double** xtrain, double* paramKexp, double* paramKM52, double* ytrain, int n_param, int n_x,
//	int Dimension, double bestSolution, double power)
//{
//
//	problem_Solving = 20;
//	 initializeproblem2A(AinvKexp, AinvKM52, xtrain, paramKexp,
//		paramKM52, ytrain, n_param, n_x, Dimension, bestSolution, power);
//
//	//-b -w -s -r -f 1 -1e2 1e2 0.35 10 25 0.9 1 0 1e-10 100 0.0  10000
//	const char** argv[19];
//	//argv[17]= { "","-s", "-r", "-f", "1", "-1e2", "1e2", "0.35", "10", "25", "0.9", "1", "0", "1e-10", "100", "0.0", "10" };
//	argv[0] = " ";
//	argv[1] = "-b";
//	argv[2] = "-b";
//	//argv[3] = "-s"; generational statistics
//	argv[3] = "-b";
//	argv[4] = "-r";
//	argv[5] = "-f";
//	argv[6] = "1";
//	argv[7] = "-1e2";
//	argv[8] = "1e2";
//	argv[9] = "0.35";
//	argv[10] = "10";
//	argv[11] = "25";
//	argv[12] = "0.9";
//	argv[13] = "1";
//	argv[14] = "0";
//	argv[15] = "-1000";
//	argv[16] = "100";
//	argv[17] = "0.0";
//	argv[18] = "10";
//
//
//
//	/*printf("\n");
//	printf("# Tau                     = %e\n", tau);
//	printf("# Population size/normal  = %d\n", base_population_size);
//	printf("# FOS element size        = %d\n", FOS_element_size);
//	printf("# Max num of populations  = %d\n", maximum_number_of_populations);
//	printf("# Dis. mult. decreaser    = %e\n", distribution_multiplier_decrease);
//	printf("# St. dev. rat. threshold = %e\n", st_dev_ratio_threshold);
//	printf("# Maximum numb. of eval.  = %lf\n", maximum_number_of_evaluations);
//	printf("# Value to reach (vtr)    = %e\n", vtr);
//	printf("# Max. no improv. stretch = %d\n", maximum_no_improvement_stretch);
//	printf("# Fitness var. tolerance  = %e\n", fitness_variance_tolerance);
//	printf("# Random seed             = %ld\n", random_seed);
//	printf("#\n");*/
//
//	int argc = 19;
//	interpretCommandLine(argc, argv);
//
//	Individual* ouput = run();
//
//	closeProblemA();
//
//	return ouput;
//}
//
//Individual* RVGOMEA1A(double** AinvKexp, double** AinvKM52, double** xtrain, double* paramKexp, double* paramKM52, double* ytrain, int n_param,
//	int n_x, int Dimension, double kappa)
//{
//	problem_Solving = 100;
//	//initializeproblem0A(Ainv, xtrain, param, ytrain, n_param, n_x, Dimension, kappa);
//	initializeproblem0A(AinvKexp, AinvKM52, xtrain, paramKexp, paramKM52, ytrain, n_param, n_x, Dimension, kappa);
//	//-b -w -s -r -f 1 -1e2 1e2 0.35 10 25 0.9 1 0 1e-10 100 0.0  10000
//	const char** argv[19];
//	//argv[17]= { "","-s", "-r", "-f", "1", "-1e2", "1e2", "0.35", "10", "25", "0.9", "1", "0", "1e-10", "100", "0.0", "10" };
//	argv[0] = " ";
//	argv[1] = "-b";
//	argv[2] = "-b";
//	//argv[3] = "-s"; generational statistics
//	argv[3] = "-b";
//	argv[4] = "-r";
//	argv[5] = "-f";
//	argv[6] = "1";
//	argv[7] = "-1e2";
//	argv[8] = "1e2";
//	argv[9] = "0.35";
//	argv[10] = "10";
//	argv[11] = "25";
//	argv[12] = "0.9";
//	argv[13] = "1";
//	argv[14] = "0";
//	argv[15] = "-1e4";
//	argv[16] = "100";
//	argv[17] = "0.0";
//	argv[18] = "10";
//
//
//
//	int argc = 19;
//	interpretCommandLine(argc, argv);
//
//	Individual* ouput = run();
//
//	closeProblemA();
//
//	return ouput;
//}*/
double* PBCGSTAB(double** A, double* b, int size)
{
	int iterations = 1000;

	double** M = newDouble(size, size);
	//*M[i,i]=A[i,i];
	//Minv[i,i]=1/A[i,i];
	//M=Minv;
	int i;
	for (int i = 0; i < size; i++)
	{

		if (A[i][i] != 0)
		{
			M[i][i] = 1 / A[i][i];
		}

	}

	double* x = newDoubleV(size);
	//r=b-A*x;
	double* Ax = multVector(A, x, size, size);
	double* r = subVector(b, Ax, size);
	//rbis=r;
	double* rbis = newDoubleV(size);
	for (i = 0; i < size; i++)
	{
		rbis[i] = r[i];
	}


	double ro = 1;
	double roant = 1;
	double w = 1;
	double alpha = 1;

	double* v = newDoubleV(size);
	double* p = newDoubleV(size);

	double res = dotVector(r, r, size);

	int its = 0;
	while (its < iterations && res > 1e-10)
	{
		//?i = (r?0, ri?1)
		ro = dotVector(rbis, r, size);
		//? = (?i/?i?1)(?/?i?1)
		double beta = (ro / roant) * (alpha / w);
		//pi = ri?1 + ?(pi?1 ? ?i?1*vi?1)
		//pi = ri?1 + ?(pi?1 ? ?i?1*vi?1)
		//temp1=?i?1vi?1
		//temp2=pi?1 ? temp1
		//temp3= ?(temp2)
		double* temp1 = multConstantV(v, w, size);
		double* temp2 = subVector(p, temp1, size);
		double* temp3 = multConstantV(temp2, beta, size);
		//pi = ri?1 +temp3
		free(p);
		p = addVector(r, temp3, size);
		//y = K^(?1)*pi
		double* y = multVector(M, p, size, size);
		//vi = Ay
		free(v);
		v = multVector(A, y, size, size);
		//? = ?i/(r?0, vi)

		alpha = ro / dotVector(rbis, v, size);
		//s = ri?1 ? ?vi
		//temp1=?vi
		//s = ri?1 ? temp1
		double* temp4 = multConstantV(v, alpha, size);
		double* s = subVector(r, temp4, size);
		//z = K^(?1)*s
		double* z = multVector(M, s, size, size);
		//t = Az
		double* t = multVector(A, z, size, size);
		//z2=K^(?1)_1*t
		double* z2 = multVector(M, t, size, size);
		//?i = (z2, K^(?1)_1s)/(z2, z2)
		w = dotVector(z2, z, size) / dotVector(z2, z2, size);
		//xi = xi?1 + ?y + ?iz
		double* temp5 = multConstantV(z, w, size);
		double* temp6 = multConstantV(y, alpha, size);
		double* temp7 = addVector(temp5, temp6, size);
		double* xtemp = addVector(x, temp7, size);
		for (i = 0; i < size; i++)
		{
			x[i] = xtemp[i];
		}

		//ri = s ? ?it
		double* temp8 = multConstantV(t, w, size);
		free(r);
		r = subVector(s, temp8, size);
		roant = ro;
		res = dotVector(s, s, size);
		its++;

		free(s);
		free(z);
		free(t);
		free(z2);
		free(y);
		free(temp1);
		free(temp2);
		free(temp3);
		free(temp4);
		free(temp5);
		free(temp6);
		free(temp7);
		free(temp8);
		free(xtemp);
	}


	free(v);
	free(p);
	free(rbis);
	free(r);
	free(Ax);

	for (i = 0; i < size; i++)
	{
		free(M[i]);
	}
	free(M);

	return x;

}

double* PBCGSTABC(double** A, double* b, int size)
{
	int iterations = 1000;
	int i;
	double** M = Cholesky(A, size);

	double* x = newDoubleV(size);
	//r=b-A*x;
	double* Ax = multVector(A, x, size, size);
	double* r = subVector(b, Ax, size);
	//rbis=r;
	double* rbis = newDoubleV(size);
	for (i = 0; i < size; i++)
	{
		rbis[i] = r[i];
	}


	double ro = 1;
	double roant = 1;
	double w = 1;
	double alpha = 1;

	double* v = newDoubleV(size);
	double* p = newDoubleV(size);

	double res = dotVector(r, r, size);

	int its = 0;
	while (its < iterations && res > 1e-10)
	{
		//?i = (r?0, ri?1)
		ro = dotVector(rbis, r, size);
		//? = (?i/?i?1)(?/?i?1)
		double beta = (ro / roant) * (alpha / w);
		//pi = ri?1 + ?(pi?1 ? ?i?1*vi?1)
		//pi = ri?1 + ?(pi?1 ? ?i?1*vi?1)
		//temp1=?i?1vi?1
		//temp2=pi?1 ? temp1
		//temp3= ?(temp2)
		double* temp1 = multConstantV(v, w, size);
		double* temp2 = subVector(p, temp1, size);
		double* temp3 = multConstantV(temp2, beta, size);
		//pi = ri?1 +temp3
		free(p);
		p = addVector(r, temp3, size);
		//y = K^(?1)*pi
		double* y = multVector(M, p, size, size);
		//vi = Ay
		free(v);
		v = multVector(A, y, size, size);
		//? = ?i/(r?0, vi)

		alpha = ro / dotVector(rbis, v, size);
		//s = ri?1 ? ?vi
		//temp1=?vi
		//s = ri?1 ? temp1
		double* temp4 = multConstantV(v, alpha, size);
		double* s = subVector(r, temp4, size);
		//z = K^(?1)*s
		double* z = multVector(M, s, size, size);
		//t = Az
		double* t = multVector(A, z, size, size);
		//z2=K^(?1)_1*t
		double* z2 = multVector(M, t, size, size);
		//?i = (z2, K^(?1)_1s)/(z2, z2)
		w = dotVector(z2, z, size) / dotVector(z2, z2, size);
		//xi = xi?1 + ?y + ?iz
		double* temp5 = multConstantV(z, w, size);
		double* temp6 = multConstantV(y, alpha, size);
		double* temp7 = addVector(temp5, temp6, size);
		double* xtemp = addVector(x, temp7, size);
		for (i = 0; i < size; i++)
		{
			x[i] = xtemp[i];
		}

		//ri = s ? ?it
		double* temp8 = multConstantV(t, w, size);
		free(r);
		r = subVector(s, temp8, size);
		roant = ro;
		res = dotVector(s, s, size);
		its++;

		free(s);
		free(z);
		free(t);
		free(z2);
		free(y);
		free(temp1);
		free(temp2);
		free(temp3);
		free(temp4);
		free(temp5);
		free(temp6);
		free(temp7);
		free(temp8);
		free(xtemp);
	}


	free(v);
	free(p);
	free(rbis);
	free(r);
	free(Ax);

	for (i = 0; i < size; i++)
	{
		free(M[i]);
	}
	free(M);


	return x;

}

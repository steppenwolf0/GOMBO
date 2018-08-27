#ifndef SOLVERSAQ_H
#define SOLVERSAQ_H
#include "problemDef.h"

Individual* RVGOMEAUCB(double** Ainv, double** xtrain, double* param, double* ytrain,
	int n_param, int n_x, int Dimension, double kappa, int kernelType, int evals);

Individual* RVGOMEAEI(double** Ainv, double** xtrain, double* param, double* ytrain, int n_param, int n_x, int Dimension,
	double bestSolution, int kernelType);

double* PBCGSTAB(double** A, double* b, int size);

double* PBCGSTABC(double** A, double* b, int size);

//Individual* RVGOMEA3A(double** AinvKexp, double** AinvKM52, double** xtrain, double* paramKexp, double* paramKM52, double* ytrain, int n_param, int n_x,
//	int Dimension, double bestSolution, double power);
//Individual* RVGOMEA1A(double** AinvKexp, double** AinvKM52, double** xtrain, double* paramKexp, double* paramKM52, double* ytrain, int n_param,
//	int n_x, int Dimension, double kappa);
#endif

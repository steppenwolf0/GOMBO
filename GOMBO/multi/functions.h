#ifndef FUNC_H
#define FUNC_H

double sphereFunctionProblemEvaluation(double *xtrain, int Dimension);
double griewankFunctionProblemEvaluation(double *xtrain, int Dimension);
double rosenbrockFunctionProblemEvaluation(double *xtrain, int Dimension);
double ackleyFunctionProblemEvaluation(double *xtrain, int Dimension);
double michalewiczFunctionProblemEvaluation(double *xtrain, int Dimension);

double cnnFunctionProblemEvaluation(double* xtrain, int Dimension);
double cnnFunctionProblemEvaluationOriginal(double* xtrain, int Dimension);
double cnnFunctionProblemEvaluationBest(double* xtrain, int Dimension, double best);
#endif
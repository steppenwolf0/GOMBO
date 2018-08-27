#ifndef SUPPORT_H
#define SUPPORT_H

double** newDouble(int dim_1M1, int dim_2M2);
double* newDoubleV(int dim);
int* newIntV(int dim);
void freeMatrix(double** matrix, int dim);
double* newDoubleVInit(int dim, double constant);
void copyMatrix(double** matrix1, double** matrix2, int size1, int size2);
void copyVector(double* vector1, double* vector2, int size);
#endif

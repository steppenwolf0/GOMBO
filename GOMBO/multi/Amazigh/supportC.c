#include <stdlib.h>
double** newDouble(int dim_1, int dim_2)
{
	double** matrixdouble = (double **) malloc( dim_1*( sizeof( double* ) ) );
    for( int i = 0; i < dim_1; i++ )
        matrixdouble[i] = (double*) malloc( dim_2*( sizeof( double ) ) );
        
    for ( int i=0;i<dim_1;i++)
            {
                for ( int j = 0; j < dim_2; j++)
                {
                	matrixdouble[i][j]=0.0;
                }
            }
    return matrixdouble;                	
}
double* newDoubleV(int dim)
{
	double* vector = (double*) malloc( dim*( sizeof( double ) ) );
        
        for ( int i=0;i<dim;i++)
            {
            	vector[i]=0.0;
            }
        return vector;
}

double* newDoubleVInit(int dim,double constant)
{
	double* vector = (double*)malloc(dim*(sizeof(double)));

	for (int i = 0; i<dim; i++)
	{
		vector[i] = constant;
	}
	return vector;
}


int* newIntV(int dim)
{
	int* vector = (int*)malloc(dim*(sizeof(int)));

	for (int i = 0; i<dim; i++)
	{
		vector[i] = 0;
	}
	return vector;
}

void freeMatrix(double** matrix, int dim)
{
	int i;
	for (i = 0; i < dim; i++)
	{
		free(matrix[i]);
	}
	free(matrix);
}

void copyVector(double* vector1, double* vector2, int size)
{
	int j;
	for (j = 0; j < size; j++)
	{
		vector1[j] = vector2[j];

	}
}

void copyMatrix(double** matrix1, double** matrix2, int size1, int size2)
{
	int i, j;
	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			matrix1[i][j] = matrix2[i][j];

		}
	}
}
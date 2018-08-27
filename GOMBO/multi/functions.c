#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "./Amazigh/supportC.h"
#include "./Amazigh/openFile.h"
#define PI 3.141592653589793238462643 
#ifdef OSisWindows
#include "./Solvers/util/ToolsWIN.h"
#else
#include "./Solvers/util/Tools.h"
#endif
double cnnFunctionProblemEvaluation(double* xtrain, int Dimension)
{
	initializeRandomNumberGenerator();
	//h1=random.randint(4, 64)
		int h1=xtrain[0]*(64-4)+4;
	//wd1=random.randint(4, 64)
		int wd1=xtrain[1]*(64-4)+4;
	//s1=math.ceil(float(545)/float(h1))
	//print(int(s1))
		int s1=ceil((545)/h1);
		printf("%d \n", s1);
	//h2=random.randint(2, s1)
		int h2=xtrain[2]*(s1-2)+2;
	//wd2=random.randint(2, s1)
		int wd2=xtrain[3]*(s1-2)+2;
	//s2=math.ceil(float(s1)/float(h2))
		int s2=ceil((double)s1/(double)h2);
		printf("%d \n", s2);
	//while (s2<2):
	//	h2=random.randint(2, s1)
	//	wd2=random.randint(2, s1)
	//	s2=math.ceil(float(s1)/float(h2))
	//print(int(s2))
	while (s2<2)
	{
		h2=randomRealUniform01()*(s1-2)+2;
		wd2=randomRealUniform01()*(s1-2)+2;
		s2=ceil((double)s1/(double)h2);
		printf("%d \n", s2);
	}
	int h3=0;
	int wd3=0;
	if (s2==2)
	{
		h3=2;
		wd3=2;
	}
	else
	{
		h3= xtrain[4]*(s2-2)+2;
	 	wd3=xtrain[5]*(s2-2)+2;
	 	int s3=ceil((double)s2/(double)h3);
		printf("%d \n", s3);
	}
		//h3=random.randint(2, s2)
		//wd3=random.randint(2, s2)
	//s3=math.ceil(float(s2)/float(h3))
	//print(int(s3))
			
	int w1=xtrain[6]*(256-2)+2;
	int w2=xtrain[7]*(256-2)+2;
	int w3=xtrain[8]*(256-2)+2;
	int w4= xtrain[9]*(256-2)+2;
	//var='sh gen3.sh 1000 '+str(w1)+' '+str(w2)+' '+str(w3)+' '+str(w4)+' '+str(h1)+' '+str(h2)+' '+str(h3)+' '+str(wd1)+' '+str(wd2)+' '+str(wd3)+' '+str(x)+' '+str(fold)+' '+str(generation)
	
	 char* tempName[100];
		sprintf(tempName, "sh ./gen3.sh 500 %d %d %d %d %d %d %d %d %d %d 0 0 0",w1,w2,w3,w4,h1,h2,h3,wd1,wd2,wd3);
	printf("%d %d %d %d %d %d %d %d %d %d \n",w1,w2,w3,w4,h1,h2,h3,wd1,wd2,wd3);
	system(tempName);
	
	double** result=openMatrix("outputVector.txt");
	double value=result[0][0];
	freeMatrix(result,1);
	return value;
}

double cnnFunctionProblemEvaluationOriginal(double* xtrain, int Dimension)
{
	//h1=random.randint(4, 128)
	int h1=xtrain[0]*(28/2-4)+4;
	//wd1=random.randint(4, 128) 
	int wd1= xtrain[1]*(28/2-4)+4;
	int s1=ceil((28)/h1);
	printf("%d \n", s1);
	//h2=random.randint(2, int(s1/2)+1)
	//wd2=random.randint(2, int(s1/2)+1)
	//int h2= xtrain[2]*(h1-2)+2;
	//int wd2=xtrain[3]*(wd1-2)+2;
	int temp=((s1/2)+1);
	int h2= xtrain[2]*(temp-2)+2;
	int wd2=xtrain[3]*(temp-2)+2;
	int s2=ceil((double)s1/(double)h2);
	printf("%d \n", s2);
	//h3=random.randint(2, int(s2/2)+1)
	//wd3=random.randint(2, int(s2/2)+1)
	//int h3= xtrain[4]*(h2-2)+2;
	//int wd3=xtrain[5]*(wd2-2)+2;
	temp=((s2/2)+1);
	int h3= xtrain[4]*(temp-2)+2;
	int wd3=xtrain[5]*(temp-2)+2;
	int s3=ceil((double)s2/(double)h3);
	printf("%d \n", s3);
	//w1=random.randint(2, 256)
	//w2=random.randint(2, 256)
	//w3=random.randint(2, 256)
	//w4=random.randint(2, 256)
	int w1=xtrain[6]*(256-2)+2;
	int w2=xtrain[7]*(256-2)+2;
	int w3=xtrain[8]*(256-2)+2;
	int w4= xtrain[9]*(256-2)+2;
	//var='sh gen3.sh 1000 '+str(w1)+' '+str(w2)+' '+str(w3)+' '+str(w4)+' '+str(h1)+' '+str(h2)+' '+str(h3)+' '+str(wd1)+' '+str(wd2)+' '+str(wd3)+' '+str(x)+' '+str(fold)+' '+str(generation)
	
	 char* tempName[100];
		sprintf(tempName, "sh ./gen3.sh 500 %d %d %d %d %d %d %d %d %d %d 0 0 0",w1,w2,w3,w4,h1,h2,h3,wd1,wd2,wd3);
	printf("%d %d %d %d %d %d %d %d %d %d \n",w1,w2,w3,w4,h1,h2,h3,wd1,wd2,wd3);
	system(tempName);
	
	double** result=openMatrix("outputVector.txt");
	double value=result[0][0];
	freeMatrix(result,1);
	return value;
}

const double point[100] =  {
	52.34, -85.34, -50.21, 56.2, 89.55, 10.07, -63.64, -77.9, 66.73, 64.45,
		-98.51, -51.52, -41.99, -17.25, 48.74, -44.86, -87, 48.13, 30.18, -53.53, 30.49, 92.69, 10.46,
		52.01, 45.3, -87.43, 0.46, -71.61, 25.69, -49.87, -74.86, -8.04, -9.43, 3.78, 82.61, -57.76,
		11.6, -71.07, 99.42, -11.39, 73.25, 84.24, 22.3, 49.72, 4.86, -53.32, 98.37, -96.25, 96.65,
		37.24, 72.7, 95.84, 26.33, 62.58, 16.86, -86.87, -90.32, 28.12, 20.44, -6.52, 98.74, -60.42,
		4.48, 64.53, 7.89, 21.54, -5.37, 33.11, -10.98, -33.08, 34.87, 0.87, -0.93, 22.94, 55.97, 2,
		43.27, -94.82, -52.26, -65.53, 30.72, -54.08, -78.03, 42.36, -87.16, -85.26, -7.77, 51.21,
		38.54, 27.73, 3.31, 22.61, 87.26, 70.41, -57.51, -19.36, 33.78, -92.23, -35.92, -53.91
};

double griewankFunctionProblemEvaluation(double *xtrain, int Dimension)
{
	
		double* parameters = newDoubleV(Dimension);

	int    i;
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = 200 * xtrain[i] - 100;
	}

	double yi, sum, prod, result;

	sum = 0;
	prod = 1.0;
	for (i = 0; i < Dimension; i++)
	{
		yi = parameters[i] - point[i];
		sum += yi*yi;

		yi = (parameters[i] - point[i]) / sqrt((double)(i + 1));
		prod *= cos(yi);
	}

	result = sum / 4000.0 - prod + 1.0;

	free(parameters);
	return result;
}

double sphereFunctionProblemEvaluation(double *xtrain, int Dimension)
{
	
	double* parameters = newDoubleV(Dimension);

	int    i;
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = 200 * xtrain[i] - 100;
	}
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = parameters[i]-point[i];
	}
	double result;

	result = 0.0;
	for (i = 0; i < Dimension; i++)
		result += parameters[i] * parameters[i];

	free(parameters);
	return result;

}

double rosenbrockFunctionProblemEvaluation(double *xtrain, int Dimension)
{
	
	double* parameters = newDoubleV(Dimension);

	int    i;
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = 200 * xtrain[i] - 100;
	}
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = parameters[i] - point[i];
	}
	double result;

	result = 0.0;
	for (i = 0; i < Dimension - 1; i++)
		result += 100 * (parameters[i + 1] - parameters[i] * parameters[i])*(parameters[i + 1] - parameters[i] * parameters[i]) + (1.0 - parameters[i])*(1.0 - parameters[i]);

	free(parameters);
	return result;
}

double ackleyFunctionProblemEvaluation(double *xtrain, int Dimension) 
{
	

	double* parameters = newDoubleV(Dimension);

	int    i;
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = 200 * xtrain[i] - 100;
	}
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = parameters[i] - point[i];
	}
	double sum1, sum2;
	double result;



	sum1 = 0.0L;
	sum2 = 0.0L;
	for (i = 0; i < Dimension; i++) {
		sum1 += (parameters[i] * parameters[i]);
		sum2 += cos(2.0L * PI * parameters[i]);
	}
	double e = 2.7182818284590452353602874713526625L;
	result = -20.0L * exp(-0.2L * sqrt(sum1 / Dimension)) - exp(sum2 / Dimension) + 20.0L + e;

	free(parameters);
	return result;
}

double michalewiczFunctionProblemEvaluation(double *xtrain, int Dimension)
{
	double* parameters = newDoubleV(Dimension);

	int    i;
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = 200 * xtrain[i] - 100;
	}
	for (i = 0; i < Dimension; i++)
	{
		parameters[i] = parameters[i] - point[i];
	}
	double result;

	result = 0.0;
	for (i = 0; i < Dimension; i++)
		result += -sin(parameters[i])*pow(sin(((i + 1)*parameters[i] * parameters[i]) / PI), 20.0);

	free(parameters);
	return result;
}

double cnnFunctionProblemEvaluationBest(double* xtrain, int Dimension, double best)
{
	int maxSize=100;
	//h1=random.randint(4, 128)
	int h1=xtrain[0]*(100/2-4)+4;
	//wd1=random.randint(4, 128) 
	int wd1= xtrain[1]*(maxSize-2)+2;
	int s1=ceil((100)/h1);
	printf("%d \n", s1);
	//h2=random.randint(2, int(s1/2)+1)
	//wd2=random.randint(2, int(s1/2)+1)
	//int h2= xtrain[2]*(h1-2)+2;
	//int wd2=xtrain[3]*(wd1-2)+2;
	int temp=((s1/2)+1);
	int h2= xtrain[2]*(temp-2)+2;
	int wd2=xtrain[3]*(maxSize-2)+2;
	int s2=ceil((double)s1/(double)h2);
	printf("%d \n", s2);
	//h3=random.randint(2, int(s2/2)+1)
	//wd3=random.randint(2, int(s2/2)+1)
	//int h3= xtrain[4]*(h2-2)+2;
	//int wd3=xtrain[5]*(wd2-2)+2;
	temp=((s2/2)+1);
	int h3= xtrain[4]*(temp-2)+2;
	int wd3=xtrain[5]*(maxSize-2)+2;
	int s3=ceil((double)s2/(double)h3);
	printf("%d \n", s3);
	//w1=random.randint(2, 256)
	//w2=random.randint(2, 256)
	//w3=random.randint(2, 256)
	//w4=random.randint(2, 256)
	int w1=xtrain[6]*(maxSize-2)+2;
	int w2=xtrain[7]*(maxSize-2)+2;
	int w3=xtrain[8]*(maxSize-2)+2;
	int w4= xtrain[9]*(maxSize-2)+2;
	//var='sh gen3.sh 1000 '+str(w1)+' '+str(w2)+' '+str(w3)+' '+str(w4)+' '+str(h1)+' '+str(h2)+' '+str(h3)+' '+str(wd1)+' '+str(wd2)+' '+str(wd3)+' '+str(x)+' '+str(fold)+' '+str(generation)
	
	 char* tempName[100];
		sprintf(tempName, "sh ./pythonBridge.sh 500 %d %d %d %d %d %d %d %d %d %d 0 0 0 %.6f",w1,w2,w3,w4,h1,h2,h3,wd1,wd2,wd3,best);
	printf("50 %d %d %d %d %d %d %d %d %d %d 0 0 0 %.6f \n",w1,w2,w3,w4,h1,h2,h3,wd1,wd2,wd3,best);
	system(tempName);
	
	
	double** result=openMatrix("outputVector.txt");
	double value=result[0][0];
	freeMatrix(result,1);
	return value;
	//return 0;
}
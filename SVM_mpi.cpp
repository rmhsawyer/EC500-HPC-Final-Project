/*	
* EC500 HPC final project: 
* Simplified Sequntial Minimazation Optimazation solver for SVM problem
* Wrote by Minghe Ren, Xinqiao Wei.
* 2018 Boston Univerisity
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <array>
#include <time.h>   
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <mpi.h>
using namespace std;

#define C 1.6
#define tolerance 0.01

void read_X(string x_file, double** x);
void read_Y(string y_file, double* y);

double f_x(int j, double** x, double* y, double* alpha, double b);
double dot_product(double** x, int i, int j);


double L_nequal(double* alpha, int i, int j);
double H_nequal(double* alpha, int i, int j);
double L_equal(double* alpha, int i, int j);
double H_equal(double* alpha, int i, int j);

double eta(double** x, int i, int j);
void alpha_j(double* alpha, int j, double H, double L);
double b_value(int i, int j, double b, double* E, double** x, double* y, double* alpha, double* alpha_old);

int main(int argc, char** argv)
{
	//we random generate 15 testing data points to run the code, each data point with label +1, -1 
	double *a;
	double *a_old;
	a= new double [15];
	a_old= new double [15];
	double b=0;
	double E[15];
	//passes
	int passes = 0;
	int max_passes = 100;
	
	int i,j;
	double L,H,theta;
	srand (time(NULL));

	double **x;
	double *y;
	x = new double *[15];
	y = new double [15];

	for (i=0;i<15;i++)
	{
	    x[i] = new double [2];
	}
	
	for(i= 0;i<15;i++)
	{
		a[i] = 0;
		a_old[i] = 0;
	}

	int world_size; // number of MPI processes
    int world_rank; // rank of MPI process
    char processor_name[MPI_MAX_PROCESSOR_NAME]; // name of processor
    int name_len; // length of name string.


    // Initialize MPI
    MPI_Init(&argc, &argv);
    
    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    // Get the rank
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  
    // Get the name of the processor
    MPI_Get_processor_name(processor_name, &name_len);

	// read in data
	printf("The data point coordinates are:\n");
	read_X("generate_data/dataset_15.txt", x);
	printf("The data point labels are:\n");
 	read_Y("generate_data/dataset_15_label.txt", y);

 	// solve the optimazation problem, the loop stops while a doesn't change 
	while(passes <max_passes)
	{	
		int num_changed_alphas = 0;

		for( i =world_rank*(15-1)/world_size;i<(world_rank+1)*(15-1)/world_size;i++)
		{	
			E[i] = f_x( i,x,y,a,b) - y[i];
			if( ( y[i] * E[i]< -tolerance && a[i] < C ) || ( y[i] * E[i]> tolerance && a[i] >0 ) )
			{
				//randomly choose j != i
				do{
					j= rand() % 15;
				}while(j ==i);

				E[j] = f_x( j,x,y,a,b) - y[j];
				a_old[i] = a[i];
				a_old[j] = a[j];
				if(y[i] != y[j])
				{
					H = H_nequal(a,i,j);
					L = L_nequal(a,i,j);
				}
				if(y[i] == y[j])
				{
					H = H_equal(a,i,j);
					L = L_equal(a,i,j);
				}

				if(L != H)
				{

					theta = eta(x,i,j);
					if(theta < 0)
					{
						a[j] =  a[j]- y[j]*(E[i]-E[j])/theta;
						alpha_j(a,j,H,L);

						if(abs(a[j]-a_old[j])>=1e-5)
						{
							a[i]= a[i]+y[i]*y[j]*(a_old[j]-a[j]);
							b = b_value(i,j,b,E,x,y,a,a_old);
							num_changed_alphas = num_changed_alphas +1;

						}
					}
				}
			}
		}

		if(num_changed_alphas ==0)
		{
			passes = passes+1;
		}
		else
		{
			passes =0;
		}
	}

	// Clean up
    MPI_Finalize();

	//Compute W
	double w_svm[2];
	printf("The calculated alphas are:\n");
	for(i= 0;i<15;i++)
	{
		w_svm[0]= w_svm[0] +a[i] * y[i]*x[i][0] ;
		w_svm[1]= w_svm[1] +a[i] * y[i]*x[i][1] ;

		cout<<a[i]<<"  ";
	}
	printf("\nThe calculated W is:\n");
	cout<<w_svm[0]<<' '<<w_svm[1]<<endl;

	//Compute b
	double b_svm;
	b_svm = y[14] - (w_svm[0]*x[14][0]+w_svm[0]* x[14][1]);
	printf("The calculated b is:\n");
	cout<<b_svm<<endl;

	return 0;
}

// Compute f(x) dual problem
double f_x(int j, double** x, double* y, double* alpha, double b){
 int m = 15;
 double sum = 0;
 for (int i = 0; i < m; i++){
  sum += alpha[i] * y[i] * dot_product(x,i,j);
 }
 sum += b;
 return sum;
}

// Compute dot product
double dot_product(double** x, int i, int j){
 double sum1, sum2;
 sum1 = x[i][0] * x[j][0];
 sum2 = x[i][1] * x[j][1];
 return sum1+sum2;
}


// Find boundary of ai and aj
double L_nequal(double *alpha, int i, int j){
 double tmp = alpha[j] - alpha[i];
 if (tmp > 0)
 {
  return tmp;
 }
 else
 {
  tmp = 0;
  return tmp;
 }
}

double H_nequal(double *alpha, int i, int j){
 double tmp = C + alpha[j] - alpha[i];
 if (tmp > C)
 {
  return C;
 }
 else
 {
  return tmp;
 }
}

double L_equal(double *alpha, int i, int j){
 double tmp = alpha[j] + alpha[i] - C;
 if (tmp > 0)
 {
  return tmp;
 }
 else
 {
 	tmp = 0;
  return tmp;
 }
}

double H_equal(double *alpha, int i, int j){
 double tmp = alpha[j] + alpha[i];
 if (tmp > C)
 {
  return C;
 }
 else
 {
  return tmp;
 }
}

// Compute theta
double eta(double** x, int i, int j){
 double t1 = 2 * dot_product(x, i, j);
 double t2 = dot_product(x, i, i);
 double t3 = dot_product(x, j, j);
 return t1 - t2  -t3;
}

// Compute aj to maximize the objective function
void alpha_j(double* alpha, int j, double H, double L){
 if (alpha[j] > H)
 {
  alpha[j] = H;
 }
 else if (alpha[j] < L)
 {
  alpha[j] = L;
 }
}

// Compute the b shreshold
double b_value(int i, int j, double b, double* E, double** x, double* y, double* alpha, double* alpha_old)
{
 double t1, t2, b1,b2;
 t1 = y[i] * (alpha[i] - alpha_old[i]) * dot_product(x, i, i);
 t2 = y[j] * (alpha[j] - alpha_old[j]) * dot_product(x, i, j);
 b1 = b - E[i] - t1 - t2;

 t1 = y[i] * (alpha[i] - alpha_old[i]) * dot_product(x, i, j);
 t2 = y[j] * (alpha[j] - alpha_old[j]) * dot_product(x, j, j);
 b2 = b - E[j] - t1 - t2;

 if (alpha[i] < C && alpha[i] > 0)
 {
  b = b1;
 }
 else if (alpha[j] < C && alpha[j] > 0)
 {
  b = b2;
 }
 else
 {
  b = (b1 + b2) / 2;
 }

 return b;
}


// Function to read in data
void read_X(string x_file, double** x){

 string line;
 ifstream input(x_file);
 int i = 0;
    for  (int i = 0; i < 15; i++)
    {
        double x1, x2;
        char separator;
        cout<<i<<" ";
        input >> x1  >> x2;
        x[i][0] = x1;
        x[i][1] = x2;
        cout << x[i][0] << ", " << x[i][1] << endl;
    }
}
void read_Y(string y_file, double* y){

 string line;
 ifstream input(y_file);
 int i = 0;
    for  (int i = 0; i < 15; i++)
    {
        double y1;
        char separator;
        cout<<i<<" ";
        input >> y1;
        y[i] = y1;

        cout << y[i] << endl;
    }
}








#include<iostream>
#include<cmath>
#include<omp.h>
#include"CGSolver.hpp"

// nvc++ -o main main.cpp CGSolver.cpp   
// nvc++ -ta=multicore -Minfo=accel -o main main.cpp CGSolver.cpp
// nvc++ -ta=tesla:managed -Minfo=accel -o main main.cpp CGSolver.cpp
// nvc++ -ta=tesla -Minfo=accel -o main main.cpp CGSolver.cpp
CGSolver::CGSolver(int size, int nnzc,int max_itrations, double tolerence, int *rowPtr,int *colPtr, double *nonZeros, double *f,double *u)
{
       
	maxItr = max_itrations;
	tol = tolerence;	
	N = size;
	Rptr = rowPtr;
	Cptr = colPtr;
	NNZ = nonZeros;
	int nonZeroscount = nnzc; 
	b = f;
	scalar_alpha = 0.0;
	scalar_beta = 0.0;
	
	x = u;
	r_prev = new double[N];
	 
	p = new double[N];
	r_new = new double[N];
	temp = new double[N];
	Ap = new double[N];
	
	for(int i = 0; i < N; i++)
	{
	  x[i] = 0.0;
	  r_prev[i] = 0.0;
	  r_new[i] = 0.0;
	  p[i] = 0.0;
	  temp[i] = 0.0;
	}
	#pragma acc data copyin(x[0:N],r_prev[0:N],r_new[0:N],p[0:N],temp[0:N],Ap[0:N],b[0:N],Rptr[0:N+1],Cptr[0:nonZeroscount],NNZ[0:nonZeroscount])
}
void CGSolver::MatVecProduct(int *rowPtr, int *colPtr, double *nonZeros, double *vec, double *out)
{
  #pragma acc parallel loop independent
  for(int i = 0; i < N; i++)
  {
    int start = rowPtr[i];
    int end   = rowPtr[i+1];
    out[i] = 0.0;
    //std::cout<<"start = "<<start<<"End = "<<end<<std::endl;
    for(int j = start; j < end; j++)
    {
      out[i] += nonZeros[j] * vec[colPtr[j]];
      //std::cout<<"out["<<i<<"] = "<<nonZeros[j]<<" * "<<vec[colPtr[j]]<<std::endl;
    }
    //std::cout<<"*****";
  }
  
}

void CGSolver::vecSubtract(double *vec1, double *vec2, double *out)
{
  #pragma acc parallel loop 
  for(int i = 0; i < N; i++)
  {
   out[i] = vec1[i] - vec2[i];
  }
}

double CGSolver::dotProduct(double *vec1, double *vec2)
{
 double sum = 0.0;
 #pragma acc parallel loop reduction(+:sum)   
 for(int i = 0; i < N; i++)
 {
   sum += (vec1[i] * vec2[i]);
 }
 return sum;
}

double CGSolver::L2Norm(double *vec1)
{
  double error = 0.0;
  #pragma acc parallel loop reduction(+:error)
  for(int i = 0; i < N; i++)
  {
   error += vec1[i]*vec1[i];
  }
  return sqrt(error);
}

void CGSolver::scalarVecProduct(double scalar, double *vec, double *temp)
{
  #pragma acc parallel loop
  for(int i = 0; i < N; i++)
  {
   temp[i] = scalar * vec[i];
  }
}
void CGSolver::vecAddition(double *vec1, double *vec2, double *out)
{
  #pragma acc parallel loop
  for(int i = 0; i < N; i++)
  {
   out[i] = vec1[i] + vec2[i];
  }
}

void CGSolver::update_prev_residue(double * new_r, double *prev_r)
{
 #pragma acc parallel loop
 for(int i = 0; i < N; i++)
 {
   prev_r[i] = new_r[i];
 }
}

void CGSolver::initVecP(double *vecP, double *vec){
#pragma acc parallel loop
for(int i = 0; i < N; i++)
 {
   vecP[i] = vec[i];
 }

}

void CGSolver::printVector(double *vec)
{
 for(int i = 0; i < N; i++)
 {
   std::cout<<"i = "<<i<<" "<<vec[i]<<std::endl;
 }
}

void CGSolver::update_x(double *x, double *p, double alpha)
{
 #pragma acc parallel loop
 for(int i = 0; i < N; i++)
 {
   x[i] = x[i] + alpha * p[i];
 }
}

void CGSolver::update_r(double *r, double *temp, double alpha)
{
 #pragma acc parallel loop
 for(int i = 0; i < N; i++)
 {
   r[i] = r[i] - alpha * temp[i];
 }
}

void CGSolver::update_p(double *p, double *r, double beta)
{
  #pragma acc parallel loop
  for(int i = 0; i < N; i++)
  {
    p[i] = r[i] + beta * p[i];
  }
}


 void CGSolver::solve()
 {
 
 int k = 0;
 initVecP(r_new,b);
 
 double error = L2Norm(r_new);
 std::cout<<"Error "<<error<<std::endl;
 
 while(k < N)
 {
   if(k > 0)
   {
     scalar_beta = dotProduct(r_new,r_new)/dotProduct(r_prev,r_prev);
   }
    
   update_p(p,r_new,scalar_beta);
   MatVecProduct(Rptr,Cptr,NNZ,p,Ap);
   
   scalar_alpha = dotProduct(r_new,r_new)/dotProduct(p,Ap);
   
   update_x(x,p,scalar_alpha);
   update_prev_residue(r_new,r_prev);
   MatVecProduct(Rptr,Cptr,NNZ,x,temp);
   vecSubtract(b,temp,r_new);
   error = L2Norm(r_new);
   if(error < tol)
   {
    break;
   }
   if(k%100 == 0)
   {
   // std::cout<<"Iteration: "<<k<<" Error "<<error<<std::endl;   
   }
   k = k + 1;     
 }
 #pragma acc data copyout(x[0:N])
//printVector(x);

}



/*
 
*/

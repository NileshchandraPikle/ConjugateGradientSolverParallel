#include<iostream>
#include<fstream>
#include"CGSolver.hpp"

int main(int argc, char *argv[])
{  
  int rows, cols,nnz;
  std::ifstream myfile;
  myfile.open ("out.txt");
  myfile >> rows >> cols>> nnz;
 
 std::cout<<"Number of non zeros "<<nnz<<std::endl;
 std::cout<<"Number of rows "<<rows<<std::endl;
 std::cout<<"Number of columns "<<cols<<std::endl;
  
 int *rptr_ld     = new int[nnz];
 int *cptr_ld     = new int[nnz];
 double *nonZeros_ld = new double[nnz];
 
 int count_diagonal = 0;
 int count_non_diagonal = 0;
 
 for(int i = 0; i < nnz; i++)
 {
  myfile>>rptr_ld[i];
  myfile>>cptr_ld[i];
  myfile>>nonZeros_ld[i];
  if(rptr_ld[i] == cptr_ld[i])
  {
    count_diagonal++;
  }
 }
 count_non_diagonal = nnz - count_diagonal;
 std::cout<<"Total Diagonal non-zero Elements = "<<count_diagonal<<std::endl;
 std::cout<<"Total Non Diagonal non-zero Elements = "<<count_non_diagonal<<std::endl;
 std::cout<<"Total non-zero Elements = "<<count_non_diagonal + count_diagonal<<"  From File: "<<nnz<< std::endl;
 /*for(int i = 0; i < nnz; i++)
 {
 std::cout<<rptr_ld[i]<<"  " << cptr_ld[i] <<"  "<< nonZeros_ld[i] <<std::endl;
 }*/
 
 int *R_COO = new int[nnz+count_non_diagonal];
 int *C_COO = new int[nnz+count_non_diagonal];
 double *N_COO = new double[nnz+count_non_diagonal];
  for(int i = 0; i < nnz; i++)
  {
     R_COO[i] = rptr_ld[i]-1;
     C_COO[i] = cptr_ld[i]-1;
     N_COO[i] = nonZeros_ld[i];
  }
 
 int k = nnz;
 for(int i = 0; i < nnz; i++)
 {
   if(rptr_ld[i] != cptr_ld[i])
   {
     R_COO[k] = cptr_ld[i]-1;
     C_COO[k] = rptr_ld[i]-1;
     N_COO[k] = nonZeros_ld[i];
     k++;
   }
 }
 
 /*std::cout<<"**********************"<<std::endl;
 for(int i = 0; i < nnz+count_non_diagonal; i++)
 {
 std::cout<<R_COO[i]<<"  " << C_COO[i] <<"  "<< N_COO[i] <<std::endl;
 }*/


 for(int i = 0; i < nnz+count_non_diagonal; i++)
 {
   for(int j = 0; j < nnz+count_non_diagonal - i -1 ; j++)
   {
     if(R_COO[j] > R_COO[j+1])
     {
       double temp1 = R_COO[j];
       R_COO[j] = R_COO[j+1];
       R_COO[j+1] = temp1;
       
       double temp2 = C_COO[j];
       C_COO[j] = C_COO[j+1];
       C_COO[j+1] = temp2;
       
       double temp3 = N_COO[j];
       N_COO[j] = N_COO[j+1];
       N_COO[j+1] = temp3;
       
     }
   }
 }

 std::cout<<"**********************"<<std::endl;
 for(int i = 0; i < nnz+count_non_diagonal; i++)
 {
 //std::cout<<R_COO[i]<<"  " << C_COO[i] <<"  "<< N_COO[i] <<std::endl;
 }
 
 /*********************** COO to CSR Coversion *****************************/
    double *N_CSR = new double[nnz+count_non_diagonal];
    int *C_CSR = new int[nnz+count_non_diagonal];
    
    int *R_CSR = new int[rows+1];
    
    for (int i = 0; i < nnz+count_non_diagonal; i++)
    {
        N_CSR[i] = N_COO[i];
        C_CSR[i] = C_COO[i];
        R_CSR[R_COO[i]+1]++;
    }
    for (int i = 0; i < rows; i++)
    {
        R_CSR[i + 1] += R_CSR[i];
    }
 
// std::cout<<"**********************"<<std::endl;
 for(int i = 0; i < rows+1; i++)
 {
 //std::cout<<R_CSR[i]<<"  " << "  "<<R_COO[i]<<"  "<<C_CSR[i] <<"  "<< N_CSR[i] <<std::endl;
 }
  
 
 
 double *x = new double[rows];
 for(int i = 0; i < rows; i++)
 {
  x[i] = 1.0;
 }
 
 double *f = new double[rows];
 
 for(int i = 0; i < rows; i++)
 {
  f[i] = 0.0;
 }
 
 for(int i = 0; i < rows; i++)
  {
    int start = R_CSR[i];
    int end   = R_CSR[i+1];
    //std::cout<<"start "<<start<<std::endl;
    //std::cout<<"end "<<end<<std::endl;
    for(int j = start; j < end; j++)
    {
      f[i] += N_CSR[j] * x[C_CSR[j]];
      //std::cout<<  N_CSR[j] * x[C_CSR[j]] <<" = "<<N_CSR[j]<< "*"<< x[C_CSR[j]]<<std::endl;
    }
   // std::cout<<"** "<<i<<" "<<f[i]<<std::endl;
  }

  // A x = b

 
 
  /*
    Input to CG Solver: 
    	1. Size of vectors or global stiffness matrix 
    	2. Maximum number of iterations
    	3. Global Matrix K
    	4. RHS vector f
  */
  int nnz_count = nnz+count_non_diagonal;
  
  double *u = new double[rows];
  for(int i = 0; i < rows; i++)
  {
   u[i] = 0.0;
  }
  CGSolver S(rows,nnz_count, 4, 10e-2, R_CSR,C_CSR, N_CSR,f,u);
  S.solve();
  for(int i = 0; i < rows; i++)
  {
   std::cout<<"u["<<i<<"] = "<<u[i]<<std::endl;
  }
  
  return 0;
}

#ifndef SRC_CGSOLVER_HPP
#define SRC_CGSOLVER_HPP

class CGSolver{
      
public:
	int maxItr;
	int N;
	int nonzeroCount;
	double tol;
	double scalar_alpha, scalar_beta;
	int *Rptr;
	int *Cptr;
	double *NNZ;
	double *x;
	double *b;
	double *temp;
	
	double *r_prev;
	
	double *p;
	double *r_new;
	double *Ap;
	   
	CGSolver(int size, int nnz_count, int max_itrations, double tolerence, int *rowPtr,int *colPtr, double *nonZeros, double *f,double *u);
	
	void MatVecProduct(int *rowPtr, int *colPtr, double *nonZeros, double *vec, double *out);
	
	double L2Norm(double *r);
	double dotProduct(double *vec1, double *vec2);
	void scalarVecProduct(double scalar, double *vec, double *temp);
	void vecAddition(double *vec1, double *vec2, double *out);
	void vecSubtract(double *vec1, double *vec2, double *out);
	void update_prev_residue(double * new_r, double *prev_r);
	void initVecP(double *vecP, double *vec);
        void solve();
        void printVector(double *vec);
        void update_x(double *x, double *p, double alpha);
        void update_r(double *r, double *temp, double alpha);
        void update_p(double *p, double *r, double beta);
        void resetTemp(double *temp);
};

#endif

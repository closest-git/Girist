#ifndef _GIP_SS_H_
#define _GIP_SS_H_

/*
	SPARSE SOLVER for GIPACK
	v0.1	cys
		8/24/2008

目标
	千万阶稀疏矩阵求解，余量降低一个数量级以上

	基于algebraic multigrid框架
	支持矩阵的显式及隐式表达
	支持LU分解，Krylov迭代等

*/
typedef enum {	
	//0xFFFF	保留为GSS中的m_type
	GIP_COLUMN_MAJOR=0x10000,GIP_ROW_MAJOR=0x20000
}GIP_MATRIX_TYPE;

typedef enum {		//ALGORITHM TYPE for Meyer decomposition
	GIP_SS_ALG_UNKNOWN=0x0,
	GIP_SS_DIRECT=0x100,GIP_SS_ITER=0x200,			//sparse solver	
}GIP_SS_ALG;

#ifdef __cplusplus
extern "C"  {
#endif


	void DSCAL( const int *, const double *, const double *, const int * );
	double DNRM2( int *, const double *, int *);
	void DGEMV( char*,const int *,const int *, const double *, const double *,const int *, const double *,const int *, const double *, double *,const int *);
	void DGEMM( char*,char*, const int *,const int *,const int *, const double *, const double *, const int *, const double *, const int *, const double *, double *, const int *);

	#define BLAS_SCAL_0( n, a, x, incx )			DSCAL( &(n), &(a), x, &(incx) );	
	#define BLAS_NRM2_0( n, x, incx )				DNRM2( &(n), x, &(incx) );															 
	#define BLAS_GEMV_0(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)	\
	{	int _M_=m,_N_=n;													\
		DGEMV( trans, &(_M_), &(_N_), &(alpha), a, &(lda), x, &(incx), &(beta), y, &(incy) );		}				
	#define BLAS_GEMM_0( transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc ) \
	{ 	int _M_=m,_N_=n,_K_=k;													\
		DGEMM( transa,transb,&(_M_),&(_N_),&(_K_),&(alpha),a,&(lda),b,&(ldb),&(beta),c,&(ldc) );		}		


typedef struct GIP_SS_tag GIP_SPARSE_SOLVER;
struct GIP_SS_tag{
	GIP_SS_ALG alg;
	int dim,nz,nz_max,*ptr,*ind,m_type,*I_temp;			//
	double *x,*rhs,*val,*temp,res;
	void *hLU;
};

	void *GIP_SS_init_( int dim,int nz_max,int *_ptr,int *_ind,double *_val,int m_type,int flag );
	void GIP_SS_solve_( GIP_SPARSE_SOLVER *hSS,int flag );
	void GIP_SS_update_( GIP_SPARSE_SOLVER *hSS,GIP_BOUNDARY_CONDITION *type,double *a );
	void GIP_SS_clear_( GIP_SPARSE_SOLVER *hSS );

/*
	KRYLOV SPARSE SOLVER
	
	KRYLOV SPARSE SOLVER原则上采用矩阵的隐式表达，暂保留(ptr,ind,val）结构
*/
typedef struct GIP_KS_tag GIP_KRYLOV_SOLVER;
struct GIP_KS_tag  {
	char enviro_str[200];
	int info,nRestart,MAX_RESTART;
	int dim,*ptr,*ind,m_type;
	double *val;
	int method,nEV,nGV,nRV,ldh,nCurEV,V_m,V_p,ldq,*i_temp;
	GIP_FLOATING *V,*H,*x,*rhs;
	double *d_temp,tole,ortho_thresh,res;
};

	GIP_KRYLOV_SOLVER* GIP_Krylov_init_( int dim,int m,int flag,int *ptr,int*ind,double *val  );
	void GIP_Krylov_update_( GIP_KRYLOV_SOLVER* hSolver,int j,int flag  );
	void GIP_Krylov_core_( GIP_KRYLOV_SOLVER* hKS,GIP_FLOATING *x_0,GIP_FLOATING *rhs,int flag  );

#ifdef __cplusplus
}
#endif

#endif
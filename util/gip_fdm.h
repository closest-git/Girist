#ifndef _GIP_FDM_H_
#define _GIP_FDM_H_

/*
	FINITE DIFFERENCE METHOD

*/
#ifdef __cplusplus
extern "C"  {
#endif
/*
	应包含线性,非线性两种方程组
*/
	typedef struct GIP_FD_EQUATION_tag GIP_FD_EQUATION;
	struct GIP_FD_EQUATION_tag{
		int dim,nz,*ptr,*ind,type;			
		GIP_FLOATING *val;
	};

	double FD_R1_0( int isBndry,GIP_FD_METHOD fdm,double *u,int pos,int h ); 
//	double FD_R2_0( int tag,GIP_FD_METHOD fdm,double *u,int pos,int h );
	double FD_LAPLAS_0( int tag,GIP_FD_METHOD fdm,double *u,int pos,int h );	

	int FDM_init_equations_m_( int M,int N,int ldu,GIP_FD_METHOD fdm,int flag,GIP_FLOATING alpha,GIP_FLOATING diag,GIP_FD_EQUATION *equation );

#ifdef __cplusplus
}
#endif

#endif
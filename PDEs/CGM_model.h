#ifndef _CGM_MODEL_H_
#define _CGM_MODEL_H_
/*
	[REF] A nonlinear primal-dual method for total variation-based image restoration (1999)
	by Tony F.Chan G.H.Golub and P.Mulet

	v0.1	8/17/2008
		cys
*/
typedef enum {		//ALGORITHM TYPE for Meyer decomposition
	CGM_ALG_UNKNOWN=0x0,
	CGM_METHOD_DUAL=0x01,CGM_METHOD_MEDIAN=0x02,								//
	CGM_SS_DIRECT=0x100,CGM_SS_ITER=0x200,			//sparse solver	
	CGM_A_AT=0x1000,CGM_A_BOUNDARY=0x2000,			//modification to A
}CGM_MODEL_ALG;

/*
	[REF] Modeling Textures with Total Variation Minimization and Oscillating Patterns in Image Processing
	[REF] Image denoising and decomposition with total variation. minimization and oscillatory functions

	v0.1	8/4/2008
		cys
*/
typedef struct CGM_MODEL_tag CGM_MODEL;
struct CGM_MODEL_tag{
	CGM_MODEL_ALG alg;
	GIP_FD_METHOD fdm;
	int M,N,ldu,shift,type,*tag,flag;
	int A_dim,A_nz,*A_ptr,*A_ind,m_type;			//
	double *u,*w,*rhs,*f,*g,*A_val,*temp,lenda,alpha,beta;
	double nn_res,g_norm;					//nonlinear residual of newton method
	void *hSS;
};

#ifdef __cplusplus
extern "C"  {
#endif
	void CGM_MODEL_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName );

#ifdef __cplusplus
}
#endif

#endif

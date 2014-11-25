#ifndef _OSV_MODEL_H_
#define _OSV_MODEL_H_
/*
	[REF] Image decomposition and restoration using total variation minimization and the H-1 norm
	by Stanley Osher, Andrs Sol, Luminita Vese

	v0.1	8/27/2008
		cys
*/
typedef enum {		//ALGORITHM TYPE for Meyer decomposition
	OSV_ALG_UNKNOWN=0x0,
	OSV_ALG_LP=0x01,								//0xFF:		space	
	OSV_SS_DIRECT=0x100,OSV_SS_ITER=0x200,			//0xFF00:	sparse solver		
}OSV_MODEL_ALG;

/*
	[REF] Modeling Textures with Total Variation Minimization and Oscillating Patterns in Image Processing
	[REF] Image denoising and decomposition with total variation. minimization and oscillatory functions

	v0.1	8/4/2008
		cys
*/
typedef struct OSV_MODEL_tag OSV_MODEL;
struct OSV_MODEL_tag{
	OSV_MODEL_ALG alg;
	GIP_FD_METHOD fdm;
	int M,N,ldu,shift,type,*tag,flag;
	int A_dim,A_nz,*A_ptr,*A_ind,m_type;			//
	double *u,*g,*k,*w,*rhs,*f,*A_val,*temp,lenda,beta;
	double nn_res,g_norm;					//nonlinear residual of newton method
	void *hSS;
};

#ifdef __cplusplus
extern "C"  {
#endif
	void OSV_MODEL_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName );

#ifdef __cplusplus
}
#endif

#endif

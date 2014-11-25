#ifndef _MEYER_IDEA_H_
#define _MEYER_IDEA_H_
/*
	IMAGE decomposition:	following the ideas of Yves Meyer, 

	v0.1	8/4/2008
		cys
*/
typedef enum {		//ALGORITHM TYPE for Meyer decomposition
	MEYER_ALG_UNKNOWN=0x0,
	MEYER_ALG_LP=0x01,								//space	0xFF
	MEYER_SS_DIRECT=0x100,MEYER_SS_ITER=0x200,		//sparse solver		
}MEYER_SPLIT_ALG;

/*
	[REF] Modeling Textures with Total Variation Minimization and Oscillating Patterns in Image Processing
	[REF] Image denoising and decomposition with total variation minimization and oscillatory functions

	v0.1	8/4/2008
		cys
*/
typedef struct MEYER_SPLIT_tag MEYER_SPLIT;
struct MEYER_SPLIT_tag{
	MEYER_SPLIT_ALG alg;
	int M,N,ldu,shift,type,*tag,flag;
	int A_dim,A_nz,*A_ptr,*A_ind,m_type;			//
	double *u,*rhs,*f,*A_val,lenda,c1,c2,mui,dt,eps,beta;	
	double *g,*g1,*g2,g_max,gx_max,gy_max;
	void *hSS;
};

#ifdef __cplusplus
extern "C"  {
#endif
	void MEYER_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName );

#ifdef __cplusplus
}
#endif

#endif

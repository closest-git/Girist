#ifndef _TOTAL_VARIATION_MODEL_H_
#define _TOTAL_VARIATION_MODEL_H_

/*
	for tat of TV_MODEL	
	RANGE:	0xFF0000
*/
#define GIP_TVM_ADAPTIVE_P	0x100000

/*
	Solver for Total Variation(TV) model	

	v0.1	5/21/2008
		cys
*/
typedef struct TV_MODEL_tag TV_MODEL;
struct TV_MODEL_tag{	
	int M,N,ldu,shift,type,*tag,flag;
	int A_dim,A_nz,*A_ptr,*A_ind,m_type;			//
	double *u,*z,*A_val,*p,alpha,beta;	
	void *hSS;			//sparse solver
};

#ifdef __cplusplus
extern "C"  {
#endif
	void TVM_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName );

#ifdef __cplusplus
}
#endif

#endif

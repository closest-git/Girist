#ifndef _CHAN_VESE_MODEL_H_
#define _CHAN_VESE_MODEL_H_
/*
	Solver for Chan¨CVese(CV) model	
	region-based models implemented via level-set techniques 

	v0.1	5/12/2008
		cys
*/
typedef enum {		//ALGORITHM TYPE for SIGNED DISTANCE FUNCTION
	CVM_ALG_0,CVM_GLOBAL_CONVEX
}GIP_CVM_ALG;

typedef struct CV_MODEL_tag CV_MODEL;
struct CV_MODEL_tag{
	int M,N,ldu,shift,type,*tag,flag;
	int A_dim,A_nz,*A_ptr,*A_ind,m_type,alg;			//
	double *u_buffer,*z,*A_val,lenda1,lenda2,c1,c2,mui,dt,eps;	
	void *hSS;
};

#ifdef __cplusplus
extern "C"  {
#endif
	void CVM_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName );

#ifdef __cplusplus
}
#endif

#endif

#ifndef _GIP_DECOMPOSITION_H_
#define _GIP_DECOMPOSITION_H_

/*
	v0.1	10/6/2008
		cys
*/



typedef struct GIP_DECOMP_tag GIP_DECOMP;
struct GIP_DECOMP_tag{
	GIP_OBJECT *hObj;
	int shift,type,flag;
	double *temp,lenda,alpha,beta;
	double nn_res,g_norm;					
};

#ifdef __cplusplus
extern "C"  {
#endif
	void GIP_DECOMP_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,int model,int ,char* sWndName );
	int _DECOM_Hs_force( int M,int N,int ldu,double*u,double *f,int flag );	
#ifdef __cplusplus
}
#endif

#endif

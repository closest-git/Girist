#ifndef _GIP_DIFFUSION_H_
#define _GIP_DIFFUSION_H_

/*
	v0.1	10/6/2008
		cys
*/



typedef struct GIP_DIFFUSION_tag GIP_DIFFUSION;
struct GIP_DIFFUSION_tag{
	int M,N,ldu,shift,type,*tag,flag;
	double *u,*f,*g,*temp,lenda,alpha,beta;
	double nn_res,g_norm;					//nonlinear residual of newton method
};

#ifdef __cplusplus
extern "C"  {
#endif
	void GIP_DIFFUSION_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName );

#ifdef __cplusplus
}
#endif

#endif

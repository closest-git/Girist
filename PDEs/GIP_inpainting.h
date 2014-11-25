#ifndef _GIP_INPAINTING_H_
#define _GIP_INPAINTING_H_

typedef double INP_FLOAT;
/*
	for flag of GAC_SOLVER
	����opencv��Լ����u,g�����ȴ洢

*/

typedef struct GIP_INPAINT_tag GIP_INPAINT;
struct GIP_INPAINT_tag{
	int M,N,ldr,type,banner,*flag;
	INP_FLOAT *T;
	double *u_buffer;//*u;
};

#ifdef __cplusplus
extern "C"  {
#endif
	void GAC_basic( GIP_IMAGE *hIMG );

#ifdef __cplusplus
}
#endif

#endif

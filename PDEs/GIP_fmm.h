#ifndef _GIP_FMM_H_
#define _GIP_FMM_H_

typedef double FMM_FLOAT;
/*
	fast marching method
	依照opencv的约定，u,g行优先存储
	利用tag的16-31位，0xFF00
	需将其移到高位			cys	5/13/2008
*/
#define GIP_FMM_BAND	0x1000
#define GIP_FMM_KNOWN	0x2000
#define GIP_FMM_SIDE_0	0x000		//未定
#define GIP_FMM_SIDE_1	0x100		//某一侧
#define GIP_FMM_SIDE_2	0x200		//另外一侧

#define GIP_FMM_OO		3.0+38F		//FLT_MAX:3.402823466e+38F        

#ifdef __cplusplus
extern "C"  {
#endif
	int GIP_FMM_update( int M,int N,FMM_FLOAT *T,int *flag,int ldr,GIP_HEAP *heap,int side,int );
	void GIP_FMM_smooth_u_( int M,int N,double *u,int ldu,int *tag,int flag );

#ifdef __cplusplus
}
#endif

#endif

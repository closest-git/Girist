#ifndef _GIP_FMM_H_
#define _GIP_FMM_H_

typedef double FMM_FLOAT;
/*
	fast marching method
	����opencv��Լ����u,g�����ȴ洢
	����tag��16-31λ��0xFF00
	�轫���Ƶ���λ			cys	5/13/2008
*/
#define GIP_FMM_BAND	0x1000
#define GIP_FMM_KNOWN	0x2000
#define GIP_FMM_SIDE_0	0x000		//δ��
#define GIP_FMM_SIDE_1	0x100		//ĳһ��
#define GIP_FMM_SIDE_2	0x200		//����һ��

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

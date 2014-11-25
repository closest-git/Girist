#ifndef _GEO_ACTIVE_CONTOUR_H_
#define _GEO_ACTIVE_CONTOUR_H_

/*
	for flag of GAC_SOLVER
	依照opencv的约定，u,g行优先存储

*/
#define GAC_SOLVER_CONVERGE 0x10

/*
	for tat of GAC_SOLVER	
	RANGE:	0xFF0000
*/
#define TAG_INIT_FIX		 0xFF			//tag的0-15位始自init之后，就不再变化
#define GIP_GAC_INTERFACE	0x10000
#define GIP_GAC_SMOOTH		0x20000

typedef struct GAC_SOLVER_tag GAC_SOLVER;
struct GAC_SOLVER_tag{
	int M,N,ldu,shift,type,*tag,flag;
	double *g,*gx,*gy,g_max,gx_max,gy_max;
	double *u_buffer;//*u;
};

#ifdef __cplusplus
extern "C"  {
#endif
	void GAC_basic( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName );

#ifdef __cplusplus
}
#endif

#endif

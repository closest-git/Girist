#ifndef _GEO_ACTIVE_CONTOUR_H_
#define _GEO_ACTIVE_CONTOUR_H_
/*
	Solver for Chan¨CVese(CV) model	
	region-based models implemented via level-set techniques 

	v0.1	5/12/2008
		cys
*/
typedef struct CVM_SOLVER_tag CVM_SOLVER;
struct CVM_SOLVER_tag{
	int M,N,ldu,shift,type,*tag,flag;
	double *g,*gx,*gy,g_max,gx_max,gy_max;
	double *u_buffer;//*u;
};

#ifdef __cplusplus
extern "C"  {
#endif
	void CVM_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName );

#ifdef __cplusplus
}
#endif

#endif

#ifndef _GIP_CONTOUR_H_
#define _GIP_CONTOUR_H_

typedef struct GIP_CONTOUR_tag GIP_CONTOUR;
struct GIP_CONTOUR_tag{
	int M,N,nz;			//行数，列数
	int *ptr,*ind;
	double *val;
};

#ifdef __cplusplus
extern "C"  {
#endif


#ifdef __cplusplus
}
#endif

#endif
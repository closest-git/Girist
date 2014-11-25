#ifndef _GIA_HANDLE_H_
#define _GIA_HANDLE_H_

/*
	handle of IMAGE ANALYSIS 
	包含各种分析工具

*/
#include "gia_basic.h"
#include "gia_region.h"

/*
	GIa_Handle还提供buffer(I_buffer,D_buffer)
*/
typedef struct GIa_Handle_tag GIa_HANDLE;
struct GIa_Handle_tag{
	int flag,nRegion;
	int *I_buffer;
	GIa_HISTOGRAM histo;
	GIa_REGION *regions;
	double area,contrast,mean,deviation,rgn_thrsh;
	double *D_buffer;
};

#ifdef __cplusplus
extern "C"  {
#endif
	void GIa_HANDLE_init( GIa_HANDLE* hIa,GIP_OBJECT *,int nz_X,int flag);
	void GIa_HANDLE_clear( GIa_HANDLE* hIa );
#ifdef __cplusplus
}
#endif

#endif

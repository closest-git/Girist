#ifndef _GIA_QUERY_H_
#define _GIA_QUERY_H_

#include "gia_handle.h"
#include "gia_region.h"

/*
*/

#ifdef __cplusplus
extern "C"  {
#endif
	int query_iris_pupil( GIa_HANDLE *hIa,GIa_REGION **hRegion,double thresh,int flag );
#ifdef __cplusplus
}
#endif

#endif

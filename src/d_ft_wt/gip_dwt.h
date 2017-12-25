#ifndef _GIP_DWT_H_
#define _GIP_DWT_H_

#include "../pch/gip_def.h"

enum{
	GIP_CONVOLV_ROW=0x10000,GIP_CONVOLV_COL=0x20000
};

typedef enum{
	GIP_DWT_HAAR=0x10,
	GIP_DWT_DB1=0x20,
	GIP_DWT_DB4=0x100,GIP_DWT_SYM=0x101,
}GIP_DWT_FILTER;

typedef struct GIP_DWT_tag GIP_DWT;
struct GIP_DWT_tag{
	GIP_OBJECT *hObj;
	GIP_DWT_FILTER filter;
	int level,nMat,max_mat,M,N;
	GIP_MATRIX *arrMat,*hp,*lp,tempMat;
};

#ifdef __cplusplus
extern "C"  {
#endif
	void GIP_DWT_clear( GIP_DWT *hDWT );
	void GIP_DWT_init( GIP_DWT *hDWT,GIP_OBJECT *hObj,int level,GIP_DWT_FILTER filter,int flag );

	void GIP_Convol_testing( );
	void GIP_Convol_( GIP_MATRIX *mC,GIP_MATRIX *mA,GIP_MATRIX *mB,int flag );

	void GIP_dwt( GIP_OBJECT *hObj,double *temp,int level,GIP_DWT_FILTER filter,int flag );
	void GIP_filter_ft( int M,int N,GIP_MATRIX *mC,GIP_MATRIX *mA,GIP_MATRIX *mB,int flag );

	void GIP_fourier_matching( int M,int N,GIP_MATRIX *mA,GIP_MATRIX *mB,int flag );
	void GIP_fourier_spectrum( int M,int N,GIP_MATRIX *mA,double *mB,int flag );

#ifdef __cplusplus
}
#endif

#endif
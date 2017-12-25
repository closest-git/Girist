#ifndef _GIP_UTIL_H_
#define _GIP_UTIL_H_

/*
	命名约定
	v0.1	2/21/2009

	GIP_*_*_	中间全部小写
*/
#include <malloc.h>

/*
	注意和GIP_MASK_3_3等的匹配
*/
#define G_MASK_3	0x100
#define G_MASK_5	0x200

typedef enum{
	G_MASK_3_3=0x100,					//3*3 operator
	G_MASK_PREWITT_R=0x100,G_MASK_PREWITT_C,G_MASK_SOBEL_R,G_MASK_SOBEL_C,G_MASK_SOBEL_45,G_MASK_SOBEL_135,	
	G_MASK_5_5=0x200,
	G_MASK_FIVE_R=G_MASK_5_5,G_MASK_FIVE_C,G_MASK_SPOT_R,G_MASK_SPOT_C,
}GIP_MASK_OPERATOR;

typedef struct{
	size_t alloc,old;
	int nAlloc,nFree;
	double peak;
//	void* temp[1000];
} GIP_HEAP_INFO;


#ifdef __cplusplus
extern "C"  {
#endif

void GIP_Error( int status, const TCHAR* head,const TCHAR* err_msg, const TCHAR* file_name, int line );

//extern GIP_UTIL_DLL GIP_HEAP_INFO heap_info;
extern GIP_HEAP_INFO heap_info;

size_t HEAP_diff( );	
void GIP_init_heapinfo( );
void GIP_exit_heapinfo( );
void* GIP_alloc( size_t size );
void* GIP_realloc( void* old,size_t size );
void* GIP_calloc( size_t num,size_t size );
void GIP_free( void* memblock  );
double GIP_peak_heap( );	

void grus_dump_init( void **hDump,char *sTitle );
void grus_dump_finish( void *hDump );

void GIP_BUFFER_Init( int maxI,int maxD );
void GIP_BUFFER_Free( );
int* GIP_BI_alloc( int len );
double* GIP_BD_alloc( int len );
void GIP_BI_free( int len );
void GIP_BD_free( int len );

void distri_ouput( int len,double *w,double* val,const TCHAR* sPath,int flag )	;
void D_ouput( int len,double* val,const TCHAR* sPath )	;
void D_input( int len,double* val,const TCHAR* sPath )	;
void I_ouput( int len,int* val,const char* sPath );
void I_input( int len,int* val,const char* sPath );
void Equation_output( int M,int N,int dim,int *ptr,int *ind,double *val,double *rhs,const char* sPath );

extern int grus_stp,*grus_signal;

void GIP_STAMP_Init( int order );
void GIP_STAMP_Clear( );
void GIP_STAMP_Next( int need );
void GIP_STAMP_Need( int need );
#define GIP_STAMP_Set( ele )	{ grus_signal[ele] = grus_stp;}
#define GIP_STAMP_Test( ele )	(grus_signal[ele] == grus_stp)
#define GIP_STAMP_Set_2( ele,stp )	{ grus_signal[ele] = stp;}
#define GIP_STAMP_Test_2( ele,stp )	(grus_signal[ele] == stp)

TCHAR *GIP_path_( TCHAR *sTitle,int no );		//ys

double  GIP_DOFF2( int dim,double *Va,double *Vb );
int  GIP_MAXOFF_POS( int dim,double *Va,double *Vb,int );

void GIP_interactive( int M,int N,double *u,int ldu,int no,int *t_step,int flag );

void GIP_MODEL_info( char *sModel,GIP_MODEL model,int flag );
void GIP_RES_quality( int M,int N,int ldu,GIP_FLOATING *u,GIP_FLOATING *u_0,double *rmse,double *snr,double*,double* )	;

void GIP_CRS2CCS( int nRow,int nCol,GIP_FLOATING *crs_data,int* crs_ptr,
			 int *crs_ind,GIP_FLOATING **data,int **colptr,int **rowind,int* temp );

int GIP_Neibors( int pos,int *neibors,int M,int N,int ldu,int flag );
int GIP_adjacent_2(int M,int N,int ldu,int *u,int r,int c,int *adj,GIP_ADJACENT_TYPE type  );

double GIP_LAPLAS_( int pos,int M,int N,int ldu,GIP_FLOATING *u,int flag );
double GIP_Mask_3_( int pos,int M,int N,int ldu,GIP_FLOATING *u,GIP_MASK_OPERATOR type,int flag );
double GIP_Mask_5_( int pos,int M,int N,int ldu,GIP_FLOATING *u,GIP_MASK_OPERATOR type,int flag );

double GIP_sub_pixel( int M,int N,int ldu,GIP_FLOATING *u,GIP_FLOATING r_x,GIP_FLOATING c_x,int flag );
void GIP_pixel_shift_( int r,int c,double dr,double dc,int len,int *r1,int *c1 );
void GIP_arr_top( GIP_FLOATING *arr,int len,int nTop,int *arrPos,int flag );

int G_ReadImage( TCHAR *sPath,GIP_OBJECT *hObj,int flag );
int G_SaveImage( TCHAR *sPath,GIP_OBJECT *hObj,int flag );

/*
void GIP_IMAGE_init( int M,int N,double *val,GIP_IMAGE *hImag );
void GIP_IMAGE_clear( GIP_IMAGE *hImag );
*/
/*
	以下包含gip_util_1.cpp的函数	
*/
int GIP_OBJ_resize( GIP_OBJECT *hObj_0,GIP_OBJECT *hObj,int type,int flag );	
int GIP_OBJ_resize_nn( GIP_OBJECT *hObj_0,GIP_OBJECT *hObj,int *I_buffer,int flag );

#ifdef __cplusplus
}
#endif

#endif
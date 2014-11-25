#ifndef _GIP_IMAG_TRANSFORM_H_
#define _GIP_IMAG_TRANSFORM_H_

enum{
	GIP_TRANS_ABS=0x10,
	GIP_TRANS_LAPLACIAN=0x100,GIP_TRANS_GRADIENT=0x101
};

typedef enum{
	GIP_SMOOTH_MEDIAN=0x10,	
	GIP_SMOOTH_GAUSSIAN=0x20,
	GIP_SMOOTH_LAPLACIAN=0x30,
	GIP_SMOOTH_ANISO=0x100,
	GIP_SMOOTH_16=0x1000
}GIP_SMOOTH_TYPE;

#ifdef __cplusplus
extern "C"  {
#endif
	void GIP_Binarize( GIP_OBJECT *hSrc,GIP_OBJECT *hDest,double thresh,int flag );	
	void GIP_OBJ_sub( GIP_OBJECT *hSrc,GIP_OBJECT *hDest,int r_0,int c_0,int dr,int dc,int flag );
	void GIP_OBJ_shift( GIP_OBJECT *hSrc,GIP_OBJECT *hDest,int dr,int dc,int flag );
	void GIP_OBJ_init( GIP_OBJECT *hObj,int M,int N,int flag );	
	GIP_OBJECT *GIP_OBJ_create( int M,int N,GIP_FLOATING *data,int flag );

	void GIP_MATRIX_downsample( GIP_MATRIX *hMat,int s_0,int s_1,int step,int flag );
	void GIP_MATRIX_symextend( GIP_MATRIX *hMat, GIP_MATRIX *hSrc,int pad,int *k_M,int *k_N,int flag );
	void GIP_MATRIX_fill( GIP_MATRIX *hMat,GIP_MATRIX *hFill,int m_0,int n_0,int flag )	;

	int GIP_lsu_( GIP_OBJECT *hObj,int I_1,int I_0,double g_t0,double g_t1,int *mask,double *g_I,double *g_rad,int flag );
	int GIP_canny_( GIP_OBJECT *hObj,int I_1,int I_0,double g_t0,double g_t1,int *mask,double *g_I,double *g_rad,int flag );
//	void GIP_canny( GIP_OBJECT *hObj,GIP_OBJECT *hBinary,int flag );

	int GIP_circle_pos( GIP_OBJECT *hObj,double c_r,double c_c,double c_R,int *arrPos,int flag );
	void GIP_Circle( GIP_OBJECT *hObj,int c_r,int c_c,double c_R,GIP_IMAG_TYPE type );
	void GIP_Parabolic( GIP_OBJECT *hObj,double a,double h,double k,GIP_IMAG_TYPE type );
	double GIP_Circle_GMeasure( GIP_OBJECT *hObj,int c_r,int c_c,double c_R,double xita, GIP_FLOATING *g_I,GIP_FLOATING*g_rad,int *mask,int flag );

	void GIP_EdgeMap( GIP_OBJECT *hObj,GIP_OBJECT *hEdge,double thrsh,int type,int flag );	
	int GIP_Edge2List( GIP_OBJECT *hEdge,int *E_ptr,int *E_pos,int *I_buffer,int flag );

	void GIP_OBJ_trans( GIP_OBJECT *hObj,int type,int flag );
	void GIP_Histo_Equal(  GIP_OBJECT *hObj,int L,int type,int *temp,int flag );
	void GIP_Contrast_Enhance(  GIP_OBJECT *hObj,int L,int type,int flag );
	void GIP_Smooth( GIP_OBJECT *hObj,GIP_SMOOTH_TYPE type,GIP_FLOATING *d_temp,double q,int flag );


#ifdef __cplusplus
}
#endif

#endif
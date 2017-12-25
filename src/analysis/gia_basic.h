#ifndef _GIA_BASIC_H_
#define _GIA_BASIC_H_

/*
	basic tool and algorithm for IMAGE ANALYSIS 

	1 Histogram processing
*/
typedef struct GIa_Histo_Entry_tag GIa_Histo_Entry;
struct GIa_Histo_Entry_tag{			//extened
	int level,nz,nBlob;			//blob is a 4-connected piece with same intensity
	double r,c;		
};

/*
	haar like rectangle feature
*/
typedef enum	{
	G_FETR_1_2,G_FETR_2_1,
	G_FETR_1_3,G_FETR_3_1,
	G_FETR_1_4,G_FETR_4_1,		
	G_FETR_2_2,
	G_FETR_POINT,
	G_FETR_T_12,G_FETR_T_21,
	G_FETR_T_13,G_FETR_T_31,
	G_FETR_T_14,G_FETR_T_41,	

	G_FETR_nTYPEs
}G_RECT_FETR;

typedef struct _weight_rectangle_tag _W_rect;
struct _weight_rectangle_tag{
	float r,c,height,width;
	float scale; 
	int type;
};
/*
	1 考虑到直方图变换的需要，[0,I_limit]可能大于[I_0,I_1]
	2 基于图像的不连续性，nPeak用处不大
*/
typedef struct GIa_Histogram_tag GIa_HISTOGRAM;
struct GIa_Histogram_tag{
	int type,I_limit,I_0,I_1,I_peak,nPeak;
	double r_peak;			//*h,
	GIa_Histo_Entry *entrys;
};

#ifdef __cplusplus
extern "C"  {
#endif
	void GIa_Histo_init( GIP_OBJECT *hObj,GIa_HISTOGRAM *histo,int type,int *temp,int flag );
	void GIa_Histo_clear( GIa_HISTOGRAM *histo );
	void GIa_Histo_attrib( GIa_HISTOGRAM *histo,double *attrib,int flag );

	void GIP_GlobalThresh_2( GIP_OBJECT *hObj,double *threshold,int flag );
	int GIP_GlobalThresh_x( GIP_OBJECT *hObj,double I_0,double I_1,int *nTh,double *threshold,int flag );

	double GIa_CHT( GIP_OBJECT *hObj,int nz,int* arrPos,double *R,double *c_r,double *c_c,int flag );

	void GIa_Circle_Detect_GPV( GIP_OBJECT *hObj,int th_1,int type,double *D_buffer,int *I_buffer,int flag );

	int G_rect2fetr( int *f_rect,int *pos,int *w,int ldx,int flag );
	int W_rect2fetr( _W_rect *w_rect,int nzRect,int *pos,double *w,int ldx,int flag );
	void _W_rect_set( _W_rect *hRect,float r,float c,float w,float h,float scale );

#ifdef __cplusplus
}
#endif

#endif
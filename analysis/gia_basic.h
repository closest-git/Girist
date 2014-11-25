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

#ifdef __cplusplus
}
#endif

#endif
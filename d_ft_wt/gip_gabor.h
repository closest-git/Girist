#ifndef _GIP_GABOR_H_
#define _GIP_GABOR_H_


enum{
	GIP_GABOR_FREQ=0x10,
	GIP_GABOR_1D=0x100,GIP_GABOR_2D=0x200,
};

typedef enum{
	GIP_GABOR_MALI=0x10,
	GIP_GABOR_FIELD=0x100,
}GIP_GABOR_TYPE;

/*
	为了方便其它函数调用，增加hMat(直接对应于M,N,filter)
*/
typedef struct GIP_GABOR_tag GIP_GABOR;
struct GIP_GABOR_tag{
	GIP_GABOR_TYPE type;
	int M,N,flag;
	double f,sigma_r,sigma_c,lenda;		//	;filter bandwidth - R;filter bandwidth - C;base wavelength
	GIP_FLOATING *filter/*,*re,*im*/;
	GIP_MATRIX *hMat;
};


typedef struct GIP_GABOR_BANK_tag GIP_GBANK;
struct GIP_GABOR_BANK_tag{
	int type;
	int nFilter;
//	double 
};

#ifdef __cplusplus
extern "C"  {
#endif
	void GIP_GABOR_clear( GIP_GABOR *hGabor );
	void GIP_GABOR_init( GIP_GABOR *hGabor,int M,int N,double dr,double dc,GIP_GABOR_TYPE type,int flag );	

	void GIP_LOG_GABOR_init( GIP_GABOR *hGabor,int M,int N,double ,double ,double ,GIP_GABOR_TYPE,int flag );

	void GIP_GABOR_filter( GIP_GABOR *hGabor,int M,int N,GIP_FLOATING *data,double *re,double *im,int flag );	

#ifdef __cplusplus
}
#endif

#endif
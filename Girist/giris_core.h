#ifndef _GIRIS_CORE_H_
#define _GIRIS_CORE_H_

#include "../analysis/gia_handle.h"
/*		
*******************缺省配置，适用于MMU*********************		
	v0.1	
		4/1/2009		cys
	1	无需 GIP_M_HISTO
	2	n_limit=可取为9000~15000，适当增大有利于提高精度
		对应n_limit的q可取为0.5~0.8
	3   control[GIRIS_T_INTN_1]=0.6，放大到0.8似乎没差别
	control parameters:
		control[GIRIS_MODE_]=1.0;
		control[GIRIS_T_GRAD]=5.0;
		control[GIRIS_T_INTN_0]=0.0;
		control[GIRIS_T_INTN_1]=0.8;
		control[GIRIS_T_MTCH_S]=0.25;
		control[GIRIS_T_MTCH_P]=0.6;
		control[GIRIS_FV_SIFT]=6;
		control[GIRIS_LG_WAVELEN]=32;
		control[GIRIS_GRAD_BIAS]=PI/18;
		control[GIRIS_SPOT_AREA]=0.03;
		control[GIRIS_T_R0]=0.05;
		control[GIRIS_SMOOTH_Q]=1.0;
*******************缺省配置，适用于MMU*********************	
*/

/*		
*******************适用于UBIRIS*********************		
	v0.1	
		4/1/2009		cys
	1	无需 GIP_M_HISTO
	2	n_limit=可取为9000~15000，适当增大有利于提高精度
	3   有特小pupil，sclera/pupil 放大到6.5
	4	_IRIS_normalize_包含smooth后，略有提升
*******************适用于UBIRIS*********************	
*/

/*		
*******************适用于CASIA*********************		
	v0.1	
		4/1/2009		cys
	1	必须采取多尺度架构 n_limit=640*480时，定位速度要慢10倍
*******************适用于CASIA*********************	
*/

/*GIRIS_MODE_ */
#define GIRIS_MODE_TRAINING			1		
#define GIRIS_MODE_VERIFY			2
#define GIRIS_MODE_MATCHING			3

/* 
	global error code
*/
#define GIRIST_WARN_INNER_BOUNDARY (-10001)		//Failed to verify the inner boundary of iris
#define GIRIST_WARN_OUTER_BOUNDARY (-10002)		//Failed to verify the outer boundary of iris
//#define GIRIST_ERROR_LOCATION (-10100)			//Failed to find localization of iris
#define GIRIST_ERROR_BOUNDARY (-10101)				//
#define GIRIST_ERROR_INIT_CONTROL (-10004)				//

enum{	
	IRIS_CENTER_FIX_=0x10000 
};
typedef enum	{
	IRIS_UPPER_EYELID=0x10,IRIS_LOWER_EYELID=0x11,
	IRIS_PUPIL_BOUNDARY=0x20,IRIS_SCLERA_BOUNDARY=0x21,
	IRIS_NOISE_EYELID=0x100,
}G_IRIS_REGION;

/*
	mask标记每个点的状态
	G_IRIS_location：mask[pos]!=1,则该点不用考虑
*/
typedef struct G_IRIS_tag G_IRIS_;
struct G_IRIS_tag{
	GIa_HANDLE *hIa;
	GIP_OBJECT *hObj,normal;//,*hMask;
	GIP_FLOATING *g_I,*g_rad;
	int shift,type,flag,nzNoise,*mask;
//	double lenda,alpha,beta;
	float p_radius,p_r,p_c,s_radius,s_r,s_c;			//the center of the pupil-boundary and sclera-boundary
	float p_match,s_match,ecc;							//匹配率,偏心率
	double ul_a,ul_h,ul_k;							//upper eyelids,lower eyelids
	double x;					
};


#ifdef __cplusplus
extern "C"  {
#endif
//	#define GIRIS_CONTROL_ITEM 0x20
//	extern double GIRIS_CONTROL[GIRIS_CONTROL_ITEM];

	int G_IRIS_init_control( TCHAR *sControlPath,double *control,int flag );

	//void G_IRIS_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_0,int model,int flag,char* sWndName );
	void G_IRIS_program( char *sTrainPath,char *sLibPath,int flag );
	void G_IRIS_matching( char *sFilePath,char *sLibPath,int flag );
//以下由界面调用,与界面有信息交互
	int G_IRIS_location( GIP_LIB *hLib,int flag );
	int G_IRIS_coding( GIP_LIB *hLib,int flag );

//	G_IRIS_ *G_IRIS_process( double *control,TCHAR *sFilePath,int model,int flag,int *ret );	
	G_IRIS_ *G_IRIS_on_entry( double *control,GIP_LIB_entry *hEntry,int model,int flag,int *ret );	
	void G_IRIS_clear( G_IRIS_ *hIris )	;
	int G_IRIS_LIB_info( GIP_LIB *hLib,TCHAR* sFilePath,int flag );	

#ifdef __cplusplus
}
#endif

#endif
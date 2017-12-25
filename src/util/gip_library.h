#ifndef _GIP_LIBRARY_H_
#define _GIP_LIBRARY_H_

#include <tchar.h>
#define MAX_ENTRY_PATH 100
 
typedef enum{
	G_ENTRY_UNKNOWN=0,
	G_ENTRY_OK=1,		
	G_ENTRY_WARNING=2,			//结果可能不对
	G_ENTRY_FAILED=-1,			//处理失败

	G_ENTRY_SELECT=0x10000,			//被选中

	G_ENTRY_ALL=0xFFFF,
}G_ENTRY_STATUS;	

typedef enum{
	G_LIB_UNKNOWN=0,
	G_LIB_UNTITLE=0x100,		//首次创建,并无对应的文件
	G_LIB_UPDATE=0x200,		

}G_LIB_STATUS;

/*
	G_LIB
		每个类的元素连续存储,以提高计算效率
	
	1 支持部分采样，
	采样空间:	即nEntry,nClass,arrEntry）对应于所有的元素
	W空间:	记录成功提取特征的元素,特征向量存储于fv_dat.与采样空间的对应关系为W_dim,*W_ptr,*W_no
*/
typedef struct GIP_LIB_entry_tag GIP_LIB_entry;
struct GIP_LIB_entry_tag{
	TCHAR sFilePath[MAX_ENTRY_PATH]/**sTitle*/;		//对应的文件
	int no,cls;		//no指在class中的编号
//	int model,nPathLen;
	float p_radius,p_r,p_c,s_radius,s_r,s_c,p_match,s_match,x;
	int status;
};

/*
	1 isUpdate,status由界面程序调用
	2 G_FVector_独立于GIP_LIB,GIP_LIB只引用hfv,不负责释放hfv
	3 需将GIP_LIB_entry::sFilePath提出，单独存储

	3/11/2009	增加(D_inter,D_intra,D_span)，用于Distance Distributions
*/
typedef struct GIP_LIB_tag GIP_LIB;	
struct GIP_LIB_tag{
	int version;					//版本号

	double user_control[GIP_USER_CONTROL];	 		//记录生成lib时的控制参数
	TCHAR sPath[MAX_PATH_LEN],sRoot[MAX_PATH_LEN];		//path,由界面程序设置,并不存储到文件中
	TCHAR *sNote;					//带注释,并存储到文件中
	int lenNote,isUpdate,status;			
	G_FVector_ *hfv;
	int nEntry,nSample,nMaxEntry,nClass,fv_n_0,fv_n,D_span;
	double CRR;						//recognition index 
//	double sep;						//recognition index 
	double L_spd,M_spd;					//Location, matching speed
//	int W_dim,*W_ptr,*W_no;
	float *D_inter,*D_intra;		//Normalized Intra-Class and Inter-Class Distance
	G_U8 *fv_dat;
	//	double *W_proj;		//投影矩阵
	GIP_LIB_entry *arrEntry;
};

#ifdef __cplusplus
extern "C"  {
#endif
	int GIP_LIB_init( GIP_LIB *hLIB,TCHAR *sLibFile,G_FVector_ *hfv,double *control,int );
	void GIP_LIB_clear( GIP_LIB *hLIB );
	int GIP_LIB_save( GIP_LIB *hLIB,TCHAR*,int flag );
	int GIP_LIB_load( GIP_LIB *hLIB,G_FVector_ *hfv,TCHAR*,int flag );
	void GIP_LIB_Append( GIP_LIB *hLib,TCHAR *sPath,int no,int cls,int flag );

	void GIP_LIB_reducing( GIP_LIB *,GIP_FLOATING *, int );
	void GIP_LIB_update( GIP_LIB *hLIB,char*,int flag );
	int GIP_LIB_extract( GIP_LIB *hLib,GIP_LIB *hLib_0,TCHAR *sPath,G_FVector_ *hfv_0,int *mark,int flag );

//	int GIP_LIB_scan( GIP_LIB *hLib,int flag )	;
	int GIP_LIB_discriminate( GIP_LIB *hLib,G_DCRIMI_ *hDcrimi,double *dist_4,int flag );
	void GIP_LIB_sort_( GIP_LIB *hLib,double *s,int *temp,int flag );

#ifdef __cplusplus
}
#endif

#endif
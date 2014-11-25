#ifndef _GIP_LIBRARY_H_
#define _GIP_LIBRARY_H_

#include <tchar.h>

typedef enum{
	G_FILE_UNKNOWN=0,
	G_IMAGE_BEGIN=0x10,
	G_IMAGE_JPG,G_IMAGE_BMP,
	G_IMAGE_END,
	G_FILE_DIRECTORY=0x100
}G_FILE_TYPE;
 
typedef enum{
	G_ENTRY_UNKNOWN=0,
	G_ENTRY_OK=1,		
	G_ENTRY_WARNING=2,			//������ܲ���
	G_ENTRY_FAILED=-1,			//����ʧ��

	G_ENTRY_ALL=10000,
}G_ENTRY_STATUS;	

typedef enum{
	G_LIB_UNKNOWN=0,
	G_LIB_UNTITLE=0x100,		//�״δ���,���޶�Ӧ���ļ�
	G_LIB_UPDATE=0x200,		

}G_LIB_STATUS;

/*
	G_LIB
		ÿ�����Ԫ�������洢,����߼���Ч��
	
	1 ֧�ֲ��ֲ�����
	�����ռ�:	��nEntry,nClass,arrEntry����Ӧ�����е�Ԫ��
	W�ռ�:	��¼�ɹ���ȡ������Ԫ��,���������洢��fv_dat.������ռ�Ķ�Ӧ��ϵΪW_dim,*W_ptr,*W_no
*/
typedef struct GIP_LIB_entry_tag GIP_LIB_entry;
struct GIP_LIB_entry_tag{
	TCHAR sFilePath[MAX_PATH_LEN]/**sTitle*/;		//��Ӧ���ļ�
	int no,cls,nPathLen;		//noָ��class�еı��
	int model;
	float p_radius,p_r,p_c,s_radius,s_r,s_c,p_match,s_match,x;
	G_ENTRY_STATUS status;
};

/*
	1 isUpdate,status�ɽ���������
	2 G_FVector_������GIP_LIB,GIP_LIBֻ����hfv,�������ͷ�hfv

	3/11/2009	����(D_inter,D_intra,D_span)������Distance Distributions
*/
typedef struct GIP_LIB_tag GIP_LIB;	
struct GIP_LIB_tag{
	int version;					//�汾��

	double user_control[GIP_USER_CONTROL];			//��¼����libʱ�Ŀ��Ʋ���
	TCHAR sPath[MAX_PATH_LEN];		//path,�ɽ����������,�����洢���ļ���
	TCHAR *sNote;					//��ע��,���洢���ļ���
	int lenNote,isUpdate,status;			
	G_FVector_ *hfv;
	int nEntry,nSample,nMaxEntry,nClass,fv_n_0,fv_n,D_span;
	double CRR;						//recognition index 
	double L_spd;					//Location and normalize speed
//	int W_dim,*W_ptr,*W_no;
	float *D_inter,*D_intra;		//Normalized Intra-Class and Inter-Class Distance
//	G_U8 *fv_dat;
	G_U8 *fv_dat;
	//	double *W_proj;		//ͶӰ����
	GIP_LIB_entry *arrEntry;
};

#ifdef __cplusplus
extern "C"  {
#endif
	int GIP_LIB_init( GIP_LIB *hLIB,TCHAR *sLibFile,G_FVector_ *hfv,double *control,int );
	void GIP_LIB_clear( GIP_LIB *hLIB );
	int GIP_LIB_save( GIP_LIB *hLIB,TCHAR*,int flag );
	int GIP_LIB_load( GIP_LIB *hLIB,G_FVector_ *hfv,TCHAR*,int flag );

	void GIP_LIB_reducing( GIP_LIB *,GIP_FLOATING *, int );
	void GIP_LIB_update( GIP_LIB *hLIB,char*,int flag );
	int GIP_LIB_extract( GIP_LIB *hLib,GIP_LIB *hLib_0,TCHAR *sPath,G_FVector_ *hfv_0,int *mark,int flag );

//	int GIP_LIB_scan( GIP_LIB *hLib,int flag )	;
	int GIP_LIB_discriminate( GIP_LIB *hLib,G_DCRIMI_ *hDcrimi,int flag );
	void GIP_LIB_sort_( GIP_LIB *hLib,double *s,int *temp,int flag );

#ifdef __cplusplus
}
#endif

#endif
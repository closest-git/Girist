#ifndef _GIP_LIBRARY_H_
#define _GIP_LIBRARY_H_

typedef enum{
	G_FILE_UNKNOWN=0,
	G_IMAGE_BEGIN=0x10,
	G_IMAGE_JPG,G_IMAGE_BMP,
	G_IMAGE_END,
	G_FILE_DIRECTORY=0x100
}G_FILE_TYPE;
/*
	G_LIB
		ÿ�����Ԫ�������洢,����߼���Ч��

	1 ֧�ֲ��ֲ�����
	�����ռ�:	��nEntry,nClass,arrEntry����Ӧ�����е�Ԫ��
	W�ռ�:	��¼�ɹ���ȡ������Ԫ��,���������洢��W_fv.������ռ�Ķ�Ӧ��ϵΪW_dim,*W_ptr,*W_no
*/
typedef struct GIP_LIB_entry_tag GIP_LIB_entry;
struct GIP_LIB_entry_tag{
	char sFilePath[MAX_PATH_LEN]/**sTitle*/;		//��Ӧ���ļ�
	int no,cls;		//noָ��class�еı��
	int model,status;
};

typedef struct GIP_LIB_tag GIP_LIB;
struct GIP_LIB_tag{
	G_FVector_ *hfv;
	int nEntry,nSample,nMaxEntry,nClass,fv_n_0,fv_n;
	int W_dim,*W_ptr,*W_no;
	char *W_fv;
//	double *W_proj;		//ͶӰ����
	GIP_LIB_entry *arrEntry;
};

#ifdef __cplusplus
extern "C"  {
#endif
	int GIP_LIB_init( GIP_LIB *hLIB,char *sLibFile,int );
	void GIP_LIB_clear( GIP_LIB *hLIB );
	void GIP_LIB_save( GIP_LIB *hLIB,char*,int flag );
	void GIP_LIB_load( GIP_LIB *hLIB,char*,int flag );

	void GIP_LIB_reducing( GIP_LIB *,GIP_FLOATING *, int );
	void GIP_LIB_update( GIP_LIB *hLIB,char*,int flag );

#ifdef __cplusplus
}
#endif

#endif
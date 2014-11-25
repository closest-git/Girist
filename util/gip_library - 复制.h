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
		每个类的元素连续存储,以提高计算效率

	1 支持部分采样，
	采样空间:	即nEntry,nClass,arrEntry）对应于所有的元素
	W空间:	记录成功提取特征的元素,特征向量存储于W_fv.与采样空间的对应关系为W_dim,*W_ptr,*W_no
*/
typedef struct GIP_LIB_entry_tag GIP_LIB_entry;
struct GIP_LIB_entry_tag{
	char sFilePath[MAX_PATH_LEN]/**sTitle*/;		//对应的文件
	int no,cls;		//no指在class中的编号
	int model,status;
};

typedef struct GIP_LIB_tag GIP_LIB;
struct GIP_LIB_tag{
	G_FVector_ *hfv;
	int nEntry,nSample,nMaxEntry,nClass,fv_n_0,fv_n;
	int W_dim,*W_ptr,*W_no;
	char *W_fv;
//	double *W_proj;		//投影矩阵
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
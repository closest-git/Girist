#ifndef _GIP_PACKAGE_H_
#define _GIP_PACKAGE_H_

#include "cv.h"
#include <stdio.h>

typedef IplImage GIP_IMAGE;

#define GIP_OPENCV_CREATEWND	0x10000

#ifdef __cplusplus
extern "C" {
#endif
	GIP_OBJECT* GIP_IMG2OBJ( GIP_IMAGE *hIMG,int flag );

	GIP_IMAGE *GIP_LoadImage( TCHAR *tFilePath,int iscolor,int flag );
	void GIP_import_IplImage( IplImage *hIMG,int *M,int *N,GIP_FLOATING **hval,int flag );
	void GIP_export_IplImage_d( IplImage *hIMG,double *data,int,int );
	void GIP_show_IplImage_d( int M,int N,char *sTitle,double *u,int ldu,int flag );
	void GIP_save_IplImage_d( int height,int width,TCHAR *sTitle,double *u,int ldu,int flag );
	void GIP_save_IplImage_I8( int height,int width,char *sTitle,char *u,int ldu,int flag );

	int GIP_IMG_Noise( IplImage* hIMG,double mean,double delta );
	double GIP_IMG_PSNR( int height,int width,double *u,int ldu,GIP_IMAGE *hIMG,double *temp );

	double GIP_IMG_delta( GIP_IMAGE *hIMG );
	void GIP_IMG_segm( int height,int width,double *u,int ldu,GIP_FLOATING *u_0,int flag  );
	GIP_IMAGE *GIP_IMG_resize( GIP_IMAGE *hIMG_0,int height,int width,int flag );
	void GIP_IMG_v( IplImage*f_image,IplImage* hImg,int alpha,int flag );

//for Berkeley Segmentation Dataset
	int fReadLine( FILE* fp,TCHAR sLine[MAX_LINE] );
	void GIP_import_BSDSseg( IplImage *hIMG,char *sSegFile,int flag );

#ifdef __cplusplus
}
#endif

#endif
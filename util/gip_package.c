/*
1	集中了对OPENCV的调用
	注意TCHAR和char之间的转换
*/
#include <malloc.h>
#include <memory.h>
#include <limits.h>
#include <FLOAT.h>
#include <stdio.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/GIP_util.h"
#include "cxerror.h"
#include "highgui.h"

static int g_nSave=0;

/*
	ISSUE-NEEDED:
		1 替换为wcstombs_s
		2 采用更稳定的函数

	注意:
		1 buf固定长度为MAX_PATH,由调用程序分配

	把wchar_t*转换称char*   
	http://topic.csdn.net/t/20050117/20/3733034.html

	v0.1	cys
		3/5/2009	
*/
char *WSTRToAnsi( wchar_t*   Msg,char *buf ){		
	int   len   =   wcstombs(NULL,   Msg,   0),ret;   
//	char*   buf   =   GIP_alloc( sizeof(char)*(len+1) ); 

	ASSERT( len<MAX_PATH );
	ret=wcstombs(buf,   Msg,   len);   
	buf[len]   =   0;   

	return   buf;   
}

/*
	注意：
		1 包含必要的字符转换

	Copyright 2008-present, Grusoft

	v0.1	cys
		3/5/2009	
*/
GIP_IMAGE *GIP_LoadImage( TCHAR *tFilePath,int iscolor,int flag )	{
	GIP_IMAGE *hIMG=GIP_NULL;
	char sPath[MAX_PATH];
	
	WSTRToAnsi( tFilePath,sPath );
	hIMG=cvLoadImage( sPath,iscolor);
	if( hIMG ==0x0 )
		GIP_ERROR( "GIP_LoadImage",-1, "the image must be a single-channel byte image" );

GIP_exit:
	return hIMG;
}

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		3/28/2008	
*/
GIP_OBJECT* GIP_IMG2OBJ( GIP_IMAGE *hIMG,int flag )	{
	int height=hIMG->height,width=hIMG->width,step=hIMG->widthStep/sizeof(uchar),i,j;
	GIP_OBJECT* hObj;
	GIP_FLOATING *val;
	uchar* data = (uchar *)hIMG->imageData;;

	hObj = GIP_OBJ_create( height,width,GIP_NULL,0x0 );
	val = hObj->data;
	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			*val = data[i*step+j];
			val++;
		}
	}

	return hObj;
}

/*
	暂定为import jpg

	注意：
	1、IplImage仅出现于此函数中
	2、除非flag指定，val缺省ldu=N

	Copyright 2008-present, Grusoft

	v0.1	cys
		3/28/2008	
*/
void GIP_import_IplImage( IplImage *hIMG,int *M,int *N,GIP_FLOATING **hval,int flag )	{
	int height=hIMG->height,width=hIMG->width,i,j;
	int channels=hIMG->nChannels,step=hIMG->widthStep/sizeof(uchar);
	GIP_FLOATING *val;
	uchar* data = (uchar *)hIMG->imageData;;

	*M=-1;	*N=-1;
	*hval=NULL;
	if( hIMG->depth != IPL_DEPTH_8U || hIMG->nChannels !=1 )
		GIP_ERROR( "GIP_import_IplImage",-1, "the image must be a single-channel byte image" );

	val=(GIP_FLOATING*)GIP_alloc( sizeof(GIP_FLOATING)*height*width );
	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			*val = data[i*step+j];
			val++;
		}
	}
	val -= height*width;
//	 = 111;	;

	*M=height;	*N=width;
	*hval=val;
GIP_exit:
	;
}


/*
	将数据输出到IplImage

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/9/2008	
*/
void GIP_export_IplImage_d( IplImage *hIMG,double *data,int lda,int flag )	{
//	CvArr *output_img;
//	CvMat output_hdr; 
	int i,j,out_val;
	int height=hIMG->height,width=hIMG->width,step=hIMG->widthStep/sizeof(uchar);
	double a,d_min=DBL_MAX,d_max=0.0;
	uchar* i_data = (uchar *)hIMG->imageData;
	if( hIMG->depth != IPL_DEPTH_8U || hIMG->nChannels !=1 )
		GIP_ERROR( "GIP_export_IplImage_d",-1, "the image must be a single-channel byte image" );
//	output_img = cvGetMat( hIMG, &output_hdr,NULL,0 );
//	CV_CALL( output_img = cvGetMat( hIMG, &output_hdr,NULL,0 ));
//    if( !CV_ARE_SIZES_EQ(hIMG,output_img) )
//        GIP_ERROR( "GIP_export_IplImage_d",CV_StsUnmatchedSizes, "All the input and output images must have the same size" );

	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			a = data[i*lda+j];
			d_min=MIN(d_min,a);	d_max=MAX(d_max,a);		
		}
	}
	ASSERT( d_max!=d_min );

	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
//			a = (data[i*lda+j]-d_min)/(d_max-d_min);
//			out_val = cvRound(a*UCHAR_MAX);
			out_val = cvRound(data[i*lda+j]);
			out_val = out_val<0 ? 0 : (out_val > UCHAR_MAX ? UCHAR_MAX : out_val);
//			ASSERT( out_val>=0 && out_val <= UCHAR_MAX );
			i_data[i*step+j] = out_val;	
//            CV_MAT_ELEM(*out,uchar,i-1,j-1) = CV_CAST_8U(out_val);		
		}
	}

GIP_exit:
		;
}

/*
	将数据输出到IplImage

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/21/2008	
*/
void GIP_show_IplImage_d( int height,int width,char *sTitle,double *u,int ldu,int flag )	{
	IplImage* hIMG;
	int i,j,out_val,step;
	double a,d_min=DBL_MAX,d_max=0.0;
	uchar* i_data;

	if( (hIMG=cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,1))==0x0 )
		GIP_ERROR( "GIP_show_IplImage_d",-1, "the image must be a single-channel byte image" );

	ASSERT( height==hIMG->height && width==hIMG->width );
	step=hIMG->widthStep/sizeof(uchar);
	i_data = (uchar *)hIMG->imageData;
	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			a = u[i*ldu+j];
			d_min=MIN(d_min,a);	d_max=MAX(d_max,a);		
		}
	}
	if( d_max==d_min )	{
		out_val = cvRound(d_max);			ASSERT( out_val>=0 && out_val <= UCHAR_MAX );
		for( i = 0; i < height; i++ )	{
			for( j = 0; j < width; j++ )	{
				i_data[i*step+j] = out_val;	
			}
		}
	}else	{
		for( i = 0; i < height; i++ )	{
			for( j = 0; j < width; j++ )	{
				a = (u[i*ldu+j]-d_min)/(d_max-d_min);
				out_val = cvRound(a*UCHAR_MAX);			ASSERT( out_val>=0 && out_val <= UCHAR_MAX );
				i_data[i*step+j] = out_val;	
			}
		}
	}
/*
	if( BIT_TEST( flag,GIP_OPENCV_CREATEWND ) )
 		cvNamedWindow( sTitle, 1 );
	cvShowImage( sTitle, hIMG );
*/
	cvReleaseImage( &hIMG );
GIP_exit:
		;
}

/*
	调用GIP_save_IplImage_d
	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/3/2009	
*/
void GIP_save_IplImage_I8( int height,int width,TCHAR *sTitle,char *u,int ldu,int flag )	{
	int i;
	double *u_d;

	u_d=GIP_alloc( sizeof(double)*height*ldu );
	for( i = 0; i < height*ldu; i++ )	{
		u_d[i] = u[i];
	}

	GIP_save_IplImage_d( height,width,sTitle,u_d,ldu,flag );

	GIP_free( u_d );
}

/*
	将数据输出到IplImage

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/21/2008	
*/
void GIP_save_IplImage_d( int height,int width,TCHAR *tFilePath,double *u,int ldu,int flag )	{
	IplImage* hIMG;
	int i,j,out_val,step;
	double a,d_min=DBL_MAX,d_max=0.0;
	uchar* i_data;
	char *sPath=GIP_alloc( sizeof(char)*MAX_PATH );
	
	g_nSave++;
//	G_PRINTF( "{GIP_save_IplImage_d} is called to save file <<<%s>>> \r\n",sTitle  );
	if( (hIMG=cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,1))==0x0 )
		GIP_ERROR( _T("GIP_save_IplImage_d"),-1, _T("the image must be a single-channel byte image") );

	ASSERT( height==hIMG->height && width==hIMG->width );
	step=hIMG->widthStep/sizeof(uchar);
	i_data = (uchar *)hIMG->imageData;
	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			a = u[i*ldu+j];
			d_min=MIN(d_min,a);	d_max=MAX(d_max,a);		
		}
	}
	if( d_max==d_min )		{
		out_val = (uchar)(d_min);			
		for( i = 0; i < height; i++ )	{
			for( j = 0; j < width; j++ )	
				i_data[i*step+j] = out_val;	
		}	
	}else	{
		for( i = 0; i < height; i++ )	{
			for( j = 0; j < width; j++ )	{
				a = (u[i*ldu+j]-d_min)/(d_max-d_min);
				out_val = cvRound(a*UCHAR_MAX);			ASSERT( out_val>=0 && out_val <= UCHAR_MAX );
				i_data[i*step+j] = out_val;	
			}
		}
	}

	WSTRToAnsi( tFilePath,sPath );
#ifdef _DEBUG
	if( !cvSaveImage(sPath,hIMG) ) 		
		GIP_ERROR( _T("GIP_save_IplImage_d"),-1, _T("Could not save") );
#endif

	cvReleaseImage( &hIMG );
GIP_exit:
	GIP_free( sPath );
}

/*
	将数据输出到IplImage
	注意:
		cvReleaseImage由调用函数负责

	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/4/2008	
*/
IplImage* GIP_set_IplImage_d( int height,int width,double *u,int ldu,int flag )	{
	IplImage* hIMG=GIP_NULL;
	int i,j,out_val,step;
	double a,d_min=DBL_MAX,d_max=0.0;
	uchar* i_data;

	if( (hIMG=cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,1))==0x0 )
		GIP_ERROR( "GIP_set_IplImage_d",-1, "the image must be a single-channel byte image" );

	ASSERT( height==hIMG->height && width==hIMG->width );
	step=hIMG->widthStep/sizeof(uchar);
	i_data = (uchar *)hIMG->imageData;
	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			a = u[i*ldu+j];
			d_min=MIN(d_min,a);	d_max=MAX(d_max,a);		
		}
	}
	if( d_max==d_min )		{
		out_val = (uchar)(d_min);			
		for( i = 0; i < height; i++ )	{
			for( j = 0; j < width; j++ )	
				i_data[i*step+j] = out_val;	
		}	
	}else	{
		for( i = 0; i < height; i++ )	{
			for( j = 0; j < width; j++ )	{
				a = (u[i*ldu+j]-d_min)/(d_max-d_min);
				out_val = cvRound(a*UCHAR_MAX);			ASSERT( out_val>=0 && out_val <= UCHAR_MAX );
				i_data[i*step+j] = out_val;	
			}
		}
	}

GIP_exit:
	return hIMG	;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		5/26/2008	
*/
int GIP_IMG_Noise( IplImage* hIMG,double mean,double delta )	{
	int i,j,height=hIMG->height,width=hIMG->width,pos;
	int	step=hIMG->widthStep/sizeof(uchar),a,i_old;
	CvMat* M_rand = cvCreateMat(height,MAX(step,width),CV_32FC1);	
	uchar* i_data = (uchar *)hIMG->imageData;
	float *noise=M_rand->data.fl,n_max,n_min,n_avg=0.0,n2_avg=0.0;

	CvRNG rng = cvRNG(0xffffffff);
//	cvRandArr( &rng,M_rand,CV_RAND_NORMAL,cvRealScalar(mean),cvRealScalar(delta) );
//	if( !cvSaveImage("noise.bmp",rng)) 		
//		GIP_ERROR( "main",-1, "Could not save" );

	n_max=FLT_MIN,			n_min=FLT_MAX;
	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			pos = i*step+j;
			i_old = i_data[pos];
			a = (int)(i_data[pos] + noise[pos]+0.5);
			n_max=MAX(n_max,noise[pos]);
			n_min=MIN(n_min,noise[pos]);
			i_data[pos] = a < 0 ? 0 : a > UCHAR_MAX ? UCHAR_MAX : a ;
			n_avg += i_old-i_data[pos];
			n2_avg += (i_old-i_data[pos])*(i_old-i_data[pos]);
		}
	}
	n_avg /= (height*width);
	n2_avg = sqrt(n2_avg/(height*width));
	G_PRINTF( "Noise:\tmax=%g,min=%g,n_avg=%.2g,n2_avg=%.2g\r\n",n_max,n_min,n_avg,n2_avg );		
GIP_exit:

	return 0x0;
}

/*
	check Peak signal to noise ratio
	ISSUE-NEEDED:
	1、PSNR-grad measure how well the derivaies match

	Copyright 2008-present, Grusoft.
	v0.1	cys
		5/22/2008	
*/
double GIP_IMG_PSNR( int height,int width,double *u,int ldu,GIP_IMAGE *hIMG,double *temp )	{
	int i,j,pos,step=hIMG->widthStep/sizeof(uchar);
	double rmse,rmse_x,rmse_y,psnr=DBL_MAX,psnr_g=DBL_MAX,off;
	uchar* i_data = (uchar *)hIMG->imageData;

	rmse=0.0,	rmse_x=0.0,	rmse_y=0.0;
	ASSERT( height==hIMG->height && width==hIMG->width );
	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			pos = i*step+j;
			off = u[i*ldu+j]-i_data[pos];
			rmse += off*off;
		}
	}
	rmse = sqrt( rmse/height/width );
	psnr = 20*log(255/rmse)/log(10.0);

	return psnr;
}

/*
	return the empirical standard deviation

	Copyright 2008-present, Grusoft.
	v0.1	cys
		5/28/2008	
*/
double GIP_IMG_delta( GIP_IMAGE *hIMG )	{
	int i,j,pos,height=hIMG->height,width=hIMG->width,step=hIMG->widthStep/sizeof(uchar);
	uchar* i_data = (uchar *)hIMG->imageData;
	double delta=0.0,u_0=0.0,a;

	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			pos = i*step+j;
			a = i_data[pos];
			u_0 += a;
		}
	}
	u_0 /=(height*width);

	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			pos = i*step+j;
			a = i_data[pos];
			delta += (a-u_0)*(a-u_0);
		}
	}
	delta /=(height*width);
	delta = sqrt( delta );

	return delta;
}

/*
	v0.1	cys
		7/15/2005
*/
int fReadLine( FILE* fp,TCHAR sLine[MAX_LINE] )	{
	TCHAR c;
	int getLine=0,pos=0;
	while( !feof(fp ) )	{
		c=_fgettc(fp);
		if( ferror(fp) )	{
			break;
		}
		getLine = 1;
		if( c==_T('\r') || c==_T('\n') )
			break;
		sLine[pos++]=c;
	}
	sLine[pos]=_T('\0');
	ASSERT( pos < MAX_LINE );

	return getLine;
}

/*
	http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/segbench/seg-format.txt
	Copyright 2008-present, Grusoft

	v0.1	cys
		6/6/2008	
*/
void GIP_import_BSDSseg( IplImage *hIMG,char *sSegFile,int flag )	{
	int height=hIMG->height,width=hIMG->width,channels=hIMG->nChannels,step=hIMG->widthStep/sizeof(uchar);
	int j,r,c_0,c_1,no,isData=0,len,nSeg=0,gray;
	uchar* data = (uchar *)hIMG->imageData,u;
	FILE *fp=_tfopen( sSegFile,"r" );
	char sLine[MAX_LINE];

	if( fp==0x0 )
		GIP_ERROR( "GIP_import_BSDSseg",-1, "failed to open segmentation file." );

	while( fReadLine( fp,sLine ) != 0 )	{
		len=(int)(_tcslen( sLine ));
		if( len==0 )
			continue;
		if( _tcsstr(sLine,"width")!=0x0 )	{
			if( _stscanf( sLine+6,"%d",&no ) !=1 || no != width )
				GIP_ERROR( "GIP_import_BSDSseg",-1, "width failed." );
		}
		else if( _tcsstr(sLine,"height")!=0x0 )	{
			if( _stscanf( sLine+7,"%d",&no ) !=1 || no != height )
				GIP_ERROR( "GIP_import_BSDSseg",-1, "height failed." );
		}
		else if( _tcsstr(sLine,"segments")!=0x0 )	{
			if( _stscanf( sLine+9,"%d",&nSeg ) !=1 )
				GIP_ERROR( "GIP_import_BSDSseg",-1, "height failed." );
		}else if( _tcsstr(sLine,"gray")!=0x0  )	{
			if( _stscanf( sLine+5,"%d",&gray ) !=1 || gray != channels )
				GIP_ERROR( "GIP_import_BSDSseg",-1, "height failed." );
		}else 	if( isData==1 )	{
			if( _stscanf( sLine,"%d %d %d %d",&no,&r,&c_0,&c_1 ) !=4 )
				GIP_ERROR( "GIP_import_BSDSseg",-1, "height failed." );
				u = (UCHAR_MAX*no*1.0/nSeg);
				for( j = c_0; j <= c_1; j++ )	{
					data[r*step+j]=u;							
				}
//			for( i = 0; i < height; i++ )	{
//			}
		}
		else if( _tcsstr(sLine,"data")!=0x0 )	{
			isData=1;
		}
	}

	fclose( fp );
GIP_exit:
	;
}

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		11/4/2008	
*/
void GIP_IMG_segm( int height,int width,double *u,int ldu,GIP_FLOATING *u_0,int flag  )	{
	int r,c,block_size=1000,level=2,alg=GIP_ALG_PYRAMIDS,step,pos,pos_1,cur;
	double threshold1,threshold2;
	CvSeq *comp,*contours=0;
	CvMemStorage *storage;
	IplImage *hIMG_0,*hIMG_1,*markers=0;
	char sTitle[80];
    CvMat* color_tab;
	uchar* u_1,*u_m ;

	if( (hIMG_0=GIP_set_IplImage_d(height,width,u,ldu,0x0) )==GIP_NULL )
		GIP_ERROR( "GIP_IMG_segm",-1, "the image must be a single-channel byte image" );
/*
	switch( alg )	{
	case GIP_ALG_PYRAMIDS:
	//	hIMG_0->width &= -(1<<level);
	//	hIMG_0->height &= -(1<<level);
		hIMG_1 = cvCloneImage( hIMG_0 );
		threshold1 =255;
		threshold2 =30;
		storage = cvCreateMemStorage ( block_size );
	//	cvPyrSegmentation( hIMG_0, hIMG_1, storage, &comp,level, threshold1+1, threshold2+1 );
		break;
	case GIP_ALG_WATERSHED:
		markers = cvCreateImage( cvGetSize(hIMG_0), IPL_DEPTH_32S, 1 );
		cvZero( markers );
 		step=hIMG_0->widthStep/sizeof(uchar);
		u_1 = (uchar *)hIMG_0->imageData;
		u_m = (uchar *)markers->imageData;
		for( r = 0; r < height; r++ )	{
			for( c = 0; c < width; c++ )	{
				pos_1=r*step+c;
				cur = u_1[pos_1];	
				if( (c>0&&cur>=u_1[pos_1-1]) || (c<width-1&&cur>=u_1[pos_1+1])
					|| (r>0&&cur>=u_1[pos_1-step]) || (r<height-1&&cur>=u_1[pos_1+step]) )	{
					continue;	
				}
				//u_m[pos_1] = 255;
			}
		} 
		cvWatershed( hIMG_0, markers );
		break;
	}
 
	_stprintf( sTitle,".\\trace\\SEG_%d.jpg",height*width );
	if( !cvSaveImage(sTitle,hIMG_1) ) 		
		GIP_ERROR( "GIP_IMG_segm",-1, "Could not save" );

	if( BIT_TEST( flag,GIP_OPENCV_CREATEWND ) )	{
		_stprintf( sTitle,".\\trace\\SEG_%d.jpg",height*width );
 		cvNamedWindow( sTitle, 1 );
		cvShowImage( sTitle, hIMG_1 );
	}

	step=hIMG_1->widthStep/sizeof(uchar);
	u_1 = (uchar *)hIMG_1->imageData;
	for( r = 0; r < height; r++ )	{
		for( c = 0; c < width; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
		//	u[pos] = u_0[pos];
			pos_1=r*step+c;
			cur = u_1[pos_1];	
			if( (c>0&&cur!=u_1[pos_1-1]) || (c<width-1&&cur!=u_1[pos_1+1])
				|| (r>0&&cur!=u_1[pos_1-step]) || (r<height-1&&cur!=u_1[pos_1+step]) )	{
				u[pos]=0.0;	
			}
		}
	}
	cvReleaseImage( &hIMG_1 );
*/	cvReleaseImage( &hIMG_0 );
GIP_exit:
		;
}

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		11/4/2008	
*/
GIP_IMAGE *GIP_IMG_resize( GIP_IMAGE *hIMG_0,int height,int width,int flag )	{
	GIP_IMAGE *hIMG=GIP_NULL;

	if( (hIMG=cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,1))==0x0 )
		GIP_ERROR( "GIP_IMG_resize",-1, "the image must be a single-channel byte image" );
	ASSERT( 0 );
	cvResize( hIMG_0, hIMG,CV_INTER_CUBIC );
//	cvSaveImage( ".\\trace\\cv_resize.jpg",hIMG);
GIP_exit:
		;
	return hIMG;
}

/*
	调用cvResize

	v0.1	cys
		4/9/2009	
*/
void GIP_cv_resize( GIP_OBJECT *hObj_0,GIP_OBJECT *hObj,int type,int flag )	{
	int i,r,c,M=hObj_0->M,N=hObj_0->N,I,pos;
	int height=hObj->M,width=hObj->N;
	GIP_IMAGE *hIMG_0=GIP_NULL,*hIMG;
	uchar* img,*img_0;
	int step;
	GIP_FLOATING *data_0=hObj_0->data,*data=hObj->data;

	if( (hIMG_0=cvCreateImage(cvSize(N,M),IPL_DEPTH_8U,1))==0x0 )
	{	ASSERT( 0 );		goto GIP_exit;			}
	step=hIMG_0->widthStep/sizeof(uchar);
	img_0 = (uchar *)hIMG_0->imageData;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,N );
		img_0[r*step+c] =  data_0[pos];				
	}
	}			
	if( (hIMG=cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,1))==0x0 )
	{	ASSERT( 0 );		goto GIP_exit;			}
	cvResize( hIMG_0, hIMG,CV_INTER_CUBIC );
	step=hIMG->widthStep/sizeof(uchar);
	img = (uchar *)hIMG->imageData;
	for( r = 0; r < height; r++ )	{		
	for( c = 0; c < width; c++ )	{
		pos = GIP_RC2POS( r,c,width );
		data[pos] =  img[r*step+c];				
	}
	}				
	cvReleaseImage( &hIMG_0 );
	cvReleaseImage( &hIMG );

GIP_exit:
	return;
}
/*
	v=f-u+alpha

	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/14/2008	
*/
void GIP_IMG_v( IplImage*f_image,IplImage* hImg,int alpha,int flag )	{
	int i,j,height=hImg->height,width=hImg->width,step=hImg->widthStep/sizeof(uchar),i_val;
	double a,d_min=DBL_MAX,d_max=0.0,*v=(double*)GIP_alloc( sizeof(double)*height*width ),s=1.0;
	uchar* i_u = (uchar *)hImg->imageData,*i_f=(uchar *)f_image->imageData;

	if( hImg->depth != IPL_DEPTH_8U || hImg->nChannels !=1 )
		GIP_ERROR( "GIP_IMG_v",-1, "the image must be a single-channel byte image" );

	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			a = i_f[i*step+j]-i_u[i*step+j]+alpha;
			v[i*width+j]=a;
			d_min=MIN(d_min,a);	d_max=MAX(d_max,a);		
		}
	}
	ASSERT( d_max!=d_min );
	if( d_max-d_min>UCHAR_MAX )	{
		s = UCHAR_MAX*1.0/(d_max-d_min);
		G_PRINTF( "IMG_v:\tscaled=%g",s );
	}
	a = (UCHAR_MAX-(d_max-d_min)*s)/2-d_min*s;
	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			i_val = cvRound( v[i*width+j]*s+a );
			ASSERT( i_val>=0 && i_val <= UCHAR_MAX );
			i_u[i*step+j] = i_val;	
		}
	}
GIP_exit:
		;

	GIP_free( v );
}

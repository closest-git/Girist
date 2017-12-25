/*
	以下为cpp的函数	
	1 利用了template的一些特性
	2 调用c++代码,例如CxImage
*/

#include <Windows.h>
#include "../pch/GIP_def.h"
#include "../pch/GIP_package.h"
#include "../util/GIP_util.h"
#include "../util/gip_imag_transform.h"
#include "../util/gip_imag_property.h"
#include "2PassScale.h"
#include "./CxImage/ximage.h"
#include "../os/gip_win_util.h"
#ifdef _USE_OPENCV_
	#include "H:\tools\OpenCV_1.1\otherlibs\highgui\highgui.h"	
#endif
/*
	nearst neibor

	v0.1	cys
		8/10/2009	
*/
extern "C" int GIP_OBJ_resize_nn( GIP_OBJECT *hObj_0,GIP_OBJECT *hObj,int *I_buffer,int flag )	{
	int i,r,c,M=hObj_0->M,N=hObj_0->N,height=hObj->M,width=hObj->N,pos,nz=0,r_1,c_1,*count=I_buffer;
	GIP_FLOATING *data_0=hObj_0->data,*data=hObj->data;
	double s_r,s_c;
	
	GIP_MEMCLEAR( data,sizeof(GIP_FLOATING)*height*width );
	GIP_MEMCLEAR( count,sizeof(int)*height*width );

	s_r = height*1.0/M;		s_c =width*1.0/N;		
	for( r=0; r<M; r++ )	{
	for( c=0; c<N; c++ )	{
		r_1=(int)(r*s_r);			//ASSERT( 0<=r_1 && r_1<height );
		c_1=(int)(c*s_c);
		pos = GIP_RC2POS( r_1,c_1,width );
		data[pos] += data_0[nz++];
		count[pos] ++;		
	}	
	}

	for( i = height*width-1; i >= 0; i-- )	{
		data[i]/=count[i];
		//ASSERT( data[i]>=0 && data[i]<=255 );
	}

	return GIP_OK;
}


/*
	调用http://www.codeproject.com/KB/graphics/2_pass_scaling.aspx?df=100&forumid=159&exp=0

	注意：
		1 hObj的内存已经分配
		2 暂假定hObj_0对应灰度图像
	ISSUE-NEEDED:
		需扩充2PassScale算法至double数据类型

	Copyright 2008-present, Grusoft

	v0.1	cys
		4/3/2009	
*/
extern "C" int GIP_OBJ_resize( GIP_OBJECT *hObj_0,GIP_OBJECT *hObj,int type,int flag )	{
	int i,r,c,M=hObj_0->M,N=hObj_0->N,I,pos;
	int height=hObj->M,width=hObj->N;
	BYTE r_,g_,b_;

//	GIP_save_IplImage_d( M,N,GIP_path_(_T("resize"),0),hObj_0->data,N,0x0 );

	COLORREF *pOld,*pNew;
	GIP_FLOATING *data_0=hObj_0->data,*data=hObj->data;

	pOld=(COLORREF*)GIP_alloc( sizeof(COLORREF)*M*N );
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,N );
		I = G_DOUBLE2INT( data_0[pos] );
		ASSERT( I>=0 && I<256 );
		pOld[pos] = RGB(I,I,I );				
	}
	}	

	for( i = type; i < type+1; i++ )	{
		pNew = pOld;
		switch( i )	{
		case 1:	{
			C2PassScale<CBoxFilter> ScaleEngine_1;
			pNew = ScaleEngine_1.AllocAndScale( pNew,N,M,width,height );
				}	break;
		case 2:	{
			C2PassScale<CBilinearFilter> ScaleEngine_2;
			pNew = ScaleEngine_2.AllocAndScale( pNew,N,M,width,height );
				}	break;
		case 3:	{
			C2PassScale<CGaussianFilter> ScaleEngine_3;
			pNew = ScaleEngine_3.AllocAndScale( pNew,N,M,width,height );
				}	break;
		case 4:	{
			C2PassScale<CHammingFilter> ScaleEngine_4;
			pNew = ScaleEngine_4.AllocAndScale( pNew,N,M,width,height );
				}	break;
		case 5:	{
			C2PassScale<CBlackmanFilter> ScaleEngine_5;
			pNew = ScaleEngine_5.AllocAndScale( pNew,N,M,width,height );
				}	break;
		default:
			ASSERT( 0 );
			break;
		}
		ASSERT (NULL != pNew);
		for( r = 0; r < height; r++ )	{		
		for( c = 0; c < width; c++ )	{
			pos = GIP_RC2POS( r,c,width );
			I = G_DOUBLE2INT( data_0[pos] );
			r_=GetRValue(pNew[pos]);			g_=GetGValue(pNew[pos]);			b_=GetBValue(pNew[pos]);
			data[pos] = r_;				
		}
		}	
		delete pNew;			pNew=NULL;
//		GIP_save_IplImage_d( hObj->M,hObj->N,GIP_path_(_T("resize"),i),hObj->data,hObj->N,0x0 );
	}
	GIP_free( pOld );

	return 0;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		5/18/2009	
*/
extern "C" int G_ReadImage( TCHAR *sPath,GIP_OBJECT *hObj,int flag )	{
	int height,width,step,type,r,c,pos,ret=GIP_OK;
	BYTE b;
	TCHAR ext[_MAX_EXT];
	double mean,devia;

	if( BIT_TEST( flag,GIP_OPENCV_LIB) )	{
		uchar *img=GIP_NULL;
		GIP_IMAGE * hIMG=GIP_NULL;
		hIMG=GIP_LoadImage( sPath,0x0,0x0 );
		if( hIMG ==0x0 )	{
			ret=-3;		goto END;		
//			GIP_ERROR( "GIP_LoadImage",-1, "the image must be a single-channel byte image" );
		}
		height=hIMG->height,width=hIMG->width;
		step=hIMG->widthStep/sizeof(uchar);
		GIP_OBJ_init( hObj,height,width,0x0 );
		pos = 0;
		img = (uchar *)hIMG->imageData;
		for( r = 0; r < height; r++ )	{
			for( c = 0; c < width; c++ )	{
				hObj->data[pos++] = img[r*step+c];
			}
		}	
		cvReleaseImage( &hIMG );
	}else	{
		_DIR_split( sPath,ext,3 );
		_tcslwr( ext );
		if ( _tcscmp(ext,_T(""))==1 )	{
			type = CxImage::GetTypeIdFromName(ext);
			CxImage image(sPath, type);
			if( !image.IsValid() )	{
				ret=-2;		goto END;		
			}

			height = image.GetHeight( ),width=image.GetWidth( );
			pos = 0;
			GIP_OBJ_init( hObj,height,width,0x0 );
			for( r=0; r<height; r++ )	{
			for( c=0; c<width; c++ )	{
				b = image.GetPixelGray( c,r );
				hObj->data[pos++] = b;
			}	
			}
		}else	{	
			ret=-1;		goto END;		
		}
	}

END:
	return ret;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		5/31/2009	
*/
extern "C" int G_SaveImage( TCHAR *sPath,GIP_OBJECT *hObj,int flag )	{
	int  M=hObj->M,N=hObj->N,r,c,pos,ret=GIP_OK;
	BYTE b;

	if( BIT_TEST( flag,GIP_OPENCV_LIB) )	{
		GIP_save_IplImage_d( hObj->M,hObj->N,sPath,hObj->data,hObj->N,0x0 );
	}else{
		CxImage image(N,M,32);
		if( !image.IsValid() )	{
			ret=-2;		goto END;		
		}

		pos = 0;
		for( r=0; r<M; r++ )	{
		for( c=0; c<N; c++ )	{
			b = (BYTE)(hObj->data[pos++]);
			image.SetPixelColor ( c,r,RGB(b,b,b) );
		}	
		}
		bool retval = image.Save( sPath, CXIMAGE_FORMAT_BMP );
	}
	

END:
	return ret;
}




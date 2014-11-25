/*
	以下为cpp的函数	
	1 利用了template的一些特性
	2 调用c++代码
*/

#include <Windows.h>
#include "../pch/GIP_def.h"
#include "../util/GIP_util.h"
#include "../pch/GIP_package.h"
#include "2PassScale.h"

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

	GIP_save_IplImage_d( M,N,GIP_path_(_T("resize"),0),hObj_0->data,N,0x0 );

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
		GIP_save_IplImage_d( hObj->M,hObj->N,GIP_path_(_T("resize"),i),hObj->data,hObj->N,0x0 );
	}
	GIP_free( pOld );

	return 0;
}




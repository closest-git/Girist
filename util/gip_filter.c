#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include "../PCH/GIP_def.h"
#include "GIP_filter.h"
#include <cxcore.h>
#include <cv.h>

/*
	[32] P.98
//		sobel:	delta=|(z7+2*z8+z9)-(z1+2*z2+z3)| + |(z3+2*z6+z9)-(z1+2*z4+z7)|
//		robert:	delta=|z9-z5| + |z8-z6|

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/25/2008	
*/
void linear_spatial_filters( int M,int N,double *val,double *g,int type ){
	int r,c,pos;
	double vx,vy,a,a_max,a_min,s;

//	memset( g,0x0,sizeof(double)*M*N );
	for( r = 0; r < M; r++ )	{		//初始化：1.0表示没有梯度变化
	for( c = 0; c < N; c++ )	{		
		pos=r*N+c;
		g[pos]=1.0;
	}
	}

	a_max = 0.0;		a_min=DBL_MAX;
/*
	for( c = 1; c < N-1; c++ )	{		//列优先存储
	for( r = 1; r < M-1; r++ )	{
		pos=c*M+r;
		vy = val[pos]-val[pos-1];
		vx = val[pos]-val[pos-M];
		a = 1.0/(1.0+vx*vx+vy*vy);
		a_max=MAX( a,a_max );		a_min=MIN( a,a_min );
		g[pos] = a;
	}
	}
*/
	for( r = 1; r < M-1; r++ )	{
	for( c = 1; c < N-1; c++ )	{		//行优先存储
		pos=r*N+c;
		vx = val[pos]-val[pos-1];
		vy = val[pos]-val[pos-N];
		a = 1.0/(1.0+vx*vx+vy*vy);
		a_max=MAX( a,a_max );		a_min=MIN( a,a_min );
		g[pos] = a;
	}
	}
	
	ASSERT( a_max != a_min );
	s = a_max-a_min;
	for( r = 1; r < M-1; r++ )	{
	for( c = 1; c < N-1; c++ )	{
		pos=r*N+c;
		g[pos] =(g[pos]-a_min) / s;		//归一化
	}
	}
}

/*
	Rearrange the quadrants of Fourier image so that the origin is at
	the image center
	src & dst arrays of equal size & type

	Copy from INTEL OPENCV
	v0.1	cys
		9/29/2008	
*/
void cvShiftDFT(CvArr * src_arr, CvArr * dst_arr )
{
    CvMat * tmp;
    CvMat q1stub, q2stub;
    CvMat q3stub, q4stub;
    CvMat d1stub, d2stub;
    CvMat d3stub, d4stub;
    CvMat * q1, * q2, * q3, * q4;
    CvMat * d1, * d2, * d3, * d4;

    CvSize size = cvGetSize(src_arr);
    CvSize dst_size = cvGetSize(dst_arr);
    int cx, cy;

    if(dst_size.width != size.width || 
       dst_size.height != size.height){
        cvError( CV_StsUnmatchedSizes, "cvShiftDFT", "Source and Destination arrays must have equal sizes", __FILE__, __LINE__ );   
    }

    if(src_arr==dst_arr){
        tmp = cvCreateMat(size.height/2, size.width/2, cvGetElemType(src_arr));
    }
    
    cx = size.width/2;
    cy = size.height/2; // image center

    q1 = cvGetSubRect( src_arr, &q1stub, cvRect(0,0,cx, cy) );
    q2 = cvGetSubRect( src_arr, &q2stub, cvRect(cx,0,cx,cy) );
    q3 = cvGetSubRect( src_arr, &q3stub, cvRect(cx,cy,cx,cy) );
    q4 = cvGetSubRect( src_arr, &q4stub, cvRect(0,cy,cx,cy) );
    d1 = cvGetSubRect( src_arr, &d1stub, cvRect(0,0,cx,cy) );
    d2 = cvGetSubRect( src_arr, &d2stub, cvRect(cx,0,cx,cy) );
    d3 = cvGetSubRect( src_arr, &d3stub, cvRect(cx,cy,cx,cy) );
    d4 = cvGetSubRect( src_arr, &d4stub, cvRect(0,cy,cx,cy) );

    if(src_arr!=dst_arr){
        if( !CV_ARE_TYPES_EQ( q1, d1 )){
            cvError( CV_StsUnmatchedFormats, "cvShiftDFT", "Source and Destination arrays must have the same format", __FILE__, __LINE__ ); 
        }
        cvCopy(q3, d1, 0);
        cvCopy(q4, d2, 0);
        cvCopy(q1, d3, 0);
        cvCopy(q2, d4, 0);
    }
    else{
        cvCopy(q3, tmp, 0);
        cvCopy(q1, q3, 0);
        cvCopy(tmp, q1, 0);
        cvCopy(q4, tmp, 0);
        cvCopy(q2, q4, 0);
        cvCopy(tmp, q2, 0);
    }
}

/*

	Copyright 2008-present, Grusoft.
	v0.1	cys
		9/28/2008	
*/
int GIP_IMG_DFT_( IplImage* hIMG,int flag )	{
	int height=hIMG->height,width=hIMG->width,i,j,step,dft_M,dft_N;
    CvMat *mS,*mD,*mI;
	double a_0,a_1,a,s,dT=0.2;
	IplImage* img_i = cvCloneImage( hIMG );
	uchar *i_data= (uchar *)img_i->imageData;

    dft_M = cvGetOptimalDFTSize( hIMG->height - 1 );
    dft_N = cvGetOptimalDFTSize( hIMG->width - 1 );
	mS = cvCreateMat(dft_M,dft_N,CV_32FC1);		    cvZero(mS);
	mD = cvCreateMat(dft_M,dft_N,CV_32FC1);			cvZero(mD);
	mI = cvCreateMat(dft_M,dft_N,CV_32FC1);			cvZero(mI);
	step=hIMG->widthStep/sizeof(uchar);

	a_0=DBL_MAX,		a_1=DBL_MIN;
	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			a=i_data[i*step+j];
			cvSetReal2D( mS,i,j,i_data[i*step+j] );
			a_0=MIN(a,a_0);
			a_1=MAX(a,a_1);
		}
	}
    cvDFT(mS,mD,CV_DXT_FORWARD,height);     //DFT 离散傅立叶变换
	for( i = 0; i < dft_M; i++ )	{
	for( j = 0; j < dft_N; j++ )	{
		a_0 = (i*PI/dft_M);		a_1=(j*PI/dft_N);
		s = exp(-(a_0*a_0+a_1*a_1)*dT);
		a=cvGetReal2D( mD,i,j );
		cvSetReal2D( mD,i,j,a*s );
	}	
	}
    cvDFT(mD,mI,CV_DXT_INV_SCALE,height);		//DFT 离散傅立叶逆变换

	a_0=DBL_MAX,		a_1=DBL_MIN;
	for( i = 0; i < height; i++ )	{
		for( j = 0; j < width; j++ )	{
			a=cvGetReal2D( mI,i,j );
			i_data[i*step+j]=a;
			a_0=MIN(a,a_0);
			a_1=MAX(a,a_1);
		}
	}

//	cvNamedWindow("inverse", 1);
//    cvShowImage("inverse", img_i);

    cvWaitKey(-1);

	return 0;
}

/*

	Copyright 2008-present, Grusoft.
	v0.1	cys
		9/28/2008	

int GIP_IMG_DFT_( IplImage* hIMG,int flag )	{
    IplImage * realInput;
    IplImage * imaginaryInput;
    IplImage * complexInput;
    int dft_M, dft_N,r,c,step;
    CvMat* dft_A, tmp;
    IplImage * image_Re;
    IplImage * image_Im;
    double m, M,dT=1.0,a0,a1,s;
	CvPoint2D64f *m_;

    realInput = cvCreateImage( cvGetSize(hIMG), IPL_DEPTH_64F, 1);
    imaginaryInput = cvCreateImage( cvGetSize(hIMG), IPL_DEPTH_64F, 1);
    complexInput = cvCreateImage( cvGetSize(hIMG), IPL_DEPTH_64F, 2);

    cvScale(hIMG, realInput, 1.0, 0.0);
    cvZero(imaginaryInput);
    cvMerge(realInput, imaginaryInput, NULL, NULL, complexInput);

    dft_M = cvGetOptimalDFTSize( hIMG->height - 1 );
    dft_N = cvGetOptimalDFTSize( hIMG->width - 1 );

    dft_A = cvCreateMat( dft_M, dft_N, CV_64FC2 );
    image_Re = cvCreateImage( cvSize(dft_N, dft_M), IPL_DEPTH_64F, 1);
    image_Im = cvCreateImage( cvSize(dft_N, dft_M), IPL_DEPTH_64F, 1);

    // copy A to dft_A and pad dft_A with zeros
    cvGetSubRect( dft_A, &tmp, cvRect(0,0, hIMG->width, hIMG->height));
    cvCopy( complexInput, &tmp, NULL );
    if( dft_A->cols > hIMG->width )
    {
        cvGetSubRect( dft_A, &tmp, cvRect(hIMG->width,0, dft_A->cols - hIMG->width, hIMG->height));
        cvZero( &tmp );
    }

    // no need to pad bottom part of dft_A with zeros because of
    // use nonzero_rows parameter in cvDFT() call below

    cvDFT( dft_A, dft_A, CV_DXT_FORWARD, complexInput->height );
    cvDFT( dft_A, dft_A, CV_DXT_INVERSE, complexInput->height );
	m_ = (CvPoint2D64f*)(dft_A->data.ptr);
	step=dft_A->step/sizeof(CvPoint2D64f);
	for( r = 0; r < dft_M; r++ )	{
	for( c = 0; c < dft_N; c++ )	{
		a0 = (r*PI/dft_M);		a1=(c*PI/dft_N);
		s = exp(-a0*a0-a1*a1)*dT;
		(m_[r*step+c]).x *= s;
		(m_[r*step+c]).y *= s;
	}	
	}
	cvNamedWindow("after...", 0);
	cvSplit( dft_A, image_Re, image_Im, 0, 0 );
    cvMinMaxLoc(image_Re, &m, &M, NULL, NULL, NULL);
    cvScale(image_Re, image_Re, 1.0/(M-m), 1.0*(-m)/(M-m));
    cvShowImage("after..", image_Re);

  //   cvNamedWindow("win", 0);
 //   cvShowImage("win", hIMG);
   cvNamedWindow("magnitude", 0);
 
    // Split Fourier in real and imaginary parts
    cvSplit( dft_A, image_Re, image_Im, 0, 0 );

    // Compute the magnitude of the spectrum Mag = sqrt(Re^2 + Im^2)
    cvPow( image_Re, image_Re, 2.0);
    cvPow( image_Im, image_Im, 2.0);
    cvAdd( image_Re, image_Im, image_Re, NULL);
    cvPow( image_Re, image_Re, 0.5 );

    // Compute log(1 + Mag)
   cvAddS( image_Re, cvScalarAll(1.0), image_Re, NULL ); // 1 + Mag
    cvLog( image_Re, image_Re ); // log(1 + Mag)


    // Rearrange the quadrants of Fourier image so that the origin is at
    // the image center
    cvShiftDFT( image_Re, image_Re );

    cvMinMaxLoc(image_Re, &m, &M, NULL, NULL, NULL);
    cvScale(image_Re, image_Re, 1.0/(M-m), 1.0*(-m)/(M-m));
    cvShowImage("magnitude", image_Re);

    cvWaitKey(-1);

    return 0;
}
*/
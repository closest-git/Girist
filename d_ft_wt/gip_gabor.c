/*
	Wiesel [2] deduced that simple cells are sensitive to specific orientations with approximate bandwidths of
		30'.
*/


#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include <math.h>
#include <limits.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/gip_imag_transform.h"
#include "../util/gip_util.h"
#include "gip_gabor.h"
#include "mkl_dfti.h"

/*
	Gabor Spatial filters 
		[REF]: Personal identification based on iris texture analysis

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/26/2008	
*/
void GIP_GABOR_init( GIP_GABOR *hGabor,int M,int N,double dr,double dc,GIP_GABOR_TYPE type,int flag )	{
	int r,c,pos;
	double f,g,m,gs,rs,cs;
	GIP_FLOATING *filter;
	GIP_MATRIX *hMat;

	memset( hGabor,0x0,sizeof(GIP_GABOR) );
//	hGabor->type = GIP_GABOR_MALI;

	hGabor->M=M;				hGabor->N=N;
	hGabor->sigma_r=dr;			hGabor->sigma_c=dc;
	hGabor->f = f = 2.0;
	hGabor->filter=GIP_alloc( sizeof(GIP_FLOATING)*M*N );		filter=hGabor->filter;

	gs = 1.0/2/PI/dr/dc;		
	for( r = 0; r < M; r++ )	{
	for( c = 0; c < N; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		rs=r/dr;		cs=c/dc;
		g=gs*exp( -0.5*(rs*rs+cs*cs) );
		m=cos(2*PI*f*sqrt(r*r+c*c) );
		filter[pos] = g*m;
	}
	}

	hGabor->hMat=GIP_alloc( sizeof(GIP_MATRIX) );			hMat=hGabor->hMat;
	hMat->M=M;						hMat->N=N;	
	hMat->data=filter;

//	GIP_save_IplImage_d( hMat->M,hMat->N,".\\trace\\gabor.jpg",hMat->data,hMat->N,0x0 );

	return ;
}	

/*
	Log Gabor filters in the frequency domain 
	基于LOG GABOR filter的特殊性，单独列出

	gabor filter参数的具体解释可参见 http://matlabserver.cs.rug.nl/edgedetectionweb/web/edgedetection_params.html

	注意:
		暂设定为sig_C=0.5*f0

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/24/2009	
*/
void GIP_LOG_GABOR_init( GIP_GABOR *hGabor,int M,int N,double lenda,double sig_R,double sig_C,GIP_GABOR_TYPE type,int flag )	{
	int N_half,i;
	double a,f0,s,radius;
	GIP_FLOATING *filter;
	GIP_MATRIX *hMat;

	memset( hGabor,0x0,sizeof(GIP_GABOR) );
	hGabor->type = type;				hGabor->flag = flag;

	hGabor->M=M;					hGabor->N=N;			
	hGabor->sigma_r=sig_R;			hGabor->sigma_c=sig_C;
	hGabor->lenda = lenda;
	if( BIT_TEST( flag,GIP_GABOR_1D ) )		ASSERT( M==1 && sig_R==0 );
/*	if( BIT_TEST( flag,GIP_GABOR_FREQ ) )	{
		hGabor->filter=GIP_alloc( sizeof(GIP_FLOATING)*M*N*2 );	
		hGabor->re=hGabor->filter;			hGabor->im=hGabor->filter+M*N;	
	}else	{
		hGabor->filter=GIP_alloc( sizeof(GIP_FLOATING)*M*N );	
		hGabor->re=hGabor->filter;			hGabor->im=GIP_NULL;	
	}*/
	hGabor->filter=GIP_alloc( sizeof(GIP_FLOATING)*M*N );
	GIP_MEMCLEAR(hGabor->filter,sizeof(GIP_FLOATING)*M*N );
	filter=hGabor->filter;
	N_half = (int)(N/2);
    f0 = 1.0/lenda;							//Centre frequency of filter.
	hGabor->sigma_c=0.5*f0;
	s = 2.0*log(0.5)*log(0.5);
	for( i = 0 ; i < N_half+1; i++ )	{	//logGabor(1:ndata/2+1) = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2)); 
		radius = i==0 ? 1 : i*0.5/N_half;				//Frequency values 0 - 0.5
		a = log(radius/f0 );
		filter[i] = exp( -a*a/s );
	}    
    filter[0] = 0;  
	
	hGabor->hMat=GIP_alloc( sizeof(GIP_MATRIX) );			hMat=hGabor->hMat;
	hMat->M=M;						hMat->N=N;	
	hMat->data=filter;

//	GIP_save_IplImage_d( hMat->M,hMat->N,".\\trace\\gabor.jpg",hMat->data,hMat->N,0x0 );

	return ;
}	

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/26/2008	
*/
void GIP_GABOR_clear( GIP_GABOR *hGabor )	{
	if( hGabor->filter != GIP_NULL )	
	{	GIP_free( hGabor->filter );	hGabor->filter=GIP_NULL; }
	if( hGabor->hMat != GIP_NULL )	{
		GIP_free( hGabor->hMat );	hGabor->hMat=GIP_NULL;
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/28/2008	
*/
void GIP_GBANK_init(  )	{

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/28/2008	
*/
void GIP_GBANK_clear(  )	{
}

/*
	ASSERT( BIT_TEST( flag,GIP_GABOR_FREQ ) )

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/24/2009	
*/
void GIP_GABOR_filter( GIP_GABOR *hGabor,int M,int N,GIP_FLOATING *data,double *re,double *im,int flag )	{
	int r,c,pos,ldu;
	GIP_FLOAT_Z *u_A;
	long Status,lengths[2],strides_in[3];
	double s,*filter=hGabor->filter;
    DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;

//	if( BIT_TEST( hGabor->flag,GIP_GABOR_1D ) )		ASSERT( M==1 );
	ASSERT( BIT_TEST( hGabor->flag,GIP_GABOR_FREQ ) );

	ldu=N;
	u_A=GIP_alloc( sizeof(GIP_FLOAT_Z)*M*ldu );			ASSERT( u_A!=GIP_NULL );	

	lengths[0] = M;			lengths[1] = N;
    strides_in[0]= 0;		strides_in[1]= ldu;		strides_in[2]= 1;
	Status = DftiCreateDescriptor( &Desc_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 2, lengths);
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiSetValue(Desc_Handle, DFTI_INPUT_STRIDES, strides_in);	//Set input strides parameters 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiCommitDescriptor( Desc_Handle );					//Commit Dfti descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			s=1.0;		//s=(r+c)%2==0 ? 1:-1;
			u_A[pos].re=(data[pos]);			u_A[pos].im=0.0;
		}
	}
	Status = DftiComputeForward( Desc_Handle, u_A);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			s=filter[c];		
			u_A[pos].re *=s;			u_A[pos].im*=s;
		}
	}	
	Status = DftiComputeBackward( Desc_Handle, u_A);					//   Compute Backward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			re[pos] = u_A[pos].re;			im[pos] = u_A[pos].im;
		}
	}

END_OF_TEST:
	Status = DftiFreeDescriptor(&Desc_Handle);							//   Free DFTI descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))
	{	G_PRINTF(" TEST FAIL\n");							}	

	GIP_free( u_A );

}
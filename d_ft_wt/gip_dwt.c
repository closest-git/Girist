#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include <math.h>
#include <limits.h>
#include "../PCH/GIP_def.h"
#include "../util/gip_imag_transform.h"
#include "../util/gip_util.h"
#include "gip_dwt.h"
#include "mkl_dfti.h"

static double _ld_db4_[8]={
	-1.059740178499728e-002,3.288301166698295e-002,3.084138183598697e-002,-1.870348117188811e-001,
	-2.798376941698385e-002,6.308807679295904e-001,7.148465705525415e-001,2.303778133088552e-001
};
static double _ld_sym4_[8]=	{			
	-7.576571478927333e-002,-2.963552764599851e-002,4.976186676320155e-001,8.037387518059161e-001,
	2.978577956052774e-001,-9.921954357684722e-002,-1.260396726203783e-002,3.222310060404270e-002	
};
static double _ld_bior68_[18]=	{			
	0,1.908831736481291e-003,-1.914286129088767e-003,-1.699063986760234e-002,
	1.193456527972926e-002,4.973290349094079e-002,-7.726317316720414e-002,-9.405920349573646e-002,
	4.207962846098268e-001,8.259229974584023e-001,4.207962846098268e-001,-9.405920349573646e-002,
	-7.726317316720414e-002,4.973290349094079e-002,1.193456527972926e-002,-1.699063986760234e-002,
	-1.914286129088767e-003,1.908831736481291e-003
};
static double _hd_bior68_[18]=	{
	0,0,0,1.442628250562444e-002,
	-1.446750489679015e-002,-7.872200106262882e-002,4.036797903033992e-002,4.178491091502746e-001,
	-7.589077294536542e-001,4.178491091502746e-001,4.036797903033992e-002,-7.872200106262882e-002,
	-1.446750489679015e-002,1.442628250562444e-002,0,0,0,0
};
static double _ld_jpeg97_[10]=	{
	0,0.02674875741080976,-0.01686411844287495,-0.07822326652898785,
	0.2668641184428723,0.6029490182363579,0.2668641184428723,-0.07822326652898785,
	-0.01686411844287495,0.02674875741080976
};
static double _hd_jpeg97_[10]=	{
	0,-0.09127176311424948,0.05754352622849957,0.5912717631142470,
	-1.115087052456994,0.5912717631142470,0.05754352622849957,-0.09127176311424948,0,0
};
/*
	[REF] WAVEFAST.m

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/11/2008	
*/
void GIP_WaveFilter( GIP_MATRIX **ml,GIP_MATRIX **mh,int type,int flag )	{
	int i;

	switch( type )	{
	case GIP_DWT_DB4:			//'db4'
		*ml = GIP_OBJ_create( 1,8,_ld_db4_,0x0 );
		*mh = GIP_OBJ_create( 1,8,GIP_NULL,0x0 );
		for( i = 0; i < 8; i++ )	{
			(*mh)->data[7-i]=cos( PI*(i))*_ld_db4_[i];
		}
/*		GIP_OBJ_init( ml,1,8,0x0 );
		GIP_OBJ_init( mh,1,8,0x0 );
		for( i = 0; i < 8; i++ )	{
			ml->data[i]=_ld_db4_[i];
			mh->data[7-i]=cos( PI*(i))*_ld_db4_[i];
		}*/
//	   lr = ld;    lr(end:-1:1) = ld;
//	   hr = cos(pi * t) .* ld;
	   break;
	case GIP_DWT_HAAR:
	case GIP_DWT_DB1:		//{'haar', 'db1'}
		*ml = GIP_OBJ_create( 1,2,GIP_NULL,0x0 );
		*mh = GIP_OBJ_create( 1,2,GIP_NULL,0x0 );
		(*ml)->data[0]=1.0/sqrt(2);			(*ml)->data[1]=1.0/sqrt(2);
		(*mh)->data[0]=-1.0/sqrt(2);		(*mh)->data[1]=1.0/sqrt(2);
//	   ld = [1 1]/sqrt(2);     hd = [-1 1]/sqrt(2);
//	   lr = ld;                hr = -hd;
		break;
	   
	   
/*	case 'sym4'
	   t = (0:7);
	   hd = ld;    hd(end:-1:1) = cos(pi * t) .* ld;
	   lr = ld;    lr(end:-1:1) = ld;
	   hr = cos(pi * t) .* ld;
	   
	case 'bior6.8'
	   t = (0:17);
	   lr = cos(pi * (t + 1)) .* hd;
	   hr = cos(pi * t) .* ld;
	   
	case 'jpeg9.7'
	   ld = [];
	   hd = [];
	   t = (0:9);
	   lr = cos(pi * (t + 1)) .* hd;
	   hr = cos(pi * t) .* ld;*/
	default:
		ASSERT( 0 );
	   break;
	}
}

/* 
	Algebraically, convolution is the same operation as multiplying the polynomials whose coefficients are the elements of u and v.

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/11/2008	
*/
void GIP_Convol_( GIP_MATRIX *mC,GIP_MATRIX *mA,GIP_MATRIX *mB,int flag )	{
	int M,N,m,n,r,c,r_0,r_1,c_0,c_1,posA,posB,pos;
	GIP_FLOATING *data,*dataA,*dataB;

	mC->M= mA->M+mB->M-1;			M=mC->M;
	mC->N= mA->N+mB->N-1;			N=mC->N;
	mC->data=GIP_alloc( sizeof(GIP_FLOATING)*M*N );		
	data=mC->data;
	dataA=mA->data;					dataB=mB->data;
	pos = 0;
	for( m=0; m <M; m++ )	{
	for( n=0; n <N; n++ )	{
		data[pos] = 0.0;
		r_0=MAX(m+1-mB->M,0);		r_1 = MIN( mA->M,m+1 );
		c_0=MAX(n+1-mB->N,0);		c_1 = MIN( mA->N,n+1 );
		for( r = r_0; r < r_1; r++)	{
		for( c = c_0; c < c_1; c++)	{
			posA=GIP_RC2POS( r,c,mA->N );
			posB=GIP_RC2POS( m-r,n-c,mB->N );
			data[pos] += dataA[posA]*dataB[posB];
		}
		}
		pos++;
	}
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/5/2009	
*/
void GIP_spectrum_( int M,int N,GIP_FLOAT_Z *u_A,int type,int flag )	{
	int r,c,pos;
	double *spectrum,re,im;

	spectrum=GIP_alloc( sizeof(double)*M*N );
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,N );
			re = u_A[pos].re;		im = u_A[pos].im;
			spectrum[pos] = 0.5*log(1.0+sqrt( re*re+im*im ));			//对数变换
		}
	}

	GIP_save_IplImage_d( M,N,GIP_path_(_T("spectrum"),flag),spectrum,N,0x0 );

	GIP_free( spectrum );
}

/* 
	fourier matching for translated rotated & scaled compensating 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/6/2009	
*/
void GIP_fourier_matching( int M,int N,GIP_MATRIX *mA,GIP_MATRIX *mB,int flag )	{
	int r,c,s_r,s_c,pos,pos_1,ldu,r1,c1,r2=-1,c2=-1;
	GIP_FLOATING *dataA,*dataB;
	GIP_FLOAT_Z *u_A,*u_B,*u_C;
	long Status,lengths[2],strides_in[3];
    DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;
	double re,im,s,a,a_max=0.0,a_mean=0.0,b_mean=0.0;

//	GIP_save_IplImage_d( mA->M,mA->N,".\\trace\\dft_A.jpg",mA->data,mA->N,0x0 );
//	GIP_save_IplImage_d( mB->M,mB->N,".\\trace\\dft_B.jpg",mB->data,mB->N,0x0 );
	dataA=mA->data;					dataB=mB->data;

	ldu=N;
	u_A=GIP_alloc( sizeof(GIP_FLOAT_Z)*M*ldu*3 ); 	
	u_B=u_A+M*ldu;		u_C=u_B+M*ldu;
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
			a_mean += dataA[pos];			b_mean += dataB[pos];
		}
	}
	a_mean/=(M*N);							b_mean/=(M*N);			
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			s=1.0;		//s=(r+c)%2==0 ? 1:-1;
//			u_A[pos].re=(dataA[pos])/M/N*s;			u_A[pos].im=0.0;
//			u_B[pos].re=(dataB[pos])/M/N*s;			u_B[pos].im=0.0;
			u_A[pos].re=(dataA[pos]-a_mean)/M/N*s;			u_A[pos].im=0.0;
			u_B[pos].re=(dataB[pos]-b_mean)/M/N*s;			u_B[pos].im=0.0;
		}
	}
	Status = DftiComputeForward( Desc_Handle, u_A);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiComputeForward( Desc_Handle, u_B);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			u_C[pos].re=u_A[pos].re*u_B[pos].re+u_A[pos].im*u_B[pos].im;
			u_C[pos].im=u_A[pos].re*u_B[pos].im-u_A[pos].im*u_B[pos].re;
		}
	}
	Status = DftiComputeBackward( Desc_Handle, u_C);					//   Compute Backward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			re = u_C[pos].re;		im = u_C[pos].im;
			ASSERT( fabs(im)<FLT_EPSILON );
			a = sqrt( re*re+im*im );
			if( a > a_max )	{
				a_max = a;		r2=r;	c2=c;
			}
		}
	}


END_OF_TEST:
	Status = DftiFreeDescriptor(&Desc_Handle);							//   Free DFTI descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))
	{	G_PRINTF(" TEST FAIL\n");							}	

	GIP_free( u_A );

}

/* 
	mC=mA*mB
	filter based on fourier transform

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/11/2008	
*/
void GIP_filter_ft( int M,int N,GIP_MATRIX *mC,GIP_MATRIX *mA,GIP_MATRIX *mB,int flag )	{
	int r,c,pos,ldu;
	GIP_FLOATING *data,*dataA,*dataB;
	GIP_FLOAT_Z *u_A,*u_B,*u_C;
	long Status,lengths[2],strides_in[3];
    DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;
	double re,im,s;

//	GIP_save_IplImage_d( mA->M,mA->N,".\\trace\\dft_A.jpg",mA->data,mA->N,0x0 );
//	GIP_save_IplImage_d( mB->M,mB->N,".\\trace\\dft_B.jpg",mB->data,mB->N,0x0 );
	mC->M = M;				mC->N = N;
//	mC->data=GIP_alloc( sizeof(GIP_FLOATING)*M*N );		
	data=mC->data;
	dataA=mA->data;					dataB=mB->data;

	ldu=N;
	u_A=GIP_alloc( sizeof(GIP_FLOAT_Z)*M*ldu*3 ); 	
	u_B=u_A+M*ldu;		u_C=u_B+M*ldu;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			s=1.0;		//s=(r+c)%2==0 ? 1:-1;
			u_A[pos].re=dataA[pos]/M/N*s;			u_A[pos].im=0.0;
			u_B[pos].re=dataB[pos]/M/N*s;			u_B[pos].im=0.0;
		}
	}

	lengths[0] = M;			lengths[1] = N;
    strides_in[0]= 0;		strides_in[1]= ldu;		strides_in[2]= 1;
	Status = DftiCreateDescriptor( &Desc_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 2, lengths);
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiSetValue(Desc_Handle, DFTI_INPUT_STRIDES, strides_in);	//Set input strides parameters 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiCommitDescriptor( Desc_Handle );					//Commit Dfti descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiComputeForward( Desc_Handle, u_A);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiComputeForward( Desc_Handle, u_B);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			u_B[pos].im *= -1;
			u_C[pos].re=u_A[pos].re*u_B[pos].re-u_A[pos].im*u_B[pos].im;
			u_C[pos].im=u_A[pos].re*u_B[pos].im+u_A[pos].im*u_B[pos].re;
			s=(r+c)%2==0 ? 1:-1;
//			u_C[pos].re *=-1;			u_C[pos].im*=-1;
		}
	}
	GIP_spectrum_( M,N,u_A,0x0,0x0 );
	GIP_spectrum_( M,N,u_B,0x0,0x1 );
	GIP_spectrum_( M,N,u_C,0x0,0x2 );
	Status = DftiComputeBackward( Desc_Handle, u_C);					//   Compute Backward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
//			re = u_A[pos].re;		im = u_A[pos].im;
			re = u_B[pos].re;		im = u_B[pos].im;
//			mC->data[pos] = sqrt( re*re+im*im );
			mC->data[pos] = 0.5*log(1.0+sqrt( re*re+im*im ));			//对数变换
//			mC->data[pos] = u_C[pos].re;
//			im = u_C[pos].im;		ASSERT( fabs(im)<FLT_EPSILON );
		}
	}
//	GIP_save_IplImage_d( mC->M,mC->N,".\\trace\\dft_0.jpg",mC->data,mC->N,0x0 );


END_OF_TEST:
	Status = DftiFreeDescriptor(&Desc_Handle);							//   Free DFTI descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))
	{	G_PRINTF(" TEST FAIL\n");							}	

	GIP_free( u_A );

}

/* 
	mB=|fourier(mA)|

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/14/2008	
*/
void GIP_fourier_spectrum( int M,int N,GIP_MATRIX *mA,double *dataB,int flag )	{
	int r,c,pos,ldu;
	GIP_FLOATING *dataA;
	GIP_FLOAT_Z *u_A;
	long Status,lengths[2],strides_in[3];
    DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;
	double re,im,s;

//	GIP_save_IplImage_d( mA->M,mA->N,".\\trace\\dft_A.jpg",mA->data,mA->N,0x0 );
	dataA=mA->data;					

	ldu=N;
	u_A=GIP_alloc( sizeof(GIP_FLOAT_Z)*M*ldu ); 	
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			//s=1.0;		
			s=(r+c)%2==0 ? 1:-1;
			u_A[pos].re=dataA[pos]/M/N*s;			u_A[pos].im=0.0;
		}
	}

	lengths[0] = M;			lengths[1] = N;
    strides_in[0]= 0;		strides_in[1]= ldu;		strides_in[2]= 1;
	Status = DftiCreateDescriptor( &Desc_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 2, lengths);
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiSetValue(Desc_Handle, DFTI_INPUT_STRIDES, strides_in);	//Set input strides parameters 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiCommitDescriptor( Desc_Handle );					//Commit Dfti descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiComputeForward( Desc_Handle, u_A);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			re = u_A[pos].re;		im = u_A[pos].im;
			dataB[pos] = sqrt( re*re+im*im );
	//		dataB[pos] = 0.5*log(1.0+sqrt( re*re+im*im ));
		}
	}

//	GIP_save_IplImage_d( mC->M,mC->N,".\\trace\\dft_0.jpg",mC->data,mC->N,0x0 );


END_OF_TEST:
	Status = DftiFreeDescriptor(&Desc_Handle);							//   Free DFTI descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))
	{	G_PRINTF(" TEST FAIL\n");							}	

	GIP_free( u_A );

}

	


/* 
%-------------------------------------------------------------------%
function y = symconv(x, h, type, fl, keep)
	Convolve the rows or columns of x with h, downsample,
	and extract the center section since symmetrically extended.
	[REF] WAVEFAST.m

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/10/2008	
*/
void _sym_conv( GIP_MATRIX *mY,GIP_MATRIX *mX,GIP_MATRIX *mH,int type,int k_M,int k_N )	{
	char sPath[]="E:\\GIPack\\dump\\conv.data";
	if( type== GIP_CONVOLV_ROW )	{
		GIP_Convol_( mY,mX,mH,0x0 );		//	   y = conv2(x, h);
//		D_ouput( mY->M*mY->N,mY->data,sPath );
		GIP_MATRIX_downsample( mY,mH->N,mH->N+2*k_N-2,2,GIP_MATRIX_COL );
		;//y = y(:, 1:2:end);
		;//y = y(:, fl / 2 + 1:fl / 2 + keep(2));
	}else	{
		int hm=mH->M,hn=mH->N;
		mH->M=hn,			mH->N=hm;
		GIP_Convol_( mY,mX,mH,0x0 );		//	   y = conv2(x, h');
		GIP_MATRIX_downsample( mY,mH->M,mH->M+2*k_M-2,2,GIP_MATRIX_ROW );
		;//y = y(1:2:end, :);
		;//y = y(fl / 2 + 1:fl / 2 + keep(1), :);
		mH->M=hm,			mH->N=hn;
	}
}

/* 
	[REF] WAVEBACK.m

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/14/2008	

void _sym_convup( GIP_MATRIX *mY,GIP_MATRIX *mX,GIP_MATRIX *f1,GIP_MATRIX *f2,int fln,int k_M,int k_N )	{
	char sPath[]="E:\\GIPack\\dump\\conv.data";
	if( type== GIP_CONVOLV_ROW )	{
		GIP_Convol_( mY,mX,mH,0x0 );		//	   y = conv2(x, h);
//		D_ouput( mY->M*mY->N,mY->data,sPath );
		GIP_MATRIX_downsample( mY,mH->N,mH->N+2*k_N-2,2,GIP_MATRIX_COL );
		;//y = y(:, 1:2:end);
		;//y = y(:, fl / 2 + 1:fl / 2 + keep(2));
	}else	{
		int hm=mH->M,hn=mH->N;
		mH->M=hn,			mH->N=hm;
		GIP_Convol_( mY,mX,mH,0x0 );		//	   y = conv2(x, h');
		GIP_MATRIX_downsample( mY,mH->M,mH->M+2*k_M-2,2,GIP_MATRIX_ROW );
		;//y = y(1:2:end, :);
		;//y = y(fl / 2 + 1:fl / 2 + keep(1), :);
		mH->M=hm,			mH->N=hn;
	}
}
*/

/*
	顺序:	
		H,V,D
	Copyright 2008-present, Grusoft.
	v0.1	cys	
		12/13/2008	
*/
void GIP_DWT_display( GIP_DWT *hDWT,int level,int flag )	{
	int M,N,lM,lN,i;
	GIP_MATRIX *mH,*mV,*mD,*arrMat=hDWT->arrMat,*hMat;

	ASSERT( level<=hDWT->level );
	M = hDWT->M;			N = hDWT->N;
	for( i = 0; i < level; i++ )	{
		mH=arrMat+3*i;			
		mV=arrMat+3*i+1;			
		mD=arrMat+3*i+2;
		lM=mH->M;				lN=mH->N;
		GIP_MATRIX_fill( hMat,mH,0,lN,0x0 );
		GIP_MATRIX_fill( hMat,mV,lM,0,0x0 );
		GIP_MATRIX_fill( hMat,mD,lM,lN,0x0 );
	}

//	GIP_save_IplImage_d( hMat->M,hMat->N,sTitle,hMat->data,hMat->N,0x0 );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys	
		12/13/2008	
*/
void GIP_DWT_clear( GIP_DWT *hDWT )	{
	int i;

	if( hDWT->arrMat != GIP_NULL )	{
		for( i = 0; i < hDWT->nMat; i++ )
			GIP_OBJ_clear( hDWT->arrMat+i );
		GIP_free(hDWT->arrMat);				
		hDWT->arrMat=GIP_NULL;	
	}
	if( hDWT->hp != GIP_NULL )	{
		GIP_OBJ_clear( hDWT->hp );
		GIP_free(hDWT->hp);				hDWT->hp=GIP_NULL;	
	}
	if( hDWT->lp != GIP_NULL )	{
		GIP_OBJ_clear( hDWT->lp );
		GIP_free(hDWT->lp);				hDWT->lp=GIP_NULL;	
	}
	GIP_OBJ_clear( &(hDWT->tempMat) );

	hDWT->hObj=GIP_NULL;
}

/*
	注意:
		struc与mat一一对应

	Copyright 2008-present, Grusoft.
	v0.1	cys	
		12/13/2008	
*/
void GIP_DWT_init( GIP_DWT *hDWT,GIP_OBJECT *hObj,int level,GIP_DWT_FILTER filter,int flag )	{

	memset( hDWT,0x0,sizeof(GIP_DWT) );

	hDWT->hObj=hObj;
	hDWT->filter=filter;
	hDWT->level=level;

	hDWT->max_mat=3*level+1;
	hDWT->arrMat=GIP_alloc( sizeof(GIP_MATRIX)*hDWT->max_mat );
	hDWT->nMat=0;			
}

/*
[**]
	宝宝19个月,活泼又可爱:)
[**]
	
	perform multi-level 2-dimensional fast wavelet transform
	[REF] WAVEFAST.m

	注意：
		1 hObj的数据不变

	Copyright 2008-present, Grusoft.
	v0.1	cys	
		12/10/2008	
*/
void GIP_dwt( GIP_OBJECT *hObj,double *temp,int level,GIP_DWT_FILTER filter,int flag )	{
	int i,k_M,k_N,rank,M=hObj->M,N=hObj->N,nMat;
	GIP_DWT dwt;
	GIP_MATRIX *arrMat,*hp,*lp,*app,*app_0,rows;//coefs,hp,lp;
	char sPath[80]="E:\\GIPack\\dump\\dwt.data";
	
	rank = G_DOUBLE2INT( log(MAX(M,N))/log(2.0) );
	if( level<0 || level>rank )
		level = rank;
	if( level==0 )
		return;
	GIP_DWT_init( &dwt,hObj,level,filter,0x0 );
	arrMat=dwt.arrMat;
	lp=dwt.lp;				hp=dwt.hp;
	GIP_WaveFilter( &lp,&hp,filter,0x0 );
	rank=lp->N;

//	app=dwt.app = GIP_OBJ_create( M,N,hObj->data,0x0 );
	app_0=hObj;			app=&(dwt.tempMat);
	nMat = 0;
	for( i = 0; i < level; i++ )	{
		GIP_MATRIX_symextend( app,app_0,rank-1,&k_M,&k_N,0x0 );					//[app, keep] = symextend(app, fl);
//		D_ouput( app->M*app->N,app->data,sPath );

		_sym_conv( &rows,app,hp,GIP_CONVOLV_ROW,k_M,k_N );	//rows = symconv(app, hp, 'row', fl, keep);
//		D_ouput( rows.M*rows.N,rows.data,sPath );
		_sym_conv( arrMat+nMat,&rows,hp,GIP_CONVOLV_COL,k_M,k_N );	//coefs = symconv(rows, hp, 'col', fl, keep);	
		nMat++;			   //c = [coefs(:)' c];    s = [size(coefs); s];
		_sym_conv( arrMat+nMat,&rows,lp,GIP_CONVOLV_COL,k_M,k_N );	//coefs = symconv(rows, lp, 'col', fl, keep);
		nMat++;				//c = [coefs(:)' c];

		_sym_conv( &rows,app,lp,GIP_CONVOLV_ROW,k_M,k_N );	//rows = symconv(app, lp, 'row', fl, keep);
		_sym_conv( arrMat+nMat,&rows,hp,GIP_CONVOLV_COL,k_M,k_N );	//coefs = symconv(rows, hp, 'col', fl, keep);
		nMat++;				//c = [coefs(:)' c];    s = [size(coefs); s];
		_sym_conv( arrMat+nMat,&rows,lp,GIP_CONVOLV_COL,k_M,k_N );	//app = symconv(rows, lp, 'col', fl, keep);
		//c = [coefs(:)' c];    s = [size(coefs); s];
		app_0 = arrMat+nMat;
		_stprintf( sPath,"E:\\GIPack\\dump\\dwt_%d.data",i+1 );
		D_ouput( app_0->M*app_0->N,app_0->data,sPath );
	}
	nMat++;	
	ASSERT( nMat<=dwt.max_mat );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/11/2008	
*/
void GIP_Convol_testing( )	{
	GIP_FLOATING a[9]={
		0.8147,0.9134,0.2785,0.9058,0.6324,0.5469,0.1270,0.0975,0.9575
	};
	GIP_FLOATING b[16]={
		0.9649,0.4854,0.9157,0.0357,0.1576,0.8003,0.7922,0.8491,
		0.9706,0.1419,0.9595,0.9340,0.9572,0.4218,0.6557,0.6787
	};
	GIP_MATRIX mA,mB,mC;
	int m,n,pos;

	GIP_OBJ_init( &mA,3,3,0x0 );
	memcpy( mA.data,a,sizeof(GIP_FLOATING)*9 );
	GIP_OBJ_init( &mB,4,4,0x0 );
	memcpy( mB.data,b,sizeof(GIP_FLOATING)*16 );

	GIP_Convol_( &mC,&mA,&mB,0x0 );
	pos=0;
	for( m = 0; m < mC.M; m++ )	{
		for( n = 0; n < mC.N; n++ )	{
			G_PRINTF( "%.4lf ",mC.data[pos] );
			pos++;
		}
		G_PRINTF( "\r\n" );
	}

	GIP_OBJ_clear( &mA );
	GIP_OBJ_clear( &mB );
	GIP_OBJ_clear( &mC );
/*
	0.7861    1.2768    1.4581    1.0007    0.2876    0.0099
    1.0024    1.8458    3.0844    2.5151    1.5196    0.2560
    1.0561    1.9824    3.5790    3.9432    2.9708    0.7587
    1.6790    2.0772    3.0052    3.7511    2.7593    1.5129
    0.9902    1.1000    2.4492    1.6082    1.7976    1.2655
    0.1215    0.1469    1.0409    0.5540    0.6941    0.6499
*/	;
}

/*		for( i=0;i<hObj->M*hObj->N;i++ )		hObj->data[i]=1.0;
		double _magic_mat[64]={64, 2,3,61,60,6,7,57, 9,55,54,12,13,51,50,16,17,47,46,20,21,43,42,24,40,26,27,37,36,30,31,33
			,32,34,35,29,28,38,39,25,41,23,22,44,45,19,18,48,49,15,14,52,53,11,10,56, 8,58,59, 5, 4,62,63, 1};
		hObj=GIP_OBJ_create( 8,8,_magic_mat,0x0 );
*/

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/22/2008	
*/
void GIP_MASK_testing( GIP_OBJECT *hIMG,int flag )	{
	GIP_OBJECT *hMask,*hEdge,*hRes;
	int i,pos,pos_0,r,c,r1,c1,r2,c2,M=hIMG->M,N=hIMG->N,MN=M*N,ldu=N,nz;
	int *m_off,I_0=0,I_1=200;
	GIP_FLOATING *m_data,*data=hIMG->data,*e_data;
	double w,w_max,gr,gc,g,g_thrsh=10;	

	hMask = GIP_OBJ_create( M,N,GIP_NULL,0x0 );
	m_data=hMask->data;
	if( 0 )	{
		GIP_Circle( hMask,M/2,N/2,74,GIP_IMAG_EDGE );
	}else	{
		for( i = 0; i < MN; i++ )	m_data[i]=0.0;
		m_data[GIP_RC2POS(M/2,N/2,N)] = 1.0;
	//	m_data[GIP_RC2POS(0,0,N)] = 1.0;
	}
	nz = 0;
	m_off= GIP_alloc( sizeof(int)*MN );
	pos=GIP_RC2POS( M/2,N/2,ldu );
	for( i = 0; i < MN; i++ )	{
		if( m_data[i]==0.0 )		{
			continue;
		}else	{
			m_data[i] = 1.0;
			m_off[nz++] = i-pos;
		}
	}

	hEdge = GIP_OBJ_create( M,N,GIP_NULL,0x0 );
	GIP_Circle( hEdge,M/2,N/2,74,GIP_IMAG_EDGE );
	e_data = hEdge->data;
/*	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,N );
		e_data[pos] = 0.0;
		if( data[pos]>I_1 || data[pos]<I_0 )		continue;
		gr = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_X,0x0 );
		gc = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_Y,0x0 );
		g = sqrt( gr*gr+gc*gc );
		if( g < g_thrsh )	
			continue;
		e_data[pos] = 1.0;
	}
	}*/

	hRes = GIP_OBJ_create( M,N,GIP_NULL,0x0 );
	if( 0 )	{
		w_max = 0.0;
		for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos_0 = GIP_RC2POS( r,c,N );
			w = 0;
			for( i = 0; i < nz; i++ )	{
				pos = pos_0+m_off[i];
				r1=GIP_POS2R( pos,ldu);			c1=GIP_POS2C( pos,ldu);
				if( !GIP_RC_VALID(r1,0,M-1) || !GIP_RC_VALID(c1,0,N-1) )
					continue;
				w += e_data[pos];
			}
			hRes->data[pos_0]=w;
			if( w>w_max )	{
				w_max = w;		r2=r;		c2=c;
			}
		}
		}
	}else	{
		GIP_filter_ft( M,N,hRes,hEdge,hMask,0x0 );
//		GIP_Convol_( GIP_MATRIX *mC,GIP_MATRIX *mA,GIP_MATRIX *mB,int flag )
		w_max = 0.0;
		for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos_0 = GIP_RC2POS( r,c,N );
			w = hRes->data[pos_0];
			if( w>w_max )	{
				w_max = w;		r2=r;		c2=c;
			}
		}
		}
	}
//	GIP_save_IplImage_d( hRes->M,hRes->N,".\\trace\\dft_res.jpg",hRes->data,hRes->N,0x0 );

	GIP_OBJ_clear( hMask );				GIP_free( hMask );			
	GIP_OBJ_clear( hEdge );				GIP_free( hEdge );	
	GIP_free( m_off );	
}
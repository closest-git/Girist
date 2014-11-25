/*
	这里指广义的decomposition，即包含图像，纹理，噪声等的分解，也含图像分割

	1 IMAGE SEGMENTATION USING ADAPTIVE FINITE ELEMENTS
		根据u动态调整网格
	2 
	v0.1	10/28/2008
		cys
*/
#include <math.h>
#include <FLOAT.h>
#include <limits.h>
#include <stdio.h>
#include "../pch/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/GIP_util.h"
#include "../util/GIP_filter.h"
#include "../util/GIP_heap.h"
#include "../util/GIP_fdm.h"
#include "../util/GIP_ss.h"
#include "GIP_decomposition.h"
#include "../fem/GIP_fem.h"
#include "mkl_dfti.h"
#include "../util/Compact_DataStruc.h"

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/30/2008			
*/
int _SEG_CV_force( int M,int N,int ldu,double*u,double *f,int flag )	{
	int n1=0,n2=0,dim=M*N,i;
	double c1=0.0,c2=0.0;

	for( i = 0; i < dim; i++ )	{
		if( u[i]>=0.0 )		
		{	c1+=u[i];		n1++;	}
		else
		{	c2+=u[i];		n2++;	}
	}
	c1/=n1;			c2/=n2;

	for( i = 0; i < dim; i++ )	{
//		u[i] = u[i]>=0.0 ? (f[i]-c1)*(f[i]-c1) :-(f[i]-c2)*(f[i]-c2);
		u[i] = (f[i]-c1)*(f[i]-c1)-(f[i]-c2)*(f[i]-c2);
	}

	return 0;
}

/*
	Original Discrete Fourier Transform

	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/6/2008			
*/
int _DFT_test_0( int M,int N,int ldu,double*f_0,int flag )	{
	int r,c,err=0,pos,s=1,u,v,isDoubleF=0;
    double  *d_temp,ksi_1,ksi_2,*f,*Fr,*Fi,real,imag,off;
    DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;
	char sTitle[80];
	FILE *fp;

	d_temp=GIP_alloc( sizeof(double)*M*ldu*5 ); 		
	f=d_temp;
	Fr=f+M*ldu;					Fi=Fr+M*ldu;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			f[pos]=s*f_0[pos];
//			s *= -1;
		}
	}

	for( u = 0; u < M; u++ )	{
	for( v = 0; v < N; v++ )	{
		real=0.0;			imag=0.0;
		for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			ksi_1 = -r*u*2.0*PI/M-c*v*2.0*PI/N;		
			real+=f[pos]*cos(ksi_1);		imag+=f[pos]*sin(ksi_1);
		}
		}	
		pos=GIP_RC2POS(u,v,ldu );
		Fr[pos]=real/M/N;		Fi[pos]=imag/M/N;
		if( isDoubleF==1 )	Fi[pos]*=-1;	
	}
	}
	for( r = 0; r < M; r++ )	{
	for( c = 0; c < N; c++ )	{
		real=0.0;			imag=0.0;
		for( u = 0; u < M; u++ )	{
		for( v = 0; v < N; v++ )	{
			pos=GIP_RC2POS(u,v,ldu );
			ksi_1 = (r*u*2.0*PI/M+c*v*2.0*PI/N);		
			if( isDoubleF==1 )	ksi_1 = -(r*u*2.0*PI/M+c*v*2.0*PI/N);	
			real+=Fr[pos]*cos(ksi_1)-Fi[pos]*sin(ksi_1);		
			imag+=Fr[pos]*sin(ksi_1)+Fi[pos]*cos(ksi_1);
		}
		}	
		pos=GIP_RC2POS(r,c,ldu );
		f[pos]=real;		
		ASSERT( fabs(imag)<FLT_EPSILON );
	}
	}

	_stprintf( sTitle,"..\\data\\dft_%d.txt",M*N );	
	fp = _tfopen( sTitle,"w" );
	off=0.0;
	s=1;		
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			fprintf( fp,"%d\t%d\t%.15lf\t%.15lf\n",r+1,c+1,Fr[pos],Fi[pos] );
			f[pos]*=s;
//			s *= -1;
			off += (f[pos]-f_0[pos])*(f[pos]-f_0[pos]);
		}
	}
	fclose( fp );
	_stprintf( sTitle,".\\trace\\dft_%d.jpg",M*N );	
	GIP_save_IplImage_d( M,N,sTitle,f,N,0x0 );
	off=sqrt(off);

	GIP_free( d_temp );
	return err;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/28/2008			
*/
int _DFT_test( int M,int N,int ldu,double*u_0,int flag )	{
	int r,c,err=0,pos,s=1,isDoubleF=0,t;
	long    Status,lengths[2],strides_in[3];
    double  *d_temp,ksi_1,ksi_2,*u,off,Scale,*Fr,*Fi;
    DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;
	char sTitle[80];
	FILE *fp;

	if( 1 )	{
		d_temp=GIP_alloc( sizeof(double)*M*ldu*3 ); 		
	}
	_stprintf( sTitle,"..\\data\\dft_%d.txt",M*N );	
	if( isDoubleF==1 )	fp = _tfopen( sTitle,"r" );
	u=d_temp;			Fr=u+M*ldu,			Fi=Fr+M*ldu;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			u[pos]=s*u_0[pos]/M/N;
			if( isDoubleF==1 )	fscanf( fp,"%d\t%d\t%lf\t%lf\n",&t,&t,Fr+pos,Fi+pos );
			s *= -1;
		}
	}
	if( isDoubleF==1 )	fclose(fp);

	lengths[0] = N;			lengths[1] = M;
    strides_in[0]= 0;		strides_in[1]= 1;		strides_in[2]= ldu;
	Status = DftiCreateDescriptor( &Desc_Handle, DFTI_DOUBLE,DFTI_REAL, 2, lengths);	//Create Dfti descriptor for 2D double precision  transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiSetValue(Desc_Handle, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT);
 	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiSetValue(Desc_Handle, DFTI_INPUT_STRIDES, strides_in);	//Set input strides parameters 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto FREE_DESCRIPTOR;
//	Scale = 1.0/(double)(M*N);											//   Set Scale number for Backward transform
//	Status = DftiSetValue(Desc_Handle, DFTI_FORWARD_SCALE, Scale);
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiCommitDescriptor( Desc_Handle );					//Commit Dfti descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiComputeForward( Desc_Handle, u);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	if( isDoubleF==1 )	{
		for( r = 0; r < M; r++ )	{		//for Pack Format
			if( r==0 || (r+1==M && (r+1)%2==0) )	{
				pos=GIP_RC2POS(r,0,ldu );			//(r,0)
//				ASSERT( fabs(u[pos]-Fr[GIP_RC2POS((r+1)/2,0,ldu )])<FLT_EPSILON );					
				s = 1;			
				for( c = 1; c < N; c++ )	{
					pos=GIP_RC2POS(r,c,ldu );			
//					if( s==1 )
//						ASSERT( fabs(u[pos]-Fr[GIP_RC2POS((r+1)/2,(c+1)/2,ldu )])<FLT_EPSILON );		
//					else
//						ASSERT( fabs(u[pos]-Fi[GIP_RC2POS((r+1)/2,(c+1)/2,ldu )])<FLT_EPSILON );		
					u[pos]=s*u[pos];
					s *= -1;
				}
				s = 1;
			}else	{
				for( c = 0; c < N; c++ )	{
					pos=GIP_RC2POS(r,c,ldu );
//					if( s==1 )
//						ASSERT( fabs(u[pos]-Fr[GIP_RC2POS((r+1)/2,c,ldu )])<FLT_EPSILON );		
//					else
//						ASSERT( fabs(u[pos]-Fi[GIP_RC2POS((r+1)/2,c,ldu )])<FLT_EPSILON );		
					u[pos]=s*u[pos];
				}
				s *= -1;
			}
		}
	}
	Status = DftiCommitDescriptor( Desc_Handle );						//   Commit Dfti descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	if( isDoubleF==1 )
		Status = DftiComputeForward( Desc_Handle, u);					//   Compute Backward transform 
	else
		Status = DftiComputeBackward( Desc_Handle, u);					//   Compute Backward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	_stprintf( sTitle,".\\trace\\dft_%d.jpg",M*N );	
	off=0.0;
	s=1;		
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			u[pos]*=s;
			s *= -1;
			off += (u[pos]-u_0[pos])*(u[pos]-u_0[pos]);
		}
	}
	GIP_save_IplImage_d( M,N,sTitle,u,N,0x0 );
	off=sqrt(off);

//	maxerr = check_result_2d_d(x_in, x_exp, m, n, strides_in);			//   Check result
FREE_DESCRIPTOR:
	if( !DftiErrorClass(Status, DFTI_NO_ERROR) )
		G_PRINTF( "%s",DftiErrorMessage(Status) );
	Status = DftiFreeDescriptor(&Desc_Handle);							//   Free DFTI descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))
	{	G_PRINTF(" TEST FAIL\n");							}
END_OF_TEST:
	if( 1 )		GIP_free( d_temp );
	if( err==0 )	{
	}
	return err;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/6/2008			
*/
int _DFT_test_c( int M,int N,int ldu,double*u_0,int flag )	{
	int r,c,err=0,pos,s=1,isDoubleF=1,t;
	long    Status,lengths[2],strides_in[3];
    GIP_FLOAT_Z  *u_c;
	double ksi_1,ksi_2,*u,off,Scale,a;
    DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;
	char sTitle[80];
	FILE *fp;

	u_c=GIP_alloc( sizeof(GIP_FLOAT_Z)*M*ldu*3 ); 		
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			u_c[pos].re=s*u_0[pos]/M/N;			u_c[pos].im=0.0;
//			s *= -1;
		}
	}

	lengths[0] = M;			lengths[1] = N;
    strides_in[0]= 0;		strides_in[1]= ldu;		strides_in[2]= 1;
	Status = DftiCreateDescriptor( &Desc_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 2, lengths);
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiSetValue(Desc_Handle, DFTI_INPUT_STRIDES, strides_in);	//Set input strides parameters 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto FREE_DESCRIPTOR;
	Status = DftiCommitDescriptor( Desc_Handle );					//Commit Dfti descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiComputeForward( Desc_Handle, u_c);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	if( isDoubleF==1 )	{
		for( r = 0; r < M; r++ )	{		//for Pack Format
			for( c = 0; c < N; c++ )	{
				pos=GIP_RC2POS(r,c,ldu );			
				u_c[pos].im *= -1;
			}
		}
	}
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	if( isDoubleF==1 )
		Status = DftiComputeForward( Desc_Handle, u_c);					//   Compute Backward transform 
	else
		Status = DftiComputeBackward( Desc_Handle, u_c);					//   Compute Backward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	_stprintf( sTitle,".\\trace\\dft_%d.jpg",M*N );	
	off=0.0;
	s=1;		
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			a = u_c[pos].re*s;
//			s *= -1;
			off += (a-u_0[pos])*(a-u_0[pos]);
			a = u_c[pos].im*s;
			off += a*a;
		}
	}
//	GIP_save_IplImage_d( M,N,sTitle,u,N,0x0 );
	off=sqrt(off);

//	maxerr = check_result_2d_d(x_in, x_exp, m, n, strides_in);			//   Check result
FREE_DESCRIPTOR:
	if( !DftiErrorClass(Status, DFTI_NO_ERROR) )
		G_PRINTF( "%s",DftiErrorMessage(Status) );
	Status = DftiFreeDescriptor(&Desc_Handle);							//   Free DFTI descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))
	{	G_PRINTF(" TEST FAIL\n");							}
END_OF_TEST:
	if( 1 )		GIP_free( u_c );
	if( err==0 )	{
	}
	return err;
}

/*
	[REF] Image Restoration and Decompostion via Bounded Total Variation and Negative Hilbert-Sobolev Spaces
	u=DFT( (f^-u^)/(1+|KSI|) )

	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/27/2008			
*/
int _DECOM_Hs_force( int M,int N,int ldu,double*u,double *f,int flag )	{
	int r,c,ret=0,pos,M_even=(M/2)*2;
	long    Status,lengths[2],strides_in[3];
    double  Scale,ksi_1,ksi_2,gauss=1.0,s;
    DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;
	GIP_FLOAT_Z  *u_c,*f_c;

	u_c=GIP_alloc( sizeof(GIP_FLOAT_Z)*M*ldu*2 ); 	f_c=u_c+M*ldu;	
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			u_c[pos].re=u[pos]/M/N;			u_c[pos].im=0.0;
			f_c[pos].re=f[pos]/M/N;			f_c[pos].im=0.0;
		}
	}

	lengths[0] = M;			lengths[1] = N;
    strides_in[0]= 0;		strides_in[1]= N;		strides_in[2]= 1;
	Status = DftiCreateDescriptor( &Desc_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 2, lengths);
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiSetValue(Desc_Handle, DFTI_INPUT_STRIDES, strides_in);	//Set input strides parameters 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto FREE_DESCRIPTOR;
	Status = DftiCommitDescriptor( Desc_Handle );					//Commit Dfti descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiComputeForward( Desc_Handle, u_c);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiComputeForward( Desc_Handle, f_c);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );			
//			ksi_1=r*1.0/M;		ksi_2=c*1.0/N;
			ksi_1=r;			ksi_2=c;
			s = (1.0+ksi_1*ksi_1+ksi_2*ksi_2);	
			s = s;
//			gauss=exp(-(ksi_1*ksi_1+ksi_2*ksi_2)/2/0.8 );
			u_c[pos].re = (f_c[pos].re-u_c[pos].re*gauss)/s*gauss;
			u_c[pos].im = -(f_c[pos].im-u_c[pos].im*gauss)/s*gauss;
		}
	}
/*	for( r = 0; r < M; r++ )	{		//for Pack Format
		ksi_1=floor(r/2)/M;
		if( r==0 || (r+1==M && (r+1)%2==0) )	{
			pos=GIP_RC2POS(r,0,ldu );			ksi_2=0.0;					//(r,0)
//			gauss=exp(-(ksi_1*ksi_1+ksi_2*ksi_2)/2/0.8 );
			u[pos]=(f[pos]-u[pos]*gauss)/(1.0+ksi_1*ksi_1+ksi_2*ksi_2)*gauss;	
			s = 1;			
			for( c = 1; c < N; c++ )	{
				pos=GIP_RC2POS(0,c,ldu );			
				ksi_2=floor(c/2)/N;
//				gauss=exp(-(ksi_1*ksi_1+ksi_2*ksi_2)/2/0.8 );
				u[pos]=s*(f[pos]-u[pos]*gauss)/(1.0+ksi_1*ksi_1+ksi_2*ksi_2)*gauss;
				s *= -1;
			}
		}else	{
			for( c = 0; c < N; c++ )	{
				pos=GIP_RC2POS(r,c,ldu );
				ksi_2=c*1.0/N;
//				gauss=exp(-(ksi_1*ksi_1+ksi_2*ksi_2)/2/0.8 );
				u[pos]=s*(f[pos]-u[pos]*gauss)/(1.0+ksi_1*ksi_1+ksi_2*ksi_2)*gauss;
			}
			s *= -1;
		}
	}*/

	DftiComputeForward( Desc_Handle, u_c);			
//	DftiComputeForward( Desc_Handle, f_c);			
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;

FREE_DESCRIPTOR:
	Status = DftiFreeDescriptor(&Desc_Handle);							//   Free DFTI descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))
	{	G_PRINTF(" TEST FAIL\n");							}
END_OF_TEST:

	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
//			ASSERT( fabs(u_c[pos].im) < FLT_EPSILON );
//			ASSERT( fabs(u_c[pos].re-(f[pos]-u[pos])) < FLT_EPSILON );
			u[pos] = u_c[pos].re;
		}
	}

	GIP_free( u_c );	

	return ret;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/12/2008			
*/
int _DECOM_filer_force( int M,int N,int ldu,double*u,double *f,int flag )	{
	int r,c,ret=0,pos;
	double w_0=1.0,w_1,w_2,s_0,s_1,s_2,*u_c,*f_u;

	u_c=GIP_alloc( sizeof(double)*(M*ldu*2) ); 	f_u=u_c+M*ldu;
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			f_u[pos]=f[pos]-u[pos];
			u_c[pos]=f_u[pos];
		}
	}

	w_1=exp(-1.0);		
	w_2=exp(-2.0);
	s_0=w_0/(w_0+w_1*4+w_2*4);
	s_1=w_1/(w_0+w_1*4+w_2*4);
	s_2=w_2/(w_0+w_1*4+w_2*4);
	for( r = 1; r < M-1; r++ )	{
		for( c = 1; c < N-1; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			u_c[pos]=f_u[pos]*s_0+(f_u[pos-1]+f_u[pos+1]+f_u[pos-ldu]+f_u[pos+ldu])*s_1
				+(f_u[pos-ldu-1]+f_u[pos-ldu+1]+f_u[pos+ldu-1]+f_u[pos+ldu+1])*s_2;
		}
	}

	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			u[pos] = u_c[pos];
		}
	}

	GIP_free( u_c );	

	return ret;
}

/*
	[REF] Image Restoration and Decompostion via Bounded Total Variation and Negative Hilbert-Sobolev Spaces
	<f-u>s=(f^-u^)(f^-u^)/(1+|KSI|)

	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/11/2008			
*/
int _DECOM_Hs_energy( int M,int N,int ldu,double*u,double *f,double *eng,int flag )	{
	int r,c,ret=0,pos,M_even=(M/2)*2;
	long    Status,lengths[2],strides_in[3];
    double  ksi_1,ksi_2,gauss=1.0,re,im,sum,s;
    DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;
	GIP_FLOAT_Z  *u_c,*f_c;

	u_c=GIP_alloc( sizeof(GIP_FLOAT_Z)*M*ldu*2 ); 	f_c=u_c+M*ldu;	
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );
			u_c[pos].re=u[pos]/M/N;			u_c[pos].im=0.0;
			f_c[pos].re=f[pos]/M/N;			f_c[pos].im=0.0;
		}
	}

	lengths[0] = M;			lengths[1] = N;
    strides_in[0]= 0;		strides_in[1]= N;		strides_in[2]= 1;
	Status = DftiCreateDescriptor( &Desc_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 2, lengths);
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiSetValue(Desc_Handle, DFTI_INPUT_STRIDES, strides_in);	//Set input strides parameters 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto FREE_DESCRIPTOR;
	Status = DftiCommitDescriptor( Desc_Handle );					//Commit Dfti descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiComputeForward( Desc_Handle, u_c);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;
	Status = DftiComputeForward( Desc_Handle, f_c);					//Compute Forward transform 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))			goto END_OF_TEST;

	sum=0.0;
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos=GIP_RC2POS(r,c,ldu );			
//			ksi_1=r*1.0/M;		ksi_2=c*1.0/N;
			ksi_1=r;			ksi_2=c;
			s = (1.0+ksi_1*ksi_1+ksi_2*ksi_2);	
			s = s;
//			gauss=exp(-(ksi_1*ksi_1+ksi_2*ksi_2)/2/0.8 );
			re = (f_c[pos].re-u_c[pos].re*gauss)/s*gauss;
			im = -(f_c[pos].im-u_c[pos].im*gauss)/s*gauss;
			sum += re*re+im*im;
		}
	}

FREE_DESCRIPTOR:
	Status = DftiFreeDescriptor(&Desc_Handle);							//   Free DFTI descriptor 
	if(! DftiErrorClass(Status, DFTI_NO_ERROR))
	{	G_PRINTF(" TEST FAIL\n");							}
END_OF_TEST:

	GIP_free( u_c );	
	*eng= sqrt(sum);
	return ret;
}

/*
	[REF] A nonlinear primal-dual method for total variation-based image restoration (1999)

	注意：

	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/28/2008	
*/
void GIP_DECOMP_init_( GIP_IMAGE *hIMG,GIP_DECOMP *hDECOMP )	{
	int M,N,MN,ret=GIP_OK,t_len;
	GIP_FLOATING *val=GIP_NULL;
	double beta,*temp;
	GIP_OBJECT *hObj=GIP_NULL;

	memset( hDECOMP,0x0,sizeof(GIP_DECOMP) );
	hDECOMP->hObj=GIP_alloc( sizeof(GIP_OBJECT) );			hObj=hDECOMP->hObj;
	GIP_import_IplImage( hIMG,&hObj->M,&hObj->N,&hObj->data,0x0 );
	M=hObj->M;			N=hObj->N;

//	hDECOMP->alg = CGM_SS_DIRECT | CGM_A_AT | CGM_METHOD_MEDIAN;				//CGM_SS_ITER;	CGM_METHOD_DUAL
	hDECOMP->lenda=1.0;
	hDECOMP->alpha = 30;			hDECOMP->beta = beta = 0.0001;
	MN=M*N;
	t_len=MAX((M+2)*(N+2)*6,N);
	hDECOMP->temp=GIP_alloc( sizeof(double)*t_len );			temp=hDECOMP->temp;
/*	
	hDECOMP->tag=(int*)GIP_alloc( sizeof(int)*(M+2)*(N+2) );
	memset( hDECOMP->tag,0x0,sizeof(int)*(M+2)*(N+2) );
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			hDECOMP->f[pos] = val[pos];
			hDECOMP->u[pos] = val[pos];
			if( c==0 )		BIT_SET( hDECOMP->tag[pos],GIP_BD_LEFT );
			if( c==N-1 )	BIT_SET( hDECOMP->tag[pos],GIP_BD_RIGHT );
			if( r==0 )		BIT_SET( hDECOMP->tag[pos],GIP_BD_DOWN );
			if( r==M-1 )	BIT_SET( hDECOMP->tag[pos],GIP_BD_UP );
		}
	}
*/

GIP_exit:
	;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void GIP_DECOMP_clear( GIP_DECOMP *hDECOMP )	{
	if( hDECOMP == GIP_NULL )	return;

	if( hDECOMP->temp != GIP_NULL )	
	{	GIP_free( hDECOMP->temp );	hDECOMP->temp=GIP_NULL; }
	if( hDECOMP->hObj != GIP_NULL )		{	
		GIP_OBJ_clear( hDECOMP->hObj );	
		GIP_free( hDECOMP->hObj );	hDECOMP->hObj=GIP_NULL; 
	}
}

/*
	A testing code 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/4/2008	
*/
void _DECOMP_SEG_threshold( GIP_DECOMP *hDECOMP,GIP_FLOATING *u,double *temp,int flag )	{
	GIP_OBJECT *hObj=hDECOMP->hObj;
	int r,c,M=hObj->M,N=hObj->N,pos;
	double thresh=10,ux,uy,du;

	for( r=0; r<M; r++ )	{
	for( c=0; c<N; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		ux=r==0 || r==M-1 ? 0.0 : (u[pos+N]-u[pos-N])/2.0;		
		uy=c==0 || c==N-1 ? 0.0 : (u[pos+1]-u[pos-1])/2.0;		
		du=sqrt( ux*ux+uy*uy );
		if( du > thresh )	{
			temp[pos]=0.0;
		}else
			temp[pos]=u[pos];
	}
	}	

	for( r=0; r<M; r++ )	{
	for( c=0; c<N; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		u[pos]=temp[pos];
	}
	}


}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/11/2008	
*/
void _DECOMP_get_v( GIP_DECOMP *hDECOMP,double *u,double *v )	{
	GIP_OBJECT *hObj=hDECOMP->hObj;
	int M=hObj->M,N=hObj->N,r,c,pos,poo=-1;
	double *f=hObj->data;
	char sTitle[80];

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,N );
			v[pos] = f[pos]-u[pos];
		}
	};

	_stprintf( sTitle,".\\trace\\v_%d.bmp",M*N );		
	GIP_save_IplImage_d( M,N,sTitle,v,N,0x0 );
}

/*
	A testing code

	注意：
		1 hIMG的内容返回时被修改
		1 hIMG_0的内容不变，用于校验结果。

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/5/2008	
*/
void GIP_DECOMP_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_0,int model,int loop_max,char* sWndName )	{
	double *temp,res,tol,dT=0.5,s,rmse,snr,v_1,v_2,*v;
	GIP_DECOMP decomp,*hDECOMP=&decomp;
	GIP_FEM_ fem,*hFEM=&fem;	
	GIP_OBJECT *hObj,obj_0;
	GIP_IMAGE *hIMG_core;
	int M,N,isTrace=1,nKey,difus_step=1,i,n_limit=90000;
	char sTitle[80],sModel[80];

	GIP_init_heapinfo( );

	if( hIMG->height*hIMG->width>=n_limit )	{
		s=sqrt(n_limit*1.0/hIMG->height/hIMG->width);
		M=hIMG->height*s;		N=hIMG->width*s;
		M &= -(1<<2);		N &= -(1<<2);
		hIMG_core = GIP_IMG_resize( hIMG_0,M,N,0x0 );
		GIP_import_IplImage( hIMG_core,&(obj_0.M),&(obj_0.N),&(obj_0.data),0x0 );
		cvReleaseImage( &hIMG_core );
		hIMG_core = GIP_IMG_resize( hIMG,M,N,0x0 );
		G_PRINTF( "RESIZE:\tM:%d->%d,N:%d->%d\r\n",hIMG->height,M,hIMG->width,N );		
	}	else	{
		hIMG_core = hIMG;
		GIP_import_IplImage( hIMG_0,&(obj_0.M),&(obj_0.N),&(obj_0.data),0x0 );
	}
	GIP_DECOMP_init_( hIMG_core,hDECOMP );
	hObj=hDECOMP->hObj;
	if( 0 )				//	!!testing!!		
//		_DFT_test_0( hObj->M,hObj->N,hObj->N,hObj->data,0x0 );
		_DFT_test_c( hObj->M,hObj->N,hObj->N,hObj->data,0x0 );
//		_DFT_test( hObj->M,hObj->N,hObj->N,hObj->data,0x0 );
	if( 0 )	{			//	!!testing!!		
		memset( hObj->data,0x0,sizeof(double)*hObj->M*hObj->N );
		hObj->data[0]=0.0;			hObj->data[1]=1.0;
		hObj->data[2]=1.0;			hObj->data[3]=2.0;
	}

	memset( hFEM,0x0,sizeof(GIP_FEM_) );
	M=hObj->M,				N=hObj->N;
	temp=hDECOMP->temp;
	res=hDECOMP->nn_res;
	tol = res*1.0e-4;

	GIP_RES_quality( hObj->M,hObj->N,hObj->N,hObj->data,obj_0.data,&rmse,&snr,&v_1,&v_2 );
	GIP_MODEL_info( sModel,model,0x0 );
	G_PRINTF( "START:\tM=%d,N=%d,model=%s,rmse=%.2g,snr=%.2g,v_1=%.2g,v_2=%.2g\r\n",M,N,sModel,rmse,snr,v_1,v_2 );		
//	memset( hDECOMP->u,0x0,sizeof(double)*M*N );
	GIP_FEM_preprocess( hFEM,hDECOMP->hObj,model,0x0 );			
	GIP_FEM_dump( hFEM,0x0 );	
	while( 1 )	{		//Bregman iteration
		GIP_FEM_core( hFEM,loop_max,0x0 );		//fem solver
		v=hFEM->D_temp;
		_DECOMP_get_v( hDECOMP,hFEM->u,v );
		break;
		for( i = 0; i < M*N; i++ )	{ 
			hObj->data[i]=hObj->data[i]+v[i];
		}
	}
	GIP_RES_quality( hObj->M,hObj->N,hObj->N,hFEM->u,obj_0.data,&rmse,&snr,&v_1,&v_2 );
	GIP_OBJ_clear( &obj_0 );

	G_PRINTF( "\r\nFINISH:\tRMSE=%g,SNR=%g,v_1=%g,v_2=%g\r\n",rmse,snr,v_1,v_2	);

GIP_exit:
	if( 0 )	{		//调用cvPyrSegmentation等
//		_DECOMP_SEG_threshold( hDECOMP,hFEM->u,hFEM->D_temp,0x0 );
		GIP_IMG_segm( hObj->M,hObj->N,hFEM->u,hObj->N,hObj->data,0x0  );
	}
	GIP_export_IplImage_d( hIMG_core,hFEM->u,hObj->N,0x0 );

	if( hIMG_core != hIMG )	{
		cvResize(  hIMG_core,hIMG,CV_INTER_CUBIC );
		cvReleaseImage( &hIMG_core );
	}else	{
	}

	GIP_FEM_clear( hFEM );	
	GIP_DECOMP_clear( hDECOMP );
	GIP_exit_heapinfo( );

}
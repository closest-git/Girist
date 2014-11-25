/*
	v0.1	10/8/2008
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
#include "GIP_diffusion.h"
#include "../fem/GIP_fem.h"

/*
	[REF] A nonlinear primal-dual method for total variation-based image restoration (1999)

	×¢Òâ£º
		1 hDIFUS->ldu=N
		2 AÓëschemeÒªÆ¥Åä

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/17/2008	
*/
void GIP_DIFUS_init_( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,GIP_DIFFUSION *hDIFUS )	{
	int M,N,MN,r,c,pos,i,ret=GIP_OK,ldu,BLK_shift,t_len;
	GIP_FLOATING *val=GIP_NULL;
	double beta,*temp;

	memset( hDIFUS,0x0,sizeof(GIP_DIFFUSION) );
	GIP_import_IplImage( hIMG,&M,&N,&val,0x0 );

//	hDIFUS->alg = CGM_SS_DIRECT | CGM_A_AT | CGM_METHOD_MEDIAN;				//CGM_SS_ITER;	CGM_METHOD_DUAL
	hDIFUS->lenda=1.0;
	hDIFUS->alpha = 30;			hDIFUS->beta = beta = 0.0001;
	hDIFUS->M=M;				hDIFUS->N=N;
	MN=M*N;
	hDIFUS->ldu=N;			hDIFUS->shift=0;
	t_len=MAX((M+2)*(N+2)*6,N);
	hDIFUS->temp=GIP_alloc( sizeof(double)*t_len );			temp=hDIFUS->temp;
	hDIFUS->f=(double*)GIP_alloc( sizeof(double)*(MN) );
	hDIFUS->g=(double*)GIP_alloc( sizeof(double)*(MN) );
	hDIFUS->u=(double*)GIP_alloc( sizeof(double)*MN );
	hDIFUS->tag=(int*)GIP_alloc( sizeof(int)*(M+2)*(N+2) );
	memset( hDIFUS->tag,0x0,sizeof(int)*(M+2)*(N+2) );
	ldu=hDIFUS->ldu;			
	BLK_shift=M*N;
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			hDIFUS->f[pos] = val[pos];
			hDIFUS->u[pos] = val[pos];
			if( c==0 )		BIT_SET( hDIFUS->tag[pos],GIP_BD_LEFT );
			if( c==N-1 )	BIT_SET( hDIFUS->tag[pos],GIP_BD_RIGHT );
			if( r==0 )		BIT_SET( hDIFUS->tag[pos],GIP_BD_DOWN );
			if( r==M-1 )	BIT_SET( hDIFUS->tag[pos],GIP_BD_UP );
		}
	}


GIP_exit:	
	GIP_free( val );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void GIP_DIFUS_clear( GIP_DIFFUSION *hDIFUS )	{
	if( hDIFUS == GIP_NULL )	return;

	hDIFUS->M=-1;		hDIFUS->N=-1;

	if( hDIFUS->f != GIP_NULL )	
	{	GIP_free( hDIFUS->f );	hDIFUS->f=GIP_NULL; }
	if( hDIFUS->u != GIP_NULL )	
	{	GIP_free( hDIFUS->u );	hDIFUS->u=GIP_NULL; }
	if( hDIFUS->g != GIP_NULL )	
	{	GIP_free( hDIFUS->g );	hDIFUS->g=GIP_NULL; }
	if( hDIFUS->tag != GIP_NULL )	
	{	GIP_free( hDIFUS->tag );	hDIFUS->tag=GIP_NULL; }
	if( hDIFUS->temp != GIP_NULL )	
	{	GIP_free( hDIFUS->temp );	hDIFUS->temp=GIP_NULL; }
	if( hDIFUS->temp != GIP_NULL )	
	{	GIP_free( hDIFUS->temp );	hDIFUS->temp=GIP_NULL; }
}

/*	
	energy

	v0.1	cys
		10/14/2008
*/
double  _DIFUS_Energy_( GIP_DIFFUSION *hDIFUS,double *u_cur,int flag )	{
	int r,c,ldu=hDIFUS->ldu,M=hDIFUS->M,N=hDIFUS->N,pos_D;
	double res=0.0,energy=0.0,a,fx,fy,a_max;
	double *f=hDIFUS->f,alpha=hDIFUS->alpha,beta=hDIFUS->beta;

	a_max = 0.0;
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos_D = GIP_RC2POS( r,c,ldu );
			fx = FD_R1_0( c==N-1,GIP_SCHEME_FORWARD,u_cur,pos_D,1 );
			fy = FD_R1_0( r==M-1,GIP_SCHEME_FORWARD,u_cur,pos_D,ldu );
			energy += sqrt(fx*fx+fy*fy+beta)*alpha;
			a = u_cur[pos_D]-f[pos_D];
			a_max = MAX( a_max,fabs(a) );
			energy += 0.5*a*a ;
		}
	}	

	return energy;
}
/*
	A testing code 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/5/2008	
*/
void GIP_DIFFUSION_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName )	{
	double *temp,res,tol,dT=0.5;
	GIP_DIFFUSION difus,*hDIFUS=&difus;
	GIP_FEM_ fem,*hFEM=&fem;	
	int M,N,isTrace=1,nKey,difus_step=1,i;
	char sTitle[80];

	GIP_init_heapinfo( );

	GIP_DIFUS_init_( hIMG,hIMG_mask,hDIFUS );
	if( 0 )	{			//	!!testing!!		
		memset( hDIFUS->u,0x0,sizeof(double)*hDIFUS->M*hDIFUS->ldu );
		hDIFUS->u[0]=0.0;		hDIFUS->u[1]=1.0;
		hDIFUS->u[2]=1.0;		hDIFUS->u[3]=2.0;
	}
	memset( hFEM,0x0,sizeof(GIP_FEM_) );
	M=hDIFUS->M,		N=hDIFUS->N;
	temp=hDIFUS->temp;
	res=hDIFUS->nn_res;
	tol = res*1.0e-4;

	G_PRINTF( "START:\tM=%d,N=%d, energy=%g\r\n",M,N,_DIFUS_Energy_(hDIFUS,hDIFUS->u,0x0) );		
//	memset( hDIFUS->u,0x0,sizeof(double)*M*N );
	GIP_FEM_preprocess( hFEM,hDIFUS->M,hDIFUS->N,hDIFUS->ldu,0x0,hDIFUS->u,GIP_SEG_CV );
	if( 1 )	{		//steady analysis
		GIP_FEM_steady( hFEM );
		GIP_show_IplImage_d( M,N,sWndName,hFEM->u,N,GIP_OPENCV_CREATEWND );
	}else	{
		for( i = 0; i < 10; i++ )	{ 
			GIP_FEM_transient( hFEM, dT );
			_stprintf( sTitle,".\\trace\\u_%d.bmp",difus_step );		
			GIP_save_IplImage_d( M,N,sTitle,hFEM->u,N,0x0 );
			GIP_show_IplImage_d( M,N,sWndName,hFEM->u,N,GIP_OPENCV_CREATEWND );
			G_PRINTF( "DIFUS_%d: energy=%g,|dU|=%g\r\n",difus_step,_DIFUS_Energy_(hDIFUS,hFEM->u,0x0)
				,hFEM->dU );
			difus_step++;
			switch( (nKey= cvWaitKey( )) )	{
				case 27:
					G_PRINTF( "IMAGE diffusion iteration break out\r\n"	);
					GIP_EXIT;
				case '+':
					break;
				case '-':
					break;
				default:
					break;
			}
		}	
	}
	G_PRINTF( "\r\nDIFFUSION test finish.\r\n"	);

GIP_exit:
	GIP_FEM_clear( hFEM );	
	GIP_DIFUS_clear( hDIFUS );
//	GIP_exit_heapinfo( );
}
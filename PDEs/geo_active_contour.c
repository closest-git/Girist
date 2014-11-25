/*
	geodesic active contours for 2-D Images

	1 finite difference scheme
	2 2-D Images(M rows,N columns)	[32] P.13
	3 Kink point
		normal:	pick an arbitrary direction			[11] P.11
		curvature:									[11] P.13
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
#include "geo_active_contour.h"
#include "GIP_fmm.h"

#define GAC_POS( r,c,ld )	((r)*(ld)+(c))

static int gac_step=0;
/*
	
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/27/2008	

int GAC_head_off( GAC_SOLVER *hGAC )	{
	return hGAC->ldu+1;
}
*/

/*
	ghost point of boundary
	Neumann boundary condition: du=0.0

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/27/2008	
*/
void GAC_set_boundary( GAC_SOLVER *hGAC )	{
	int M=hGAC->M,N=hGAC->N,i,ldu=hGAC->ldu;
	double *u_buffer=hGAC->u_buffer;
	for( i = 1; i <= M; i++ )	{
		u_buffer[i] = u_buffer[i+2*ldu];
		u_buffer[M*ldu+i] = u_buffer[i+(M-2)*ldu];
	}
	for( i = 1; i <= N; i++ )	{
		u_buffer[i*ldu] = u_buffer[i*ldu+2];
		u_buffer[M+i*ldu] = u_buffer[(M-2)+i*ldu];
	}
}

/*
	根据u函数确定interface,然后扩充singed distance function到整个定义域[0:M+1,0:N+1]
	注意:
		u函数应尽量光滑，即|Du|=>1.0

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/2/2008	

void GAC_smooth_u_( GAC_SOLVER *hGAC )	{
	int M=hGAC->M,N=hGAC->N,i,r,c,ldu=hGAC->ldu,*tag=hGAC->tag,pos,pos_1;
	int neibor[]={-1,1,-ldu,ldu},shift=hGAC->shift,side,border_flag[]={GIP_BD_LEFT,GIP_BD_RIGHT,GIP_BD_DOWN,GIP_BD_UP};
	double *u_buffer=hGAC->u_buffer;
	GIP_HEAP heap_1,heap_2,*heap;

	GIP_HEAP_init( &heap_1,(M+2)*(N+2) );
	GIP_HEAP_init( &heap_2,(M+2)*(N+2) );

//fast marching method/*	
	for( r = 0; r < M+2; r++ )	{	//initialization fast marching method	
		for( c = 0; c < N+2; c++ )	{
			pos = GAC_POS(r,c,ldu);	
			if(	BIT_TEST(tag[pos],GIP_GAC_INTERFACE) )			{	
				BIT_SET(tag[pos],GIP_FMM_KNOWN);	
				for( i = GIP_NEIBOR_LEFT; i <= GIP_NEIBOR_UP; i++ )	{
					if( BIT_TEST(tag[pos],border_flag[i]) )
						continue;
					pos_1 = pos + neibor[i];
					if( BIT_TEST(tag[pos_1],GIP_GAC_INTERFACE) || BIT_TEST(tag[pos_1],GIP_FMM_BAND) )	
						continue;
					if( u_buffer[pos_1]>0.0 )
					{	side = GIP_FMM_SIDE_2;	heap=&heap_2;	}
					if( u_buffer[pos_1]<0.0 )
					{	side = GIP_FMM_SIDE_1;	heap=&heap_1;	}
					BIT_SET( tag[pos_1],side );
					BIT_SET(tag[pos_1],GIP_FMM_BAND);
					if( u_buffer[pos_1] < 0.0 )		u_buffer[pos_1] *= -1;
					GIP_HEAP_insert( heap,pos_1,u_buffer[pos_1] );
				}
			}
		}
	}
	for( r = 0; r < M+2; r++ )	{	//initialization fast marching method	
	for( c = 0; c < N+2; c++ )	{
		pos = GAC_POS(r,c,ldu);	
		if( BIT_TEST(tag[pos],GIP_FMM_KNOWN) )
			continue;
		else if( BIT_TEST(tag[pos],GIP_FMM_BAND)  )	{
			continue;
		}else	{
			if( u_buffer[pos] < 0.0 )
			{	BIT_SET( tag[pos],GIP_FMM_SIDE_1 );			u_buffer[pos]=1.0e6;	}
			else if ( u_buffer[pos] > 0.0 )
			{	BIT_SET( tag[pos],GIP_FMM_SIDE_2 );			u_buffer[pos]=1.0e6;	}
		}
	}
	}

	GIP_FMM_update( M,N,hGAC->u_buffer,hGAC->tag,hGAC->ldu,&heap_2,GIP_FMM_SIDE_2,0x0 );
	GIP_FMM_update( M,N,hGAC->u_buffer,hGAC->tag,hGAC->ldu,&heap_1,GIP_FMM_SIDE_1,0x0 );
	GIP_HEAP_clear( &heap_1 );
	GIP_HEAP_clear( &heap_2 );
	for( r = 0; r < M+2; r++ )	{	//initialization fast marching method	
	for( c = 0; c < N+2; c++ )	{
		pos = GAC_POS(r,c,ldu);	
		if( BIT_TEST( tag[pos],GIP_FMM_SIDE_1 ) )
		{	u_buffer[pos] *= -1;	}		
	}
	}
//verify |delta U|=1.0

	if( 1 )	{
		double average=0.0,delta,ux,uy,*u=hGAC->u_buffer+shift;
		for( r = 0; r < M-1; r++ )	{
			for( c = 1; c < N-1; c++ )	{
				ux=(u[c+1]-u[c-1])/2;
				uy=(u[c+ldu]-u[c-ldu])/2;
				delta = sqrt(ux*ux+uy*uy);
			 	average += delta;
			}
			u+=ldu;
		}
		average /=(M*N);
		G_PRINTF( "Average |detaU|=%g\r\n",average );
	}
	if( 0 )	{
		char sTitle[80];
		_stprintf( sTitle,"u_%d",++g_img_no );
		GIP_show_IplImage_d( M,N,sTitle,(hGAC->u_buffer+hGAC->shift),hGAC->ldu,0x0 );
	}
}
*/
/*
	|Du|:	norm of hamilton operator on U. 

	ISSUE:
		check the small detaU

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/23/2008	
*/
double GAC_norm_Hu_( int M,int N,double *u,int ldu,double *d )	{
	int r,c,nKink=0,band_nz=0;
	double average=0.0,delta,ux,uy,d_max=0.0,d_min=DBL_MAX,*u_0=u;
	double band_width=2,band_average=0.0;

	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			ux=(u[c+1]-u[c-1])/2;
			uy=(u[c+ldu]-u[c-ldu])/2;
			delta = sqrt(ux*ux+uy*uy);
			if( delta < GIP_ZERO_ROUNDOFF )	{		//kink point:  
				nKink++;
			}else	{
				d_max=MAX( d_max,delta );
				d_min=MIN( d_min,delta );
			}
			if( d !=0x0 )	{
				d[GAC_POS( r,c,ldu )] = delta;
			}
			if( fabs(u[c])< band_width )	{
				band_nz++;
				band_average += delta;
			}

			average += delta;
		}
		u+=ldu;
	}
	u=u_0;
	average /=(M*N-nKink);
	band_average /= band_nz;

	G_PRINTF( "%d:	|detaU|	Band Average=%g,Average=%g,max=%g,min=%g,kink=%g\r\n"
		,gac_step,band_average,average,d_max,d_min,nKink*1.0/(M*N) );

	return average;
}

/*
	根据u函数确定interface,然后扩充singed distance function到整个定义域[0:M+1,0:N+1]
	whether a point is inside or outside is immediately apparent from the level set functoin

	注意:
		u函数应尽量光滑，即|Du|=>1.0

	Copyright 2008-present, Grusoft.
	v0.2	cys
		4/22/2008	
*/
void GAC_smooth_u_( GAC_SOLVER *hGAC )	{
	int M=hGAC->M,N=hGAC->N,i,r,c,ldu=hGAC->ldu,*tag=hGAC->tag,pos,pos_1;
	int neibor[]={-1,1,-ldu,ldu},shift=hGAC->shift,side,border_flag[]={GIP_BD_LEFT,GIP_BD_RIGHT,GIP_BD_DOWN,GIP_BD_UP};
	double *u_buffer=hGAC->u_buffer;
	GIP_HEAP heap_1,heap_2,*heap;

	GIP_HEAP_init( &heap_1,(M+2)*(N+2) );
	GIP_HEAP_init( &heap_2,(M+2)*(N+2) );

//fast marching method/*	
	for( r = 0; r < M+2; r++ )	{	//initialization fast marching method	
		for( c = 0; c < N+2; c++ )	{
			pos = GAC_POS(r,c,ldu);	
			if(	BIT_TEST(tag[pos],GIP_GAC_INTERFACE) )			{	
				BIT_SET(tag[pos],GIP_FMM_KNOWN);	
				GIP_HEAP_insert( &heap_1,pos,u_buffer[pos] );
				GIP_HEAP_insert( &heap_2,pos,u_buffer[pos] );
				for( i = GIP_NEIBOR_LEFT; i <= GIP_NEIBOR_UP; i++ )	{
					if( BIT_TEST(tag[pos],border_flag[i]) )
						continue;
					pos_1 = pos + neibor[i];
					if( BIT_TEST(tag[pos_1],GIP_GAC_INTERFACE) || BIT_TEST(tag[pos_1],GIP_FMM_BAND) )	
						continue;
					if( u_buffer[pos_1]>0.0 )
					{	side = GIP_FMM_SIDE_2;	heap=&heap_2;	}
					if( u_buffer[pos_1]<0.0 )
					{	side = GIP_FMM_SIDE_1;	heap=&heap_1;	}
					BIT_SET( tag[pos_1],side );
			//		BIT_SET(tag[pos_1],GIP_FMM_BAND);
			//		if( u_buffer[pos_1] < 0.0 )		u_buffer[pos_1] *= -1;
				}
			}
		}
	}
	for( r = 0; r < M+2; r++ )	{	//initialization fast marching method	
	for( c = 0; c < N+2; c++ )	{
		pos = GAC_POS(r,c,ldu);	
		if( BIT_TEST(tag[pos],GIP_FMM_KNOWN) )
			continue;
		else if( BIT_TEST(tag[pos],GIP_FMM_BAND)  )	{
			continue;
		}else	{
			if( u_buffer[pos] < 0.0 )
			{	BIT_SET( tag[pos],GIP_FMM_SIDE_1 );			u_buffer[pos]=1.0e6;	}
			else if ( u_buffer[pos] > 0.0 )
			{	BIT_SET( tag[pos],GIP_FMM_SIDE_2 );			u_buffer[pos]=1.0e6;	}
		}
	}
	}

	GIP_FMM_update( M,N,hGAC->u_buffer,hGAC->tag,hGAC->ldu,&heap_1,GIP_FMM_SIDE_1,0x0 );
	GIP_FMM_update( M,N,hGAC->u_buffer,hGAC->tag,hGAC->ldu,&heap_2,GIP_FMM_SIDE_2,0x0 );
	GIP_HEAP_clear( &heap_1 );
	GIP_HEAP_clear( &heap_2 );
	for( r = 0; r < M+2; r++ )	{	//initialization fast marching method	
	for( c = 0; c < N+2; c++ )	{
		pos = GAC_POS(r,c,ldu);	
		if( BIT_TEST( tag[pos],GIP_FMM_SIDE_1 ) )
		{	u_buffer[pos] *= -1;	}		
	}
	}
//	GAC_norm_Hu_( M,N,hGAC->u_buffer+hGAC->shift,hGAC->ldu );		//verify |delta U|=1.0
}

/*
	u_buffer定义为[0:M+1,0:N+1],其中u0为矩形{r=1,r=M,c=1,c=N}

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/2/2008	
*/
void GAC_Init_u0( GAC_SOLVER *hGAC,GIP_IMAGE *hMask,int type )	{
	int M=hGAC->M,N=hGAC->N,r,c,pos,isOut,ldu=hGAC->ldu;
	int r0=1,r1=M,c0=1,c1=N;		//
	double *u_buffer=hGAC->u_buffer,ar,ac,dist;

	if( hMask==0x0 )	{
		for( r = 0; r < M+2; r++ )	{
			for( c = 0; c < N+2; c++ )	{
				pos = GAC_POS( r,c,ldu );
				isOut = 0;
				ar = MIN( fabs(r-r0),fabs(r-r1) );
				ac = MIN( fabs(c-c0),fabs(c-c1) );
				if( r<r0 || r>r1 )
				{	isOut=1;	dist = ar;	}
				else if( c<c0 || c>c1 )
				{	isOut=1;	dist = ac;	}
				else
					dist = MIN(ar,ac);
				u_buffer[pos] = (isOut==1) ? dist : -dist;
				if( u_buffer[pos]==0.0 )			{	
					BIT_SET( hGAC->tag[pos],GIP_GAC_INTERFACE );	
				}
			}
		}
	}else	{		//based on mask imag
		int height=hMask->height,width=hMask->width,step=hMask->widthStep/sizeof(uchar);
		uchar *i_data = (uchar *)hMask->imageData,a;
		ASSERT( height==M && width==N && hMask->depth==IPL_DEPTH_8U && hMask->nChannels==1 );
		for( r = 0; r < M+2; r++ )	{		//暂定点全部在边界之外
			for( c = 0; c < N+2; c++ )	{
				pos = GAC_POS( r,c,ldu );
				u_buffer[pos]=1.0;				
			}
		}
		for( r = 0; r < height; r++ )	{
			for( c = 0; c < width; c++ )	{
				a = i_data[r*step+c];
				pos = GAC_POS( r+1,c+1,ldu );	
				if( a == 0 )	{
					u_buffer[pos]=0.0;		
					BIT_SET( hGAC->tag[pos],GIP_GAC_INTERFACE );
				}else if( a != UCHAR_MAX )	{		//边界内部的点
					u_buffer[pos]=-1.0;
				}
			}
		}
			
	}

}

/*
	computes the reinitialized distance function in a band 
	that is one-layer thick on each side of the interface.
	[REF]: The Fast Construction of Extension Velocities in Level Set Methods, D. Adalsteinsson and J. A. Sethian, J. Comput. Phys. 148, 2 (1999).
	[REF]: A Study of Numerical Methods for the Level Set Approach
	
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/21/2008	
	v0.2	cys
		5/6/2008	
*/
void GAC_ReInit_u0( GAC_SOLVER *hGAC,int type )	{
	int M=hGAC->M,N=hGAC->N,r,c,pos,isOut,ldu=hGAC->ldu,*tag=hGAC->tag,fix;
	int neibor[]={-1,1,-ldu,ldu},nCross,i,sdf_alg=GIP_SDF_INTERSECTION;
	double *u_buffer=hGAC->u_buffer,a,b,cross[4],s,t;

	for( r = 0; r < M+2; r++ )	{
	for( c = 0; c < N+2; c++ )	{
		pos = GAC_POS( r,c,ldu );
		tag[pos] = tag[pos]&TAG_INIT_FIX;
	}
	}

	for( r = 1; r < M+1; r++ )	{
		for( c = 1; c < N+1; c++ )	{
			pos = GAC_POS( r,c,ldu );
			a = u_buffer[pos];
			if( a==0.0 )	{
				BIT_SET( tag[pos],GIP_GAC_INTERFACE );	
			}else {
				nCross = 0;			s=DBL_MAX;		t=DBL_MAX;
				for( i = GIP_NEIBOR_LEFT; i <= GIP_NEIBOR_UP; i++ )	{
					b = u_buffer[pos+neibor[i]];
					if( a*b<0.0 )			
					{	cross[i]=fabs(a)/( fabs(a)+fabs(b) );	nCross++;	}
					else
						cross[i]=DBL_MAX;
				}
				if( nCross == 0 )	
					continue;
				BIT_SET( tag[pos],GIP_GAC_INTERFACE );	
				switch( sdf_alg )	{
				case GIP_SDF_INTERSECTION:
					s = MIN(cross[0],cross[1]);			t = MIN(cross[2],cross[3]);		//X,Y	
					if( s>1.0 || t>1.0 )	{	
						ASSERT( nCross==1 || nCross==2 );		//Figure 4b 4c 4e
						u_buffer[pos] =	MIN(s,t);		//Figure 4a 4d
					}	else	{
						ASSERT( s<1.0 && t < 1.0 );		//Figure 4b 4c 4e
						u_buffer[pos] =	s*t/sqrt(s*s+t*t);
					}
				break;
				case GIP_SDF_PROJECT_1:	{
					double *u,ux,uy,uxx,uyy,uxy,b,c,alpha;
					u = u_buffer+pos;
					ux=(u[1]-u[-1])/2;				uy=(u[ldu]-u[-ldu])/2;
					uxx=(u[1]-2*u[0]+u[-1]);		uyy=(u[ldu]-2*u[0]+u[-ldu]);
					uxy=((u[1+ldu]-u[1-ldu])/2-(u[-1+ldu]-u[-1-ldu])/2)/2;	
					b = ux*ux+uy*uy;
					c = (ux*ux*uxx+2*ux*uy*uxy+uy*uy*uyy)/2.0;
					ASSERT( (b*b-4*a*c)>=0.0 );
					alpha = (-b+sqrt(b*b-4*a*c))/2.0/c;
					u_buffer[pos] = fabs(alpha*sqrt(b));										}
				break;
				}
				}
				u_buffer[pos] = a > 0.0 ? u_buffer[pos] : -u_buffer[pos];	
//				ASSERT( -1.0 <= u_buffer[pos] && u_buffer[pos]<= 1.0 );
		}
	}
}

/*
	Init edge dector function from Intensity Image
	
	1	edge dector function    from [10] P.153
	2	边界处理：暂时空出

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/19/2008	
*/
void GAC_init( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,GAC_SOLVER *hGAC )	{
	int M,N,MN,r,c,pos;
	GIP_FLOATING *val=GIP_NULL;
	double *g,a;

	gac_step=0;

	memset( hGAC,0x0,sizeof(GAC_SOLVER) );
	GIP_import_IplImage( hIMG,&M,&N,&val,0x0 );

	hGAC->M=M;		hGAC->N=N;
	MN=M*N;
	hGAC->g = GIP_alloc( sizeof(double)*MN );		g=hGAC->g;
	hGAC->gx = GIP_alloc( sizeof(double)*MN );
	hGAC->gy = GIP_alloc( sizeof(double)*MN );

	linear_spatial_filters( M,N,val,g,0 );
	GIP_save_IplImage_d( M,N,"g.bmp",g,N,GIP_OPENCV_CREATEWND );
	memset( hGAC->gx,0x0,sizeof(double)*MN );
	memset( hGAC->gy,0x0,sizeof(double)*MN );
	for( r = 1; r < M-1; r++ )	{
	for( c = 1; c < N-1; c++ )	{
		pos=GAC_POS( r,c,N );		//需被替换为GAC_POS
		a = g[pos];
		hGAC->gx[pos] = (a-g[pos-1]);		hGAC->gy[pos] = (a-g[pos-N]);
		hGAC->g_max=MAX( hGAC->g_max,fabs(a) );
		hGAC->gx_max=MAX( hGAC->gx_max,fabs(hGAC->gx[pos]) );
		hGAC->gy_max=MAX( hGAC->gy_max,fabs(hGAC->gy[pos]) );
	}
	}
	GIP_save_IplImage_d( M,N,"gx.bmp",hGAC->gx,N,GIP_OPENCV_CREATEWND );
	GIP_save_IplImage_d( M,N,"gy.bmp",hGAC->gy,N,GIP_OPENCV_CREATEWND );
//	GIP_show_IplImage_d( M,N,"g",g,N,GIP_OPENCV_CREATEWND );
//	GIP_show_IplImage_d( M,N,"gx",hGAC->gx,N,GIP_OPENCV_CREATEWND );
//	GIP_show_IplImage_d( M,N,"gy",hGAC->gy,N,GIP_OPENCV_CREATEWND );

	hGAC->ldu=N+2;			hGAC->shift=N+3;
	hGAC->u_buffer=(double*)GIP_alloc( sizeof(double)*(M+2)*(N+2) );

	hGAC->tag=(int*)GIP_alloc( sizeof(int)*(M+2)*(N+2) );
	memset( hGAC->tag,0x0,sizeof(int)*(M+2)*(N+2) );
	for( r = 0; r < M+2; r++ )	{
		for( c = 0; c < N+2; c++ )	{
			pos = GAC_POS( r,c,hGAC->ldu );
			if( c==0 )		BIT_SET( hGAC->tag[pos],GIP_BD_LEFT );
			if( c==N+1 )	BIT_SET( hGAC->tag[pos],GIP_BD_RIGHT );
			if( r==0 )		BIT_SET( hGAC->tag[pos],GIP_BD_DOWN );
			if( r==M+1 )	BIT_SET( hGAC->tag[pos],GIP_BD_UP );
		}
	}
	GAC_Init_u0( hGAC,hIMG_mask,0x0 );
	GAC_smooth_u_( hGAC );
GIP_exit:	
	GIP_free( val );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void GAC_clear( GAC_SOLVER *hGAC )	{
	if( hGAC == GIP_NULL )	return;

	hGAC->M=-1;		hGAC->N=-1;
	if( hGAC->g != GIP_NULL )	
	{	GIP_free( hGAC->g );	hGAC->g=GIP_NULL; }
	if( hGAC->gx != GIP_NULL )	
	{	GIP_free( hGAC->gx );	hGAC->gx=GIP_NULL; }
	if( hGAC->gy != GIP_NULL )	
	{	GIP_free( hGAC->gy );	hGAC->gy=GIP_NULL; }
	if( hGAC->u_buffer != GIP_NULL )	
	{	GIP_free( hGAC->u_buffer );	hGAC->u_buffer=GIP_NULL; }
	if( hGAC->tag != GIP_NULL )	
	{	GIP_free( hGAC->tag );	hGAC->tag=GIP_NULL; }
}

/*
	bandwith is typically 5 or 6
*/

/*
	discrete scheme from [10] P.264
	1 based on signed distance function

	ISSUE-NEEDED:
		1、即时更新u值，会出现什么情况？
		2、normH==1.0时，K!=(uxx+uyy)?

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/19/2008	
*/
void GAC_step( GAC_SOLVER *hGAC,double dT,double *temp )	{
	int M=hGAC->M,N=hGAC->N,MN=M*N,i,r,c,g_pos,pos,ldu=hGAC->ldu,*tag=hGAC->tag;
	double mc_motion,advection,du,K,dU,K_min=DBL_MAX,K_max=-K_min;
	double *u,*gx=hGAC->gx,*gy=hGAC->gy,ux,uy,uxx,uyy,uxy,normH;

	for( r = 0; r < M; r++ )	{
		u = hGAC->u_buffer+hGAC->shift+r*ldu;
		for( c = 0; c < N; c++ )	{
			pos=GAC_POS(r+1,c+1,hGAC->ldu);			ASSERT( pos==u-hGAC->u_buffer );
			g_pos=GAC_POS(r,c,N);	

			uxx=(u[1]-2*u[0]+u[-1]);		uyy=(u[ldu]-2*u[0]+u[-ldu]);
			ux=(u[1]-u[-1])/2;				uy=(u[ldu]-u[-ldu])/2;
			normH=sqrt(ux*ux+uy*uy);		//norm of hamilton operator. 
			if( normH == 0.0 )	{			//kink point
				K = 1.0;
			}else	{
				uxy = ((u[1+ldu]-u[1-ldu])/2-(u[-1+ldu]-u[-1-ldu])/2)/2;	
				K=(ux*ux*uyy-2*ux*uy*uxy+uy*uy*uxx)/(normH*normH*normH);
/*				K=(uxx+uyy);	
				if( normH==1.0 )	
					if( K>0.0)	ASSERT( 1.0-GIP_ZERO_ROUNDOFF<K/(uxx+uyy) && K/(uxx+uyy)<1.0+GIP_ZERO_ROUNDOFF );			//only valid for singed distance function
				else	{
				}*/
			}
			K_max=MAX( K_max,K );		K_min=MIN( K_min,K );
//			ASSERT( -1.0<=K && K<=1.0 );
			mc_motion = K*normH;	
			advection=0.0;		du=0.0;
			if( gx[g_pos] < 0.0 )
			{	du=(u[1]-u[0]);	}
			else if( gx[g_pos] > 0.0 )
			{	du=(u[0]-u[-1]);	}
			advection+=gx[g_pos]*du;
			if( gy[g_pos] < 0.0 )
			{	du=(u[ldu]-u[0]);	}
			else if( gy[g_pos] > 0.0 )
			{	du=(u[0]-u[-ldu]);	}
			advection+=gy[g_pos]*du;
			dU = (hGAC->g[g_pos]*mc_motion+advection);
		//	dU = mc_motion;			//mean curvature motion测试
		//	dU = advection;			//advection测试
			temp[pos] = dU*dT;
		//	*u += dU*dT;
			u++;
		}
	}

	u=hGAC->u_buffer;
	for( r = 0; r < M; r++ )	{
//		u = hGAC->u_buffer+hGAC->shift+r*ldu;
		for( c = 0; c < N; c++ )	{
			pos=GAC_POS(r+1,c+1,hGAC->ldu);			
			u[pos] += temp[pos];
		}
	}

	G_PRINTF( "\tK_max=%g,K_min=%g\r\n",K_max,K_min );

//	for( i = 0; i < MN; i++ )	u[i] += temp[i]*dT;


}

/*
	A basic code 
	1、stability condition 
		[11] P.44 (4.7)	mean curvature motion : dT<0.5		

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void GAC_basic( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName )	{
	double t=0,dT=0.25,t_final=dT*1000,average;
	double *temp;
	GAC_SOLVER gac,*hGAC=&gac;
	int M,N,t_len,isTrace=1,nKey,t_step=1;
	char sTitle[80];

	GAC_init( hIMG,hIMG_mask,&gac );
	dT=1.0/(hGAC->gx_max+hGAC->gy_max+hGAC->g_max*4);

	M=gac.M,		N=gac.N;
	t_len=MAX((M+2)*(N+2),N);
	temp=GIP_alloc( sizeof(double)*t_len );
//	GIP_CONTOUR_1( M,N,gac.u_buffer+gac.shift,temp+gac.shift,gac.ldu,0x0 );
//	GIP_export_IplImage_d( hIMG,(temp+gac.shift),gac.ldu,0x0 );
//	goto GIP_exit;
	while( t < t_final )	{
		if( isTrace==1 && gac_step%t_step==0 )	{
			_stprintf( sTitle,".\\trace\\u_%d.bmp",gac_step );		
			GIP_save_IplImage_d( M,N,sTitle,(gac.u_buffer+gac.shift),gac.ldu,0x0 );
//			_stprintf( sTitle,"u_%d.bmp",gac_step );		
//			GIP_show_IplImage_d( M,N,sTitle,(gac.u_buffer+gac.shift),gac.ldu,GIP_OPENCV_CREATEWND );
//			_stprintf( sTitle,".\\trace\\interface_%d.bmp",step );		
			GIP_CONTOUR_1( M,N,gac.u_buffer+gac.shift,temp+gac.shift,gac.ldu,0x0 );
			GIP_show_IplImage_d( M,N,sWndName,(temp+gac.shift),gac.ldu,0x0 );
			nKey = cvWaitKey( );
			switch( nKey )	{
				case 27:
					G_PRINTF( "GAC iteration break out\r\n"	);
					GIP_EXIT;
				case '+':
					t_step *= 10;
					break;
				case '-':
					t_step /= 10;
					break;
				default:
					break;
			}
//			GIP_save_IplImage_d( M,N,sTitle,(temp+gac.shift),gac.ldu,0x0 );
		}
		average = GAC_norm_Hu_( M,N,hGAC->u_buffer+hGAC->shift,hGAC->ldu,temp+hGAC->shift );		
		_stprintf( sTitle,".\\trace\\d_%d.bmp",gac_step );		
		GIP_save_IplImage_d( M,N,sTitle,(temp+gac.shift),gac.ldu,0x0 );
		GAC_step( &gac,dT,temp );
 		if( BIT_TEST(gac.flag,GAC_SOLVER_CONVERGE) )
			break;
		gac_step++;
		//if( average > 2 )	{
			GAC_ReInit_u0( &gac,0x0 );
			GAC_smooth_u_( &gac );
		//}
		t += dT;
	}
	G_PRINTF( "GAC iteration finish\r\n"	);
	GIP_free( temp );
GIP_exit:
	GAC_clear( &gac );
}
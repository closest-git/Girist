#include <math.h>
#include "../pch/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/GIP_util.h"
#include "../util/GIP_heap.h"
#include "GIP_fmm.h"

#define ZERO_ROUND_OFF 1.0E-9

/*
	一些方向
	1、The specific causlity relationship:	why information propagate from smaller values to larger values?
	2、higher-order accuracy	[11] P.73	[12] P.96
	3、Osher提到：Rapid And Accurate Computation Of The Distance Function Using Grids
		also a good introduction for viscosity solution
	4、Huygen's principle and dijkstra's method
*/

/*
	Selecting the largest possible solution to the quadratic equation
		(max(dxT,0))^2+(max(dyT,0))^2=1.0     
		公式参见[13] P.3;	cvinpaint.cpp也采用了该公式，[14]中的公式有误
		

	Ref:	
		1 An Image Inpainting Technique Based on the Fast Marching Method.
		2 cvinpaint.cpp	line 233

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/5/2008	
*/
FMM_FLOAT _fmm_solve( FMM_FLOAT *T,int *flag,int pos_1,int pos_2 )	{
	FMM_FLOAT sol=1.0e6,T_1,T_2,r,s;

	T_1=T[pos_1];		T_2=T[pos_2];		
	if( BIT_TEST( flag[pos_1],GIP_FMM_KNOWN) )	{
		if( BIT_TEST( flag[pos_2],GIP_FMM_KNOWN) )	{
//			ASSERT( 2.0-(T_1-T_2)*(T_1-T_2) >= 0.0 );
			if( fabs(T_1-T_2) > 1.0 )	{
				sol = MIN(T_1,T_2) + 1.0;
			}else	{
				r = sqrt(2.0-(T_1-T_2)*(T_1-T_2) );
				s = (T_1+T_2-r)/2.0;
				sol = s+r;
			/*	if( s>=T_1 && s>=T_2 )
					sol = s;
				else	{
					sol=s+r;
				}*/
			}
		}
		else
			sol = 1+T_1;			
	}else if( BIT_TEST( flag[pos_2],GIP_FMM_KNOWN) ){
		sol = 1+T_2;
	}else{
//		ASSERT( 0 );
//		GIP_ERROR( "_fmm_solve",-1, "Both T_1 and T_2 are unknown" );;
		sol=1.0+min(T_1,T_2);
	}
//	if( BIT_TEST( flag[pos_1],GIP_FMM_KNOWN) )	ASSERT( T_1<=sol+ZERO_ROUND_OFF && sol<=T_1+1.0+ZERO_ROUND_OFF );
//	if( BIT_TEST( flag[pos_2],GIP_FMM_KNOWN) )	ASSERT( T_2<=sol+ZERO_ROUND_OFF && sol<=T_2+1.0+ZERO_ROUND_OFF );
GIP_exit:
	return sol;
}

/*
	FAST MARCHING METHOD
	[ref]	cvinpaint.cpp(629;icvNSInpaintFMM) 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/7/2008	
*/
int GIP_FMM_update( int M,int N,FMM_FLOAT *T,int *tag,int ldr,GIP_HEAP *heap,int side,int flag )	{
	int nzBand=0,pos,pos_1,i,ret=GIP_OK;
	int INP_neibor[]={-1,1,-ldr,ldr},ta_,border_flag[]={GIP_BD_LEFT,GIP_BD_RIGHT,GIP_BD_DOWN,GIP_BD_UP};
	FMM_FLOAT a,a1,a2,a3,a4;

	nzBand = heap->size;
/*	pos = 0;
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			if( BIT_TEST(flag[pos],GIP_FMM_BAND) )		{	
				GIP_HEAP_insert( &heap,pos_1,0.0 );
				nzBand++;	
			}else	{
				BIT_SET(flag[pos],GIP_FMM_INSIDE);
				T[pos]=GIP_FMM_OO;
			}
		}
		pos += ldr;
	}*/

	while( GIP_HEAP_pop( heap,&pos,&a )!=GIP_HEAP_EMPTY )	{
		ASSERT( a==T[pos] );
		nzBand--;
		BIT_SET( tag[pos],GIP_FMM_KNOWN );
		for( i = GIP_NEIBOR_LEFT; i <= GIP_NEIBOR_UP; i++ )	{
			if( BIT_TEST(tag[pos],border_flag[i]) )
				continue;
			pos_1 = pos + INP_neibor[i];
			ta_=tag[pos_1];
			if( BIT_TEST(ta_,GIP_FMM_KNOWN) )	continue;
//			if( BIT_TEST(ta_,GIP_FMM_BAND) )		continue;
			if( BIT_TEST(ta_,side) )		{
				a1 = ( BIT_TEST(ta_,GIP_BD_LEFT) || BIT_TEST(ta_,GIP_BD_DOWN) ) ? DBL_MAX :
					_fmm_solve( T,tag,pos_1+INP_neibor[GIP_NEIBOR_LEFT],pos_1+INP_neibor[GIP_NEIBOR_DOWN] );
				a2 = ( BIT_TEST(ta_,GIP_BD_RIGHT) || BIT_TEST(ta_,GIP_BD_DOWN) ) ? DBL_MAX :
					_fmm_solve( T,tag,pos_1+INP_neibor[GIP_NEIBOR_RIGHT],pos_1+INP_neibor[GIP_NEIBOR_DOWN] );
				a3 = ( BIT_TEST(ta_,GIP_BD_LEFT) || BIT_TEST(ta_,GIP_BD_UP) ) ? DBL_MAX :
					_fmm_solve( T,tag,pos_1+INP_neibor[GIP_NEIBOR_LEFT],pos_1+INP_neibor[GIP_NEIBOR_UP] );
				a4 = ( BIT_TEST(ta_,GIP_BD_RIGHT) || BIT_TEST(ta_,GIP_BD_UP) ) ? DBL_MAX :
					_fmm_solve( T,tag,pos_1+INP_neibor[GIP_NEIBOR_RIGHT],pos_1+INP_neibor[GIP_NEIBOR_UP] );
				a = MIN( MIN(a1,a2),MIN(a3,a4) );
				if( a == DBL_MAX )
					GIP_ERROR( "GIP_FMM_update",GIP_ERR_FMMVAL,"Strange value from fmm_solve" );
				ASSERT( a <= T[pos_1]+GIP_ZERO_ROUNDOFF );		//[11] P.71 tentative values should only decrease	
				if( BIT_TEST(ta_,GIP_FMM_BAND) )	{
					if( a < T[pos_1]-GIP_ZERO_ROUNDOFF )	{
						T[pos_1] = a;
						GIP_HEAP_update( heap,pos_1,a );
					}
				}else	{
					T[pos_1] = a;
					GIP_HEAP_insert( heap,pos_1,T[pos_1] );
					BIT_SET( tag[pos_1],GIP_FMM_BAND );
				}
				nzBand++;
			}
		}
	}
GIP_exit:	
	return ret;
}


/*
	基于FMM来光滑u函数，即|Du|=>1.0

	Copyright 2008-present, Grusoft.
	v0.1	cys
		5/30/2008	
*/
void GIP_FMM_smooth_u_( int M,int N,double *u,int ldu,int *tag,int flag )	{
	int i,r,c,pos,pos_1;
	int neibor[]={-1,1,-ldu,ldu},side,border_flag[]={GIP_BD_LEFT,GIP_BD_RIGHT,GIP_BD_DOWN,GIP_BD_UP};
	GIP_HEAP heap_1,heap_2,*heap;

	GIP_HEAP_init( &heap_1,M*N );
	GIP_HEAP_init( &heap_2,M*N );

//fast marching method/*	
	for( r = 0; r < M; r++ )	{	//initialization fast marching method	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS(r,c,ldu);	
			if(	BIT_TEST(tag[pos],GIP_INTERFACE) )			{	
				BIT_SET(tag[pos],GIP_FMM_KNOWN);	
				GIP_HEAP_insert( &heap_1,pos,u[pos] );
				GIP_HEAP_insert( &heap_2,pos,u[pos] );
				for( i = GIP_NEIBOR_LEFT; i <= GIP_NEIBOR_UP; i++ )	{
					if( BIT_TEST(tag[pos],border_flag[i]) )
						continue;
					pos_1 = pos + neibor[i];
					if( BIT_TEST(tag[pos_1],GIP_INTERFACE) || BIT_TEST(tag[pos_1],GIP_FMM_BAND) )	
						continue;
					if( u[pos_1]>0.0 )
					{	side = GIP_FMM_SIDE_2;	heap=&heap_2;	}
					if( u[pos_1]<0.0 )
					{	side = GIP_FMM_SIDE_1;	heap=&heap_1;	}
					BIT_SET( tag[pos_1],side );
			//		BIT_SET(tag[pos_1],GIP_FMM_BAND);
			//		if( u[pos_1] < 0.0 )		u[pos_1] *= -1;
				}
			}
		}
	}
	for( r = 0; r < M; r++ )	{	//initialization fast marching method	
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS(r,c,ldu);	
		if( BIT_TEST(tag[pos],GIP_FMM_KNOWN) )
			continue;
		else if( BIT_TEST(tag[pos],GIP_FMM_BAND)  )	{
			continue;
		}else	{
			if( u[pos] < 0.0 )
			{	BIT_SET( tag[pos],GIP_FMM_SIDE_1 );			u[pos]=1.0e6;	}
			else if ( u[pos] > 0.0 )
			{	BIT_SET( tag[pos],GIP_FMM_SIDE_2 );			u[pos]=1.0e6;	}
		}
	}
	}

	GIP_FMM_update( M,N,u,tag,ldu,&heap_1,GIP_FMM_SIDE_1,0x0 );
	GIP_FMM_update( M,N,u,tag,ldu,&heap_2,GIP_FMM_SIDE_2,0x0 );
	GIP_HEAP_clear( &heap_1 );
	GIP_HEAP_clear( &heap_2 );
	for( r = 0; r < M; r++ )	{	//initialization fast marching method	
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS(r,c,ldu);	
		if( BIT_TEST( tag[pos],GIP_FMM_SIDE_1 ) )
		{	u[pos] *= -1;	}		
	}
	}

}
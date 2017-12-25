#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include "../PCH/GIP_def.h"
#include "../util/gip_util.h"
#include "GIP_contour.h"
/*
	暂时以稀疏矩阵格式存储 contour(level set)
*/

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void GIP_CONTOUR_init( GIP_CONTOUR *hContour,int M,int N,int nz )	{
	hContour->M=M;			hContour->N=N;
	hContour->nz=nz;
	if( hContour->ptr !=GIP_NULL )		GIP_free( hContour->ptr );
	hContour->ptr=GIP_alloc( sizeof(int)*(N+1) );
	if( hContour->ind !=GIP_NULL )		GIP_free( hContour->ind );
	hContour->ind=GIP_alloc( sizeof(int)*nz );
	if( hContour->val !=GIP_NULL )		GIP_free( hContour->val );
	hContour->val=GIP_alloc( sizeof(double)*nz );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void GIP_CONTOUR_clear( GIP_CONTOUR *hContour )	{
	hContour->M=-1;			hContour->N=-1;
	hContour->nz=-1;
	if( hContour->ptr !=GIP_NULL )		GIP_free( hContour->ptr );
	hContour->ptr=GIP_NULL;
	if( hContour->ind !=GIP_NULL )		GIP_free( hContour->ind );
	hContour->ind=GIP_NULL;
	if( hContour->val !=GIP_NULL )		GIP_free( hContour->val );
	hContour->val=GIP_NULL;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void GIP_CONTOUR_0( int M,int N,double *u,GIP_CONTOUR *hContour,int flag )	{
	int r,c,pos,nz=0,*ptr,*ind;
	double a,*contour;

	pos = 0;
	for( c = 0; c < N; c++ )	{		//列优先存储
		if( u[pos++]==0.0 )		nz++;
		for( r = 1; r < M; r++,pos++ )	{
			a = u[pos];
			if( a == 0.0 )
				nz++;
			else if( a*u[pos-1] < 0.0 )
				nz++;
		}
	}

	GIP_CONTOUR_init( hContour,M,N,nz );
	ptr = hContour->ptr;		ind = hContour->ind;
	contour = hContour->val;
	ptr[0]=0;
	pos = 0;		nz = 0;
	for( c = 0; c < N; c++ )	{		//列优先存储
		if( u[pos++]==0.0 )		{
			ind[nz]=0;
			contour[nz] = 1.0;
			nz++;

		}
		for( r = 1; r < M; r++,pos++ )	{
			a = u[pos];
			if( a == 0.0 || a*u[pos-1] < 0.0 )	{
				ind[nz]=r;
				contour[nz] = 1.0;
				nz++;
			}
		}
		ptr[c+1]=nz;
	}
	ASSERT( nz==hContour->nz );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/21/2008	
*/
void GIP_CONTOUR_1( int M,int N,double *u,double *contour,int ldu,int flag )	{
	int r,c,pos;
	double a;

	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos = r*ldu+c;
			contour[pos] = 1.0;
		}
	}   
	for( r = 0; r < M; r++ )	{
		for( c = 0; c < N; c++ )	{
			pos = r*ldu+c;
			a = u[pos];
			if( a==0.0 )		{
				contour[pos]= 0.0;
			}else if( c>0 && a*u[pos-1]<0.0 )	{
				contour[pos]= 0.0;		contour[pos-1]= 0.0;
			}else if( r > 0 && a*u[pos-ldu]<0.0 )	{
				contour[pos]= 0.0;		contour[pos-ldu]= 0.0;
			}
		}
	}
}

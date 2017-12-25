#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include <limits.h>
#include <math.h>
#include "../PCH/GIP_def.h"
#include "gip_imag_property.h"

/*
	ISSUE-NEEDED:
		由于眼角阴影的关系，r_0偏低
		
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/3/2008	
*/
void GIP_centre( int M,int N,int ldu,GIP_FLOATING *u,double *temp,double *c_r,double *c_c,int flag )	{
	int r,c,i,pos;
	double r_0,c_0,r_min,c_min,*r_prof,*c_prof,a;
	
	r_prof=temp,		c_prof=r_prof+M;
	*c_r=-1,			*c_c=-1;

	for( i = 0; i < M+N; i++ )	r_prof[i]=0.0;

	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
		a =u[pos];
		c_prof[c] += a;			r_prof[r]+=a;
	}
	}
	ASSERT( M>=3 && N>=3 );
	c_min=INT_MAX;		c_0=-1;
	for( i = 1; i < N-1; i++ )	{		//采用3-平均
		a = c_prof[i-1]+c_prof[i]+c_prof[i+1];
		if( a<c_min )	
		{	c_0=i;	c_min=a;		}
	}
	r_min=INT_MAX;		r_0=-1;
	for( i = 1; i < M-1; i++ )	{
		a = r_prof[i-1]+r_prof[i]+r_prof[i+1];
		if( a<r_min )	
		{	r_0=i;	r_min=a;		}
	}

	*c_r=r_0,		*c_c=c_0;
}


/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/25/2008	
*/
void GIP_mean_deviation_abs( int M,int N,int ldu,GIP_FLOATING *u,double *mean,double *devia )	{
	int r,c,pos;
	double avg;

	avg = 0.0;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
		avg +=fabs(u[pos]);
	}
	}
	*mean = avg/M/N;

	avg = 0.0;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
		avg +=fabs( fabs(u[pos])-*mean);
	}
	}
	*devia = avg/M/N;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/25/2008	
*/
void GIP_mean_deviation( int M,int N,int ldu,GIP_FLOATING *u,double *mean,double *devia )	{
	int r,c,pos;
	double avg;

	avg = 0.0;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
		avg +=fabs(u[pos]);
	}
	}
	*mean = avg/M/N;

	avg = 0.0;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
		avg +=fabs( fabs(u[pos])-*mean);
	}
	}
	*devia = avg/M/N;
}

/*
	Integral Images for Computing Local Means and Variances	
	[ref]: Efficient Implementation of Local Adaptive Thresholding Techniques Using Integral Images

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/20/2009	
*/
void GIP_local_mean_deviation( int M,int N,int ldu,int window,GIP_FLOATING *u,double *mean,double *devia )	{
}
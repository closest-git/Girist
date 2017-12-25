#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include <limits.h>
#include "../PCH/GIP_def.h"
#include "GIP_fdm.h"
#include "GIP_util.h"

/*
	rank_1 finite difference
	boundary condition:		extend u by reflection outside the domain
*/
double FD_R1_0( int isBndry,GIP_FD_METHOD fdm,double *u,int pos,int h )	{
	int scheme=GIP_FD_SCHEME(fdm);
	double a=0.0;
	switch(scheme)	{
		case GIP_SCHEME_FORWARD:
			a = (isBndry==1) ? u[pos-h]-u[pos] : u[pos+h]-u[pos];
			break;
		case GIP_SCHEME_BACKWARD:
			a = (isBndry==1) ? u[pos]-u[pos+h] : u[pos]-u[pos-h];
			break;
		default:
			ASSERT(0);
			break;
	}
	
	return a;
}



/*
	rank_2 finite difference
	boundary condition:		extend u by reflection outside the domain

	v0.1	cys
		8/28/2008

double FD_R2_0( int tag,GIP_FD_METHOD fdm,double *u,int pos,int h )	{
	double u1,u2,u3,u4;
	u1=BIT_TEST( tag,GIP_BD_LEFT ) ? u[pos+1] : u[pos-1];
	u2=BIT_TEST( tag,GIP_BD_RIGHT ) ? u[pos-1] : u[pos+1];
	u3=BIT_TEST( tag,GIP_BD_DOWN ) ? u[pos+h] : u[pos-h];
	u4=BIT_TEST( tag,GIP_BD_UP ) ? u[pos-h] : u[pos+h];

	return u1+u2-u[pos]*2;
}
*/

/*
	laplacian operator 
	boundary condition:		extend u by reflection outside the domain

	v0.1	cys
		9/4/2008
*/
double FD_LAPLAS_0( int tag,GIP_FD_METHOD fdm,double *u,int pos,int h )	{
	double u1,u2,u3,u4;
	u1=BIT_TEST( tag,GIP_BD_LEFT ) ? u[pos+1] : u[pos-1];
	u2=BIT_TEST( tag,GIP_BD_RIGHT ) ? u[pos-1] : u[pos+1];
	u3=BIT_TEST( tag,GIP_BD_DOWN ) ? u[pos+h] : u[pos-h];
	u4=BIT_TEST( tag,GIP_BD_UP ) ? u[pos-h] : u[pos+h];

	return u1+u2+u3+u4-u[pos]*4;
}

/*
	| 0 1 2 |
	| 3 * 5 |
	| 6 7 8 |

	注意：
		1	列优先排序，并包含(r,c)本身
		2	w<0表示插值

	v0.1	cys
		9/8/2008
*/
#define FDM_STENCIL_1 9
#define _STENCIL_STEP_1 3
static const int _STENCIL_OPPO_R[FDM_STENCIL_1]={
	6,7,8,INT_MIN,INT_MIN,INT_MIN,0,1,2
};
static const int _STENCIL_OPPO_C[FDM_STENCIL_1]={
	2,INT_MIN,0,5,INT_MIN,3,8,INT_MIN,6
};
static const int _STENCIL_OPPO_B[FDM_STENCIL_1]={
	8,INT_MIN,6,INT_MIN,INT_MIN,INT_MIN,2,INT_MIN,0
};
int FDM_stencil_1( int M,int N,int r,int c,int *nb_r,int *nb_c,GIP_FLOATING *nb_w )	{
	int oppo,i;

	for( i = 0; i < FDM_STENCIL_1; i ++ )	{
		nb_r[i]=r+(i/_STENCIL_STEP_1-1);
		nb_c[i]=c+(i%_STENCIL_STEP_1-1);
		nb_w[i]=1;
	}
	for( i = 0; i < FDM_STENCIL_1; i ++ )	{		//boundary contidion
		if( !GIP_RC_VALID( nb_r[i],0,M-1) && !GIP_RC_VALID( nb_c[i],0,N-1) )	{
			oppo = _STENCIL_OPPO_B[i];
		}else if( !GIP_RC_VALID( nb_r[i],0,M-1) )	{
			oppo = _STENCIL_OPPO_R[i];
		}
		else if( !GIP_RC_VALID( nb_c[i],0,N-1) )	{
			oppo = _STENCIL_OPPO_C[i];
		}
		else
			continue;
		nb_w[i] --;
		ASSERT( GIP_RC_VALID( nb_r[oppo],0,M-1) && GIP_RC_VALID( nb_c[oppo],0,N-1) );
		nb_w[oppo] ++;
	}

	return 0;
}

/*
	生成基于median scheme的矩阵,CRS格式

	*[注]	9/7/2008 我的34周岁生日，听"英雄的黎明"，"故乡的原风景"

	v0.1	cys
		9/7/2008
*/
int FDM_init_equations_m_( int M,int N,int ldu,GIP_FD_METHOD fdm,int flag,GIP_FLOATING alpha,GIP_FLOATING diag,GIP_FD_EQUATION *equation )	{
	int dim,nz_limit,nnz,r,c,pos,rank,set;
	int *ptr,*ind,ret=-1,nb_r[FDM_STENCIL_1],nb_c[FDM_STENCIL_1],i;
	GIP_FLOATING *val,s,nb_w[FDM_STENCIL_1];

//	ASSERT( BIT_TEST(alg,CGM_METHOD_MEDIAN)) );

	rank=1;			
	
	dim=M*N;
	nz_limit=dim*9;
	ptr=(int*)GIP_alloc( sizeof(int)*(dim+1) );
	ind=(int*)GIP_alloc( sizeof(int)*(nz_limit) );
	val=(double*)GIP_alloc( sizeof(GIP_FLOATING)*(nz_limit) );

	ptr[0]=0;			nnz=0;
	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
		switch( rank )	{
		case 1:		
			set=FDM_STENCIL_1-1;
			s=alpha*2.0/rank/rank;
			FDM_stencil_1( M,N,r,c,nb_r,nb_c,nb_w );
			for( i = 0; i < FDM_STENCIL_1; i++ )	{
				if( nb_w[i]==0.0 )
					continue;
				ind[nnz]= GIP_RC2POS(nb_r[i],nb_c[i],ldu );	
				if( ind[nnz]==pos )			//diagnol element
					val[nnz]= -s+diag;
				else
					val[nnz]= s/set*nb_w[i];

				nnz++;
			}
/*			if( r>0 )	{
				val[nnz]=(r==M-1) ? s/set*2 : s/set;
				ind[nnz]=GIP_RC2POS( r-1,c,ldu );	
				nnz++;
			}
			if( c>0 )	{		
				val[nnz]=(c==N-1) ? s/set*2 : s/set;
				ind[nnz]=GIP_RC2POS( r,c-1,ldu );	
				nnz++;
			}
			val[nnz]=-s+diag;		ind[nnz++]=pos;	
			if( c<N-1 )			{
				val[nnz]=(c==0) ? s/set*2 : s/set;
				ind[nnz]=GIP_RC2POS( r,c+1,ldu );
				nnz++;
			}
			if( r<M-1 )			{
				val[nnz]=(r==0) ? s/set*2 : s/set;
				ind[nnz]=GIP_RC2POS( r+1,c,ldu );	
				nnz++;								
			}*/
			break;
		default:
			break;
		}
		ptr[pos+1]=nnz;
	} 
	}

	ASSERT( nnz<=dim*9 );
	equation->dim=dim;			equation->nz=nnz;
	equation->ptr=ptr;			equation->ind=ind;				equation->val=val;
	ret=0;

	return ret;
}
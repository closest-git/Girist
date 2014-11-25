/*
	IMAGE decomposition:	following the ideas of Yves Meyer, 

	Fixed point iteration with implicit format
	Boundary condition
		u:	reflection outside the domain
		g:	Dirichlet boundary condition
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
#include "meyer_idea.h"
#include "GIP_fmm.h"
#include "../../lib/gss_6_demo.h"

static int meyer_step=0;

/*
	u定义为[0:M+1,0:N+1],其中u0为矩形{r=1,r=M,c=1,c=N}

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/2/2008	
*/
void MEYER_Init_u0( MEYER_SPLIT *hMEYER,GIP_IMAGE *hMask,int type )	{
	int M=hMEYER->M,N=hMEYER->N,r,c,pos,isOut,ldu=hMEYER->ldu;
	int r0=10,r1=M-10,c0=20,c1=N-20;		//
	double *u=hMEYER->u,ar,ac,dist;

	if( hMask==0x0 )	{
		for( r = 0; r < M; r++ )	{
			for( c = 0; c < N; c++ )	{
				pos = GIP_RC2POS( r,c,ldu );
				isOut = 0;
				ar = MIN( fabs(r-r0),fabs(r-r1) );
				ac = MIN( fabs(c-c0),fabs(c-c1) );
				if( r<r0 || r>r1 )
				{	isOut=1;	dist = ar;	}
				else if( c<c0 || c>c1 )
				{	isOut=1;	dist = ac;	}
				else
					dist = MIN(ar,ac);
				u[pos] = (isOut==1) ? dist : -dist;
				if( u[pos]==0.0 )			{	
					BIT_SET( hMEYER->tag[pos],GIP_INTERFACE );	
				}
			}
		}
	}else	{		//based on mask imag
		int height=hMask->height,width=hMask->width,step=hMask->widthStep/sizeof(uchar);
		uchar *i_data = (uchar *)hMask->imageData,a;
		ASSERT( height==M && width==N && hMask->depth==IPL_DEPTH_8U && hMask->nChannels==1 );
		for( r = 0; r < M+2; r++ )	{		//暂定点全部在边界之外
			for( c = 0; c < N+2; c++ )	{
				pos = GIP_RC2POS( r,c,ldu );
				u[pos]=1.0;				
			}
		}
		for( r = 0; r < height; r++ )	{
			for( c = 0; c < width; c++ )	{
				a = i_data[r*step+c];
				pos = GIP_RC2POS( r+1,c+1,ldu );	
				if( a == 0 )	{
					u[pos]=0.0;		
					BIT_SET( hMEYER->tag[pos],GIP_INTERFACE );
				}else if( a != UCHAR_MAX )	{		//边界内部的点
					u[pos]=-1.0;
				}
			}
		}
			
	}

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/14/2008	
*/
void MEYER_output( MEYER_SPLIT *hMEYER,int flag )	{
	int M=hMEYER->M,N=hMEYER->N,MN=M*N,i,j,ldu=hMEYER->ldu,*tag=hMEYER->tag,ret=GIP_OK,t;
	int dim=hMEYER->A_dim,*ptr=hMEYER->A_ptr,*ind=hMEYER->A_ind,nnz,r,c,pos,*crs_ptr,*crs_ind,*temp;
	double *rhs=hMEYER->rhs,*crs_val,*u=hMEYER->u,*val=hMEYER->A_val;
	FILE *fp=NULL;
	char sPath[80];
	
	nnz = ptr[dim];
	crs_ptr=GIP_alloc( sizeof(int)*(dim+1+nnz) );		crs_ind=crs_ptr+dim+1;
	crs_val=GIP_alloc( sizeof(double)*nnz );
	temp=GIP_alloc( sizeof(int)*(dim+1) );
	memset( temp,0x0,sizeof(int)*(dim+1) );
	for( i = 0; i < nnz; i++ )	{
		r = ind[i];
		temp[r]++;
	}
	crs_ptr[dim]=nnz;
	for( i = dim-1;	i >= 0; i-- )	{
		crs_ptr[i] = crs_ptr[i+1]-temp[i];
		temp[i] = crs_ptr[i];
	}
	ASSERT( crs_ptr[0] == 0 );

	for( i = 0;	i < dim; i++ )	{
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			r = ind[j];
			crs_ind[temp[r]] = i;		crs_val[temp[r]] = val[j];
			temp[r]++;
		}
	}
	for( i = 0;	i < dim; i++ )	ASSERT( temp[i]==crs_ptr[i+1] );

	_stprintf( sPath,"H:\\GIPack\\data\\MEYER_%d.mtx",dim );
	G_PRINTF( "<<%s>> outputing...\r",sPath);
	fp = _tfopen( sPath,"w" );
	fprintf( fp,"no\tu\trhs\n" );
	for( i = 0;	i < dim; i++ )	{
		fprintf( fp,"%d\t%g\t%g\n",i+1,u[i],rhs[i] );
	}
	fprintf( fp,"************************\r\n",dim,nnz );
	fprintf( fp,"%d %d\n",dim,nnz );
	for( i = 0;	i < dim; i++ )	{
		for( j = crs_ptr[i]; j <crs_ptr[i+1]; j++ )	{
			c = crs_ind[j];
			fprintf( fp,"%d\t%d\t%lf\n",i+1,c+1,crs_val[j] );
		}
	}
	fclose( fp );
	G_PRINTF( "<<%s>> outputing finish.\n",sPath);

	GIP_free(crs_ptr);			
	GIP_free(crs_val);
}

/*
	[REF] Modeling Textures with Total Variation Minimization and Oscillating Patterns in Image Processing

	注意：
		1、hMEYER->ldu=N

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/5/2008	
*/
void MEYER_init_Lp( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,MEYER_SPLIT *hMEYER )	{
	int M,N,MN,r,c,pos,A_nz,i,setting[32],ret=GIP_OK,ldu,BLK_shift;
	GIP_FLOATING *val=GIP_NULL;
	double *g,a,s,*f,fx,fy;

	meyer_step=0;

	memset( hMEYER,0x0,sizeof(MEYER_SPLIT) );
	GIP_import_IplImage( hIMG,&M,&N,&val,0x0 );

	hMEYER->alg = MEYER_ALG_LP | MEYER_SS_DIRECT;			//MEYER_SS_ITER
	hMEYER->lenda=0.1;
	hMEYER->mui =0.1;			hMEYER->beta = 1.0e-12;
	hMEYER->M=M;				hMEYER->N=N;
	MN=M*N;
	hMEYER->ldu=N;			hMEYER->shift=0;
	hMEYER->f=(double*)GIP_alloc( sizeof(double)*(MN) );
	hMEYER->u=(double*)GIP_alloc( sizeof(double)*MN*3 );
	hMEYER->g1=hMEYER->u+MN;		hMEYER->g2=hMEYER->g1+MN;
	hMEYER->tag=(int*)GIP_alloc( sizeof(int)*(M+2)*(N+2) );
	memset( hMEYER->tag,0x0,sizeof(int)*(M+2)*(N+2) );
	ldu=hMEYER->ldu;			
	s = 1.0/(2*hMEYER->lenda);
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			hMEYER->f[pos] = val[pos];
			hMEYER->u[pos] = val[pos];
			if( c==0 )		BIT_SET( hMEYER->tag[pos],GIP_BD_LEFT );
			if( c==N-1 )	BIT_SET( hMEYER->tag[pos],GIP_BD_RIGHT );
			if( r==0 )		BIT_SET( hMEYER->tag[pos],GIP_BD_DOWN );
			if( r==M-1 )	BIT_SET( hMEYER->tag[pos],GIP_BD_UP );
			if( c==0 || c==N-1 || r==0 || r==M-1 )	{
				hMEYER->g1[pos] = 0.0;		hMEYER->g2[pos] = 0.0;		
			}else	{
				fx=(val[GIP_RC2POS(r,c+1,ldu)]-val[GIP_RC2POS(r,c-1,ldu)])/2.0;
				fy=(val[GIP_RC2POS(r+1,c,ldu)]-val[GIP_RC2POS(r-1,c,ldu)])/2.0;
				hMEYER->g1[pos] = -s*fx/sqrt(fx*fx+fy*fy+hMEYER->beta);				//0.0;	
				hMEYER->g2[pos] = -s*fy/sqrt(fx*fx+fy*fy+hMEYER->beta);				//0.0;	
			}			
		}
	}

	hMEYER->A_dim=MN*3;			
	hMEYER->m_type=0;		//unsym
	A_nz=hMEYER->m_type==11 || hMEYER->m_type==12 ? hMEYER->A_dim*6 : hMEYER->A_dim*11;			//取上限
	hMEYER->rhs=(int*)GIP_alloc( sizeof(double)*(hMEYER->A_dim) );
	hMEYER->A_ptr=(int*)GIP_alloc( sizeof(int)*(hMEYER->A_dim+1) );
	hMEYER->A_ind=(int*)GIP_alloc( sizeof(int)*(A_nz) );
	hMEYER->A_val=(double*)GIP_alloc( sizeof(double)*(A_nz) );
	hMEYER->A_ptr[0]=0;			A_nz=0;
	BLK_shift=M*N;
	ASSERT( hMEYER->m_type!=11 && hMEYER->m_type!=12 );
//	U Block		GRID2MAT
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			if( r>0 )	hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu );
			if( c>0 )	hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu );
			hMEYER->A_ind[A_nz++]=pos;
			if( c<N-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu );
			if( r<M-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu );
			if( c>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu )+BLK_shift;
			if( c<N-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu )+BLK_shift;
			if( r>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu )+BLK_shift*2;
			if( r<M-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu )+BLK_shift*2;
			hMEYER->A_ptr[pos+1]=A_nz;
		}
	}
//	g1 Block		GRID2MAT
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu )+BLK_shift;
			if( c>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu );		//u-g1
			if( c<N-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu );
			if( c>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu )+BLK_shift;
			hMEYER->A_ind[A_nz++]=pos;
			if( c<N-1 )				hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu )+BLK_shift;
			if( r>0 && c>0 && r<M-1&&c<N-1 )	{
			if( r>0 && c>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c-1,ldu )+BLK_shift*2;
			if( r>0 && c<N-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c+1,ldu )+BLK_shift*2;
			if( r<M-1&&c>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c-1,ldu )+BLK_shift*2;
			if( r<M-1&&c<N-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c+1,ldu )+BLK_shift*2;
			}
/*			if( r>0 && c>0 )hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c-1,ldu )+BLK_shift*2;
			if( r>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu )+BLK_shift*2;
			if( c>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu )+BLK_shift*2;
							hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c,ldu )+BLK_shift*2;
			if( c<N-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu )+BLK_shift*2;
			if( r<M-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu )+BLK_shift*2;
			if( r<M-1&&c<N-1 )	hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c+1,ldu )+BLK_shift*2;*/
			hMEYER->A_ptr[pos+1]=A_nz;
		}
	}
//	g2 Block		GRID2MAT
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu )+BLK_shift*2;
			if( r>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu );		//u-g1
			if( r<M-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu );
			if( r>0 && c>0 && r<M-1&&c<N-1 )	{
			if( r>0 && c>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c-1,ldu )+BLK_shift;
			if( r>0 && c<N-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c+1,ldu )+BLK_shift;
			if( r<M-1&&c>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c-1,ldu )+BLK_shift;
			if( r<M-1&&c<N-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c+1,ldu )+BLK_shift;
			}
/*			if( r>0 && c>0 )hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c-1,ldu )+BLK_shift;
			if( r>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu )+BLK_shift;
			if( c>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu )+BLK_shift;
							hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c,ldu )+BLK_shift;
			if( c<N-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu )+BLK_shift;
			if( r<M-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu )+BLK_shift;
			if( r<M-1&&c<N-1 )	hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c+1,ldu )+BLK_shift;*/
			if( r>0 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu )+BLK_shift*2;
			hMEYER->A_ind[A_nz++]=pos;
			if( r<M-1 )		hMEYER->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu )+BLK_shift*2;
			hMEYER->A_ptr[pos+1]=A_nz;
		}
	}
	hMEYER->A_nz = A_nz;			ASSERT( A_nz<=hMEYER->A_dim*9 );

	f=hMEYER->f;
	memset( hMEYER->rhs,0x0,sizeof(double)*hMEYER->A_dim );
	memcpy( hMEYER->rhs,f,sizeof(double)*MN );
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
		if( c>0 && c<N-1 )		
//		{	hMEYER->rhs[pos+MN]=(f[pos]-f[pos-1]);		}
		{	hMEYER->rhs[pos+MN]=-(f[pos+1]-f[pos-1])/2;		}
		if( r>0 && r<M-1 )		
//		{	hMEYER->rhs[pos+MN*2]=(f[pos]-f[pos-ldu]);		}
		{	hMEYER->rhs[pos+MN*2]=-(f[pos+ldu]-f[pos-ldu])/2;		}
	}
	}

	hMEYER->hSS=GIP_NULL;
	if( BIT_TEST( hMEYER->alg,MEYER_SS_DIRECT) )	{
		for( i = 0; i < 32; i++ )	setting[i]=0.0;
		if( (ret=GSS_init_id(hMEYER->A_dim,hMEYER->A_dim,hMEYER->A_ptr,hMEYER->A_ind,GIP_NULL,hMEYER->m_type,setting) ) != GRUS_OK )	
		{	GIP_ERROR( "TVM_init",ret,"\tERROR at init GSS solver. ERROR CODE:%d\r\n" );		}
		hMEYER->hSS = GSS_symbol_id( hMEYER->A_dim,hMEYER->A_dim,hMEYER->A_ptr,hMEYER->A_ind,hMEYER->A_val );	
		if( hMEYER->hSS == NULL )	
		{	GIP_ERROR( "TVM_init",-2,"\tERROR at SYMBOLIC ANALYSIS.\r\n" );		}
	}

/*
	MEYER_Init_u0( hMEYER,hIMG_mask,0x0 );
	GIP_FMM_smooth_u_( hMEYER->M,hMEYER->N,hMEYER->u,hMEYER->ldu,hMEYER->tag,0x0 );
*/
GIP_exit:	
	GIP_free( val );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void MEYER_clear( MEYER_SPLIT *hMEYER )	{
	if( hMEYER == GIP_NULL )	return;

	hMEYER->M=-1;		hMEYER->N=-1;
	if( hMEYER->A_ptr != GIP_NULL )	
	{	GIP_free( hMEYER->A_ptr );	hMEYER->A_ptr=GIP_NULL; }
	if( hMEYER->A_ind != GIP_NULL )	
	{	GIP_free( hMEYER->A_ind );	hMEYER->A_ind=GIP_NULL; }
	if( hMEYER->A_val != GIP_NULL )	
	{	GIP_free( hMEYER->A_val );	hMEYER->A_val=GIP_NULL; }
	if( hMEYER->u != GIP_NULL )	
	{	GIP_free( hMEYER->u );	hMEYER->u=GIP_NULL; }
	if( hMEYER->f != GIP_NULL )	
	{	GIP_free( hMEYER->f );	hMEYER->f=GIP_NULL; }
	if( hMEYER->tag != GIP_NULL )	
	{	GIP_free( hMEYER->tag );	hMEYER->tag=GIP_NULL; }
}

/*

	v0.1	cys
		5/30/2008	
*/
void MEYER_Du( MEYER_SPLIT *hMEYER,double *Va,double *Du )	{
	int M=hMEYER->M,N=hMEYER->N,ldu=hMEYER->ldu,r,c,pos;
	double ux,uy,g,beta=hMEYER->beta;

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			ux = c==N-1 ? 0.0 : Va[pos+1]-Va[pos];			//delta+
			uy = r==M-1 ? 0.0 : Va[pos+ldu]-Va[pos];		//delta+
			g=sqrt(ux*ux+uy*uy+beta );
			Du[pos]=1.0/g;
		}
	}
}

/*
	The matrix U:U is symmetric positive definite and block tridiagonal with positive diagonal entries and 
	negative off-diagonal entries and has rows which sum to 1.

	v0.1	cys
		8/10/2008	
*/
void MEYER_verify_UU( MEYER_SPLIT *hMEYER,int flag )	{
	int dim=hMEYER->A_dim,*ptr=hMEYER->A_ptr,*ind=hMEYER->A_ind,isCheckRes=1,M=hMEYER->M,N=hMEYER->N,MN;
	int isFind,i,j,k,r;
	double *val=hMEYER->A_val,a,sum;

	MN=M*N;
	for( i = 0; i < MN; i++ )	{
		sum=0.0;
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{			
			r = ind[j];
			if( r>=MN )	{
				continue;
			}
			if( r==i )	ASSERT( val[j]>=1.0 );
			if( r > i )	{
				isFind = 0;
				for( k = ptr[r]; k < ptr[r+1];	k++ )	{
					if( ind[k]==i && val[k]==val[j] )
						isFind=1;
				}
				ASSERT( isFind==1 );
			}
			sum += val[j];
		}
		ASSERT( fabs(sum-1.0)<FLT_EPSILON );
	}
}

/*
	[REF] Modeling Textures with Total Variation Minimization and Oscillating Patterns in Image Processing

	v0.1	cys
		8/7/2008	
*/
void MEYER_update_La( MEYER_SPLIT *hMEYER,double *Va,double *rhs )	{
	int dim=hMEYER->A_dim,*ptr=hMEYER->A_ptr,*ind=hMEYER->A_ind,ldu=hMEYER->ldu,M=hMEYER->M,N=hMEYER->N;
	int i,r,c,pos,r_D,c_D,A_nz,pos_D,alg=hMEYER->alg,BLK_shift,nz_D;
	double *val=hMEYER->A_val,*u,*g1=hMEYER->g1,*g2=hMEYER->g2,*f=hMEYER->f,mui=hMEYER->mui,s,
		a,a_max,a_min,beta=hMEYER->beta,x,y;
	
	BLK_shift=M*N;
	u = hMEYER->u;
	A_nz=0;
	s = 1.0/(hMEYER->lenda*2);
//	U column		
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
//		U:U
			pos_D=GIP_RC2POS( r,c,ldu );
			a=0.0;
			if( r>0 )	{	//=(i,j-1)
				pos = GIP_RC2POS( r-1,c,ldu );
				x=c>0&&c<N-1?(u[pos+1]-u[pos-1])/2.0:0;		y=u[pos_D]-u[pos];
				a += s/sqrt(x*x+y*y+beta);		val[A_nz++]=-s/sqrt(x*x+y*y+beta);
			}
			if( c>0 )	{	//=(i-1,j)
				pos = GIP_RC2POS( r,c-1,ldu );
				x=u[pos_D]-u[pos];				y=r>0&&r<M-1?(u[pos+ldu]-u[pos-ldu])/2.0:0;		
				a += s/sqrt(x*x+y*y+beta);		val[A_nz++]=-s/sqrt(x*x+y*y+beta);
			}
			nz_D=A_nz++;
			if( c<N-1 )		{	//=(i+1,j)
				pos = GIP_RC2POS( r,c+1,ldu );
				x=u[pos_D]-u[pos];				y=r>0&&r<M-1?(u[pos_D+ldu]-u[pos_D-ldu])/2.0:0;
				a += s/sqrt(x*x+y*y+beta);		val[A_nz++]=-s/sqrt(x*x+y*y+beta);
			}
			if( r<M-1 )		{	//=(i,j+1)
				pos = GIP_RC2POS( r+1,c,ldu );
				x=c>0&&c<N-1?(u[pos_D+1]-u[pos_D-1])/2.0:0;	y=u[pos_D]-u[pos];
				a += s/sqrt(x*x+y*y+beta);		val[A_nz++]=-s/sqrt(x*x+y*y+beta);
			}
			val[nz_D]=1.0+a;
//		G1:U
			if( c>0 )		
			{	pos = GIP_RC2POS( r,c-1,ldu );		val[A_nz++]=1.0/2;		}
			if( c<N-1 )		
			{	pos = GIP_RC2POS( r,c+1,ldu );		val[A_nz++]=-1.0/2;		}
//		G2:U
			if( r>0 )		
			{	pos = GIP_RC2POS( r-1,c,ldu );		val[A_nz++]=1.0/2;		}
			if( r<M-1 )		
			{	pos = GIP_RC2POS( r+1,c,ldu );		val[A_nz++]=-1.0/2;		}
			ASSERT( A_nz==hMEYER->A_ptr[GIP_RC2POS( r,c,ldu )+1] );
		}
	}
//	G1 column		
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
//		U:G1
			pos_D = GIP_RC2POS( r,c,ldu );
			if( c>0 )		
			{	pos = GIP_RC2POS( r,c-1,ldu );		val[A_nz++]=-1.0/2;		}
			else
			{	a=r>0&&r<M-1?(g2[pos_D+ldu]-g2[pos_D-ldu])/2.0:0.0;		rhs[pos_D]-=f[pos_D]-u[pos_D]-a;	}
			if( c<N-1 )		
			{	pos = GIP_RC2POS( r,c+1,ldu );		val[A_nz++]=1.0/2;		}
			else
			{	a=r>0&&r<M-1?(g2[pos_D+ldu]-g2[pos_D-ldu])/2.0:0.0;		rhs[pos_D]-=f[pos_D]-u[pos_D]-a;	}
//		G1:G1
			if( c>0 )		
			{	pos = GIP_RC2POS( r,c-1,ldu );		val[A_nz++]=-1.0;		}
			hMEYER->A_val[A_nz++]=mui*s/sqrt(g1[pos]*g1[pos]+g2[pos]*g2[pos]+beta)+2.0;
			if( c<N-1 )		
			{	pos = GIP_RC2POS( r,c+1,ldu );		val[A_nz++]=-1.0;		}
//		G2:G1
			if( r>0 && c>0 && r<M-1&&c<N-1 )	{
			if( r>0 && c>0 )
			{	pos = GIP_RC2POS( r-1,c-1,ldu );	val[A_nz++]=-1.0/4.0;		}
			if( r>0 && c<N-1 )		
			{	pos = GIP_RC2POS( r-1,c+1,ldu );	val[A_nz++]=1.0/4.0;		}
			if( r<M-1 && c>0 )
			{	pos = GIP_RC2POS( r+1,c-1,ldu );	val[A_nz++]=1.0/4.0;		}
			if( r<M-1&&c<N-1 )	
			{	pos = GIP_RC2POS( r+1,c+1,ldu );	val[A_nz++]=-1.0/4.0;		}
			}
			ASSERT( A_nz==hMEYER->A_ptr[GIP_RC2POS( r,c,ldu )+BLK_shift+1] );
		}
	}
//	G2 column		
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu )+BLK_shift*2;
//		U:G2
			if( r>0 )		
			{	pos = GIP_RC2POS( r-1,c,ldu );		val[A_nz++]=-1.0/2;		}
			if( r<M-1 )		
			{	pos = GIP_RC2POS( r+1,c,ldu );		val[A_nz++]=1.0/2;		}
//		G1:G2
			if( r>0 && c>0 && r<M-1&&c<N-1 )	{
			if( r>0 && c>0 )
			{	pos = GIP_RC2POS( r-1,c-1,ldu );	val[A_nz++]=-1.0/4.0;		}
			if( r>0 && c<N-1 )		
			{	pos = GIP_RC2POS( r-1,c+1,ldu );	val[A_nz++]=1.0/4.0;		}
			if( r<M-1 && c>0 )
			{	pos = GIP_RC2POS( r+1,c-1,ldu );	val[A_nz++]=1.0/4.0;		}
			if( r<M-1&&c<N-1 )	
			{	pos = GIP_RC2POS( r+1,c+1,ldu );	val[A_nz++]=-1.0/4.0;		}
			}
//		G2:G2
			if( r>0 )		
			{	pos = GIP_RC2POS( r-1,c,ldu );		val[A_nz++]=-1.0;		}
			pos = GIP_RC2POS( r,c,ldu );
			hMEYER->A_val[A_nz++]=mui*s/sqrt(g1[pos]*g1[pos]+g2[pos]*g2[pos]+beta)+2.0;
			if( r<M-1 )		
			{	pos = GIP_RC2POS( r+1,c,ldu );		val[A_nz++]=-1.0;		}

			ASSERT( A_nz==hMEYER->A_ptr[GIP_RC2POS( r,c,ldu )+BLK_shift*2+1] );
		}
	}

	MEYER_verify_UU( hMEYER,0x0 );
	a_max=DBL_MIN;		a_min=DBL_MAX;
	for( i = 0; i < A_nz; i++ )		{
		a = fabs(val[i]); 
		a_max=MAX(a_max,a);		a_min=MIN(a_min,a);	
	}

	if( dim<100 )
		MEYER_output( hMEYER,0x0 );
}

/*
	rhs=u+dt*delta*mu(lenda2-lenda1)		[17]

	v0.1	cys
		5/30/2008	

void MEYER_update_rhs( MEYER_SPLIT *hMEYER,double *rhs,double *temp )	{
	int M=hMEYER->M,N=hMEYER->N,MN=M*N,r,c,pos,ldu=hMEYER->ldu,alg=hMEYER->alg;
	double *u=hMEYER->u,*z=hMEYER->z,lenda1,lenda2,s,a,b,d,c1,c2,eps,delta;
	double eps_3=1.0e-2,eps_4=eps_3/sqrt(2.0),silen,alpha,s_max=0.0;

	lenda1=hMEYER->lenda,		lenda2=hMEYER->lenda;
	c1=hMEYER->c1,				c2=hMEYER->c2;
	eps=hMEYER->eps;

	s = hMEYER->dt;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{	
		d=0.0;
		pos = GIP_RC2POS( r,c,ldu );
		a = (z[pos]-c1);			b = (z[pos]-c2);
//		if( u[pos]>=0.0 )	{
			d -= lenda1*a*a;
//		}else	{
			d += lenda2*b*b;
//		}
		s_max = MAX( s_max,a*a-b*b );
		a = u[pos];
		delta = eps/( (eps*eps+a*a) * GIP_PI );
		d *= delta*s;
		rhs[pos] = a+d*s;
	}
	}


	return;
}
*/

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/8/2008	
*/
void MEYER_step( MEYER_SPLIT *hMEYER,double dT,double *temp )	{
	int M=hMEYER->M,N=hMEYER->N,MN=M*N,i,j,r,c,g_pos,pos,ldu=hMEYER->ldu,*tag=hMEYER->tag,ret=GIP_OK;
	int dim=hMEYER->A_dim,*ptr=hMEYER->A_ptr,*ind=hMEYER->A_ind,isCheckRes=1,MN_size;
	double *rhs=hMEYER->rhs,*u,u_max,u_min,rhs_norm=0.0,*val=hMEYER->A_val,res,a,*t_rhs;

	MN_size=sizeof(double)*MN;
	u=hMEYER->u;
	t_rhs = temp+dim;
	memcpy( t_rhs,rhs,sizeof(double)*dim );
	MEYER_update_La( hMEYER,u,t_rhs );
//use GSS to solver A*u{n+1}=f(u{n})
//	if( meyer_step%5==0 )	{
		ret = GSS_numeric_id( hMEYER->A_dim,hMEYER->A_dim,hMEYER->A_ptr,hMEYER->A_ind,hMEYER->A_val,hMEYER->hSS );
		if( ret != GRUS_OK )	{
			hMEYER->hSS=NULL;		//必须设置为NULL,GSS已自动释放内存
			GIP_ERROR( "TVM_step",ret,"\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n" );
		}
//	}

	if( isCheckRes==1 )		{
		memcpy( temp,t_rhs,sizeof(double)*dim );
		for( i = 0; i < dim; i++ )	
		{	rhs_norm += (temp[i]*temp[i]);		}
		rhs_norm = sqrt(rhs_norm);	
	}
	GSS_solve_id( hMEYER->hSS,hMEYER->A_dim,hMEYER->A_dim,hMEYER->A_ptr,hMEYER->A_ind,hMEYER->A_val,t_rhs );
	memcpy( hMEYER->u,t_rhs,sizeof(double)*dim );
	if( isCheckRes==1 )	{
		res = 0.0;
		u_max=DBL_MIN,		u_min=DBL_MAX;
		for( i = 0; i < dim; i++ )	{
			a = fabs(u[i]);
			u_max=MAX( u_max,a );			u_min=MIN( u_min,a );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				r = ind[j];
				temp[r] -= val[j]*u[i];
			}		
		}
		for( i = 0; i < dim; i++ )	
		{	res += (temp[i]*temp[i]);		}
		res = sqrt(res);	
		G_PRINTF( "\t\t|Ax-b|=%g,|b|=%g,u_max=%g,u_min=%g\r\n",res,rhs_norm,u_max,u_min );	
	}

GIP_exit:
	;
}

/*
	[REF] Modeling Textures with Total Variation Minimization and Oscillating Patterns in Image Processing P.9

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/8/2008	
*/
void MEYER_step_iter( MEYER_SPLIT *hMEYER,double *temp )	{
	int M=hMEYER->M,N=hMEYER->N,MN=M*N,i,j,ldu=hMEYER->ldu,*tag=hMEYER->tag,ret=GIP_OK,t;
	int dim=hMEYER->A_dim,*ptr=hMEYER->A_ptr,*ind=hMEYER->A_ind,isCheckRes=1,MN_size,r,c,pos;
	double *rhs=hMEYER->rhs,*u,u_max,u_min,rhs_norm=0.0,*val=hMEYER->A_val,a,d,*x;

	MN_size=sizeof(double)*MN;
	u=hMEYER->u;
	x = temp+dim;
	memcpy( x,rhs,sizeof(double)*dim );
	MEYER_update_La( hMEYER,u,x );	

	u_max=DBL_MIN,		u_min=DBL_MAX;
	for( i = 0; i < dim; i++ )	{
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			r = ind[j];
			if( r==i )	{
				d = val[j];
			}else	{
				if( r==7953 )	{
					t=0;
				}
				x[r] -= val[j]*u[i];
			}
		}
	}
	for( i = 0; i < dim; i++ )	{
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			r = ind[j];
			if( r!=i )	continue;
			d = val[j];
			x[i] /= d;
			a = fabs(x[i]);
			u_max=MAX( u_max,a );			u_min=MIN( u_min,a );
		}
	}
	memcpy( hMEYER->u,x,sizeof(double)*dim );
//	memcpy( hMEYER->g1,x+MN,MN_size );
//	memcpy( hMEYER->g2,x+MN*2,MN_size );

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/8/2008	
*/
double MEYER_energy_Lp( MEYER_SPLIT *hMEYER,double *temp )	{
	int ldu=hMEYER->ldu,M=hMEYER->M,N=hMEYER->N,r,c,alg=hMEYER->alg,pos,poo=-1;
	double energy=0.0,beta=hMEYER->beta,mui=hMEYER->mui,lenda=hMEYER->lenda,a,a_max,a_min,ux,uy,g1x,g2y;
	double *u=hMEYER->u,*f=hMEYER->f,*g1=hMEYER->g1,*g2=hMEYER->g2,u_oo=0,g1_oo=0,g2_oo=0;

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			ux = c==N-1 ? 0.0 : u[pos+1]-u[pos];			//delta+
			uy = r==M-1 ? 0.0 : u[pos+ldu]-u[pos];		//delta+
			energy+=sqrt(ux*ux+uy*uy+beta);
			g1x = c==N-1 ? 0.0 : g1[pos+1]-g1[pos];
			g2y = r==M-1 ? 0.0 : g2[pos+ldu]-g2[pos];
			a = f[pos]-u[pos]-g1x-g2y;
			energy += lenda*a*a;
			energy += mui*sqrt(g1[pos]*g1[pos]+g2[pos]*g2[pos]);
			
			u_oo=MAX( u_oo,fabs(f[pos]-u[pos]) );
			if( g1_oo<fabs(g1[pos]) )	{
				g1_oo=fabs(g1[pos]);		poo=pos;
			}
//			g1_oo=MAX( g1_oo,fabs(g1[pos]) );
			g2_oo=MAX( g2_oo,fabs(g2[pos]) );
		}
	}

	G_PRINTF( "MEYER_%d: energy=%g,|f-u|=%g,|g1|=%g(%d),|g2|=%g\r\n",meyer_step,energy,u_oo,g1_oo,poo,g2_oo	);
	return energy;
}

/*	
	off=Vb-L(Va)*Va
	Normalized residual.[REF]:	ACCELERATION METHODS FOR TOTAL VARIATION-BASED IMAGE DENOISING

	v0.1	cys
		5/23/2008

double  MEYER_NormalRes( MEYER_SPLIT *hMEYER,double *temp )	{
	int i,j,M=hTVM->M,N=hTVM->N,dim=hTVM->A_dim,*ptr=hTVM->A_ptr,*ind=hTVM->A_ind;
	double off=0.0,a,*val=hTVM->A_val,d;

	MEYER_update_La( hTVM,temp );

	memcpy( temp,hMEYER->u,MN_size );
	memcpy( temp+MN,hMEYER->g1,MN_size );
	memcpy( temp+MN*2,hMEYER->g2,MN_size );
	rhs=temp+dim;
	memset( rhs,0x0,sizeof(double)*dim );
	memcpy( rhs,hMEYER->z,MN_size );

	u_max=DBL_MIN,		u_min=DBL_MAX;
	for( i = 0; i < dim; i++ )	{
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			r = ind[j];
			if( r==i )	{
				d = val[j];
			}else	{
				if( r==7953 )	{
					t=0;
				}
				rhs[r] -= val[j]*temp[i];
			}
		}
	}
	for( i = 0; i < dim; i++ )	{
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			r = ind[j];
			if( r!=i )	continue;
			d = val[j];
			rhs[i] /= d;
			a = fabs(rhs[i]);
			u_max=MAX( u_max,a );			u_min=MIN( u_min,a );
		}
	}
	for( i = 0; i < dim; i++ )	{
		a = Vb[i];
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			a -= val[j]*Va[ind[j]];
			if( ind[j]==i )
			{	d=val[j];		}
		}
		a /= d;
		off += a*a;
	}
	off = sqrt(off);

	return off;
}
*/

/*
	v=div(g1,g2)

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/11/2008	
*/
void MEYER_get_v( MEYER_SPLIT *hMEYER,double *v )	{
	int ldu=hMEYER->ldu,M=hMEYER->M,N=hMEYER->N,r,c,pos,poo=-1;
	double g1x,g2y,*g1=hMEYER->g1,*g2=hMEYER->g2;

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			g1x = c==N-1 ? 0.0 : g1[pos+1]-g1[pos];
			g2y = r==M-1 ? 0.0 : g2[pos+ldu]-g2[pos];
			v[pos] = g1x+g2y;
		}
	}	;
}

/*
	A testing code 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/5/2008	
*/
void MEYER_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName )	{
	double *temp,energy;
	MEYER_SPLIT meyer,*hMEYER=&meyer;
	int M,N,t_len,isTrace=1,nKey,t_step=1;
	char sTitle[80];

	MEYER_init_Lp( hIMG,hIMG_mask,hMEYER );
	G_PRINTF( "MEYER:\tmui=%g,lenda=%g,alg=0x%x\r\n",hMEYER->mui,hMEYER->lenda,hMEYER->alg );		

	M=meyer.M,		N=meyer.N;
	t_len=MAX((M+2)*(N+2)*6,N);
	temp=GIP_alloc( sizeof(double)*t_len );
//	tol = TVM_NormalRes( hMEYER,temp );
//	tol *= 1.0e-4;
//	GIP_show_IplImage_d( M,N,sWndName,meyer.z,meyer.ldu,0x0 );
	while( !BIT_TEST(meyer.flag,GIP_SOLVER_CONVERGE) )	{
		if( isTrace==1 && meyer_step%t_step==0 )	{
			_stprintf( sTitle,".\\trace\\u_%d.bmp",meyer_step );		
			GIP_save_IplImage_d( M,N,sTitle,(meyer.u+meyer.shift),meyer.ldu,0x0 );
			MEYER_get_v( hMEYER,temp );
			_stprintf( sTitle,".\\trace\\v_%d.bmp",meyer_step );		
			GIP_save_IplImage_d( M,N,sTitle,temp,meyer.ldu,0x0 );
//			GIP_show_IplImage_d( M,N,sWndName,temp,meyer.ldu,0x0 );
			nKey = cvWaitKey( );
			switch( nKey )	{
				case 27:
					G_PRINTF( "IMAGE decomposition iteration break out\r\n"	);
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
//			GIP_save_IplImage_d( M,N,sTitle,(temp+meyer.shift),meyer.ldu,0x0 );
		}
		if( BIT_TEST( hMEYER->alg,MEYER_SS_DIRECT ) )
			MEYER_step( &meyer,0.0,temp );
		else if( BIT_TEST( hMEYER->alg,MEYER_SS_ITER ) )
			MEYER_step_iter( &meyer,temp );
		energy=MEYER_energy_Lp( &meyer,temp ); 
		if( BIT_TEST(meyer.flag,GIP_SOLVER_CONVERGE) )
			break;
		meyer_step++;
		//if( average > 2 )	{
		//	MEYER_ReInit_u0( &meyer,0x0 );
		//	MEYER_smooth_u_( &meyer );
		//}
	}
	G_PRINTF( "GAC iteration finish\r\n"	);
	GIP_free( temp );
GIP_exit:
	MEYER_clear( &meyer );
}
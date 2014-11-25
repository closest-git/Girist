/*
	1 Fourth order PDE[REF] 
		[REF]:	Fourth order partial differential equations on general geometries
	2 lenda is not too large
		[REF]:	Image decomposition and restoration using total variation minimization and the H-1 norm
	3 staircase reduction
		[REF]:	Image decomposition combining staircase reduction and texture extraction

	v0.1	8/4/2008
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
#include "OSV_model.h"
#include "GIP_fmm.h"
#include "../../lib/gss_6_demo.h"

static int osv_step=0;

/*
	验证矩阵的某些特性

	v0.1	cys
		8/13/2008	
*/
void OSV_verify_mat( OSV_MODEL *hOSV,double *temp,double alpha,double *F )	{
	int dim=hOSV->A_dim,*ptr=hOSV->A_ptr,*ind=hOSV->A_ind,isCheckRes=1,M=hOSV->M,N=hOSV->N,MN;
	int isFind,i,j,k,r;
	double *val=hOSV->A_val,a,sum,off;

	MN=M*N;
/*	for( i = 0; i < dim; i++ )	{
		sum=0.0;
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{			
			r = ind[j];
			if( r > i )	{
				isFind = 0;
				for( k = ptr[r]; k < ptr[r+1];	k++ )	{
					if( ind[k]==i )	{
						if( F[i] != 0.0 )
						{	off = val[j]/alpha-val[k]/F[i];			ASSERT( off==0.0 );		}
						isFind=1;				
						break;
					}
				}
				ASSERT( isFind==1 );
			}
		}
	}
*/
	for( i = 0; i < 2*MN; i++ )			temp[i]=0;
	for( i = MN*2; i < dim; i++ )	{
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{			
			r = ind[j];
			if( r > 2*MN )	continue;
			temp[r]++;
		}
	}
	for( i = 0; i < 2*MN; i++ )	ASSERT( temp[i]==3 );
}

/*
	[REF] Image decomposition and restoration using total variation minimization and the H-1 norm
		公式(2.5)

	temp[4*BLK_shift]

	v0.1	cys
		8/28/2008	
*/
double OSV_get_g( OSV_MODEL *hOSV,double *g_arr,int flag )		{
	int ldu=hOSV->ldu,M=hOSV->M,N=hOSV->N,i,j,r,c,pos,BLK_shift;
	double *g,g_norm=0,fx,fy,a,*u=hOSV->u,*f=hOSV->f,lenda=hOSV->lenda,beta=hOSV->beta,*w1_,*w2_,*v_;
	GIP_FD_METHOD fdm_w=GIP_SCHEME_BACKWARD;

	ASSERT( BIT_TEST(hOSV->fdm,GIP_SCHEME_FORWARD) );

	BLK_shift=M*N;
	w1_=hOSV->temp;		w2_=w1_+BLK_shift;		v_=w2_+BLK_shift;
//	g=GIP_alloc( sizeof(double)*BLK_shift );
	if( g_arr==GIP_NULL )
		g=v_+BLK_shift;	
	else
		g=g_arr;

	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			fx = FD_R1_0( c==N-1,hOSV->fdm,u,pos,1 );
			fy = FD_R1_0( r==M-1,hOSV->fdm,u,pos,ldu ); 
			a = sqrt( fx*fx+fy*fy+beta );
			w1_[pos]=fx/a;		w2_[pos]=fy/a;
		}
	}
	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
//			fx = FD_R1_0( c==N-1,hOSV->fdm,w1_,pos,1 );
//			fy = FD_R1_0( r==M-1,hOSV->fdm,w2_,pos,ldu ); 
			fx = FD_R1_0( c==0,fdm_w,w1_,pos,1 );
			fy = FD_R1_0( r==0,fdm_w,w2_,pos,ldu );
			v_[pos] = fx+fy;
		}
	}
//	off_max=0.0;
	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			if( pos==22 )
				a=0.0;
			a = FD_LAPLAS_0( hOSV->tag[pos],hOSV->fdm,v_,pos,ldu );
			g[pos]=2*lenda*(u[pos]-f[pos])+a;
			g_norm += g[pos]*g[pos];
		}
	}
	g_norm = sqrt(g_norm);

//	GIP_free(g);

	return g_norm;
}

/*
	[REF] Modeling Textures with Total Variation Minimization and Oscillating Patterns in Image Processing

	注意：
		1 A与scheme要匹配

	v0.1	cys
		8/7/2008	
*/
void OSV_update_equations( OSV_MODEL *hOSV )	{
	int dim=hOSV->A_dim,*ptr=hOSV->A_ptr,*ind=hOSV->A_ind,ldu=hOSV->ldu,M=hOSV->M,N=hOSV->N;
	int i,j,r,c,pos,r_D,A_nz,pos_D,alg=hOSV->alg,BLK_shift,nz_D,scheme;
	double *val=hOSV->A_val,*u=hOSV->u,*k=hOSV->k,*w=hOSV->w,*f=hOSV->f,*rhs=hOSV->rhs,*temp=hOSV->temp,lenda=hOSV->lenda,beta=hOSV->beta;
	double s,a,b,a_max,a_min,*ksi,fx,fy,res,s1,s2,u_1,u_N;
	GIP_FD_METHOD fdm_w=GIP_SCHEME_BACKWARD;

	ASSERT( BIT_TEST(hOSV->fdm,GIP_SCHEME_FORWARD) );
	
	scheme=GIP_FD_SCHEME(hOSV->fdm);			ASSERT( scheme==GIP_SCHEME_FORWARD );
	BLK_shift=M*N;
	ksi=temp+dim;		
//	F_=temp+BLK_shift;
	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			fx = FD_R1_0( c==N-1,hOSV->fdm,u,pos,1 );
			fy = FD_R1_0( r==M-1,hOSV->fdm,u,pos,ldu );
			ksi[pos] = sqrt( fx*fx+fy*fy+beta );
		}
	}
	A_nz=0;
//	|U|				
	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			val[A_nz++]=2*lenda;
			if( r>0 )	val[A_nz++] = r==M-1 ? 2.0 : 1.0;
			if( c>0 )	val[A_nz++] = c==N-1 ? 2.0 : 1.0;
			val[A_nz++]=-4.0;
			if( c<N-1 )	val[A_nz++] = c==0 ? 2.0 : 1.0;
			if( r<M-1 )	val[A_nz++] = r==0 ? 2.0 : 1.0;
 			ASSERT( ptr[pos+1]==A_nz );
		} 
	}
//	| V |			
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			val[A_nz++]=1.0;
/*			if( c==N-1 )		val[A_nz++]=-1.0;
			val[A_nz++]=1.0;
			if( c<N-1 )		val[A_nz++]=-1.0;
			if( r==M-1 )		val[A_nz++]=-1.0;
			val[A_nz++]=1.0;
			if( r<M-1 )		val[A_nz++]=-1.0; */
			if( c>0 )		val[A_nz++]=1.0;
			val[A_nz++]=-1.0;
			if( c==0 )		val[A_nz++]=1.0;
			if( r>0 )		val[A_nz++]=1.0;
			val[A_nz++]=-1.0;
			if( r==0 )		val[A_nz++]=1.0; 
			ASSERT( ptr[pos+BLK_shift+1]==A_nz );
		}
	}
//	| W1 |			
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			s = w[pos]/ksi[pos];
			if( r==M-1 )		val[A_nz++]=s*(u[pos-ldu]-u[pos]);
			if( c==N-1 )		val[A_nz++]=s*(u[pos-1]-u[pos])-1;
			u_1=c==N-1?u[pos-1]:u[pos+1];
			u_N=r==M-1?u[pos-ldu]:u[pos+ldu];
			val[A_nz++]=s*(2*u[pos]-u_1-u_N)+1;
			if( c<N-1 )		val[A_nz++]=s*(u[pos+1]-u[pos])-1;
			if( r<M-1 )		val[A_nz++]=s*(u[pos+ldu]-u[pos]);
			val[A_nz++]=ksi[pos];
			ASSERT( ptr[pos+BLK_shift*2+1]==A_nz );
		}
	}
//	| W2 |			
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			s = w[pos+BLK_shift]/ksi[pos];
			if( r==M-1 )		val[A_nz++]=s*(u[pos-ldu]-u[pos])-1;
			if( c==N-1 )		val[A_nz++]=s*(u[pos-1]-u[pos]);
			u_1=c==N-1?u[pos-1]:u[pos+1];
			u_N=r==M-1?u[pos-ldu]:u[pos+ldu];
			val[A_nz++]=s*(2*u[pos]-u_1-u_N)+1;
			if( c<N-1 )		val[A_nz++]=s*(u[pos+1]-u[pos]);
			if( r<M-1 )		val[A_nz++]=s*(u[pos+ldu]-u[pos])-1;
			val[A_nz++]=ksi[pos];
			ASSERT( ptr[pos+BLK_shift*3+1]==A_nz );
		}
	}

//	OSV_verify_UU( hOSV,temp,alpha,F_ );

	a_max=DBL_MIN;		a_min=DBL_MAX;
	for( i = 0; i < A_nz; i++ )		{
		a = fabs(val[i]); 
		a_max=MAX(a_max,a);		a_min=MIN(a_min,a);	
	}

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			if( pos==21 )	{
				a=0.0;
			}
			a=FD_LAPLAS_0( hOSV->tag[pos],hOSV->fdm,k,pos,ldu );
			rhs[pos]=-( 2*lenda*(u[pos]-f[pos])+a );			
//			fx = FD_R1_0( c==N-1,hOSV->fdm,w,pos,1 );
//			fy = FD_R1_0( r==M-1,hOSV->fdm,w+BLK_shift,pos,ldu );
			fx = FD_R1_0( c==0,fdm_w,w,pos,1 );
			fy = FD_R1_0( r==0,fdm_w,w+BLK_shift,pos,ldu );
			rhs[pos+BLK_shift]=-( k[pos]-fx-fy );		//k-w1x-w2y
			fx = FD_R1_0( c==N-1,hOSV->fdm,u,pos,1 );
			fy = FD_R1_0( r==M-1,hOSV->fdm,u,pos,ldu );
			rhs[pos+BLK_shift*2]=-(-fx+w[pos]*sqrt(fx*fx+fy*fy+beta) );				
			rhs[pos+BLK_shift*3]=-(-fy+w[pos+BLK_shift]*sqrt(fx*fx+fy*fy+beta) );				
		}
	}

	res = 0.0;
	for( i = 0; i < dim; i++ )	res += rhs[i]*rhs[i];
	hOSV->nn_res = sqrt(res);

	if( 0 )	{
		char sPath[80];
		_stprintf( sPath,"H:\\GIPack\\data\\OSV_%d.mtx",dim );
		Equation_output( M,N,dim,ptr,ind,val,rhs,sPath );
	}

}

/*
	[REF] A nonlinear primal-dual method for total variation-based image restoration (1999)

	注意：
		1 hOSV->ldu=N
		2 A与scheme要匹配,A采用压缩行格式存储

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/17/2008	
*/
void OSV_init_( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,OSV_MODEL *hOSV )	{
	int M,N,MN,r,c,pos,A_nz,i,setting[32],ret=GIP_OK,ldu,BLK_shift,t_len,scheme;
	GIP_FLOATING *val=GIP_NULL;
	double *g,a,s,*f,fx,fy,beta,*temp;
	GIP_FD_METHOD fdm_w;

	osv_step=0;

	memset( hOSV,0x0,sizeof(OSV_MODEL) );
	GIP_import_IplImage( hIMG,&M,&N,&val,0x0 );

	hOSV->alg = OSV_SS_DIRECT;			//OSV_SS_ITER;			
	hOSV->fdm = GIP_SCHEME_FORWARD;			fdm_w=GIP_SCHEME_BACKWARD;
	hOSV->lenda=0.01;
	hOSV->beta = beta = 1.0e-4;
	hOSV->M=M;				hOSV->N=N;
	MN=M*N;
	hOSV->ldu=N;			hOSV->shift=0;
	t_len=MAX((M+2)*(N+2)*8,N);
	hOSV->temp=(double*)GIP_alloc( sizeof(double)*t_len );			temp=hOSV->temp;
	hOSV->f=(double*)GIP_alloc( sizeof(double)*(MN) );
	hOSV->g=(double*)GIP_alloc( sizeof(double)*(MN) );
	hOSV->u=(double*)GIP_alloc( sizeof(double)*MN*4 );
	hOSV->k=hOSV->u+MN;
	hOSV->w=hOSV->u+MN*2;
	hOSV->tag=(int*)GIP_alloc( sizeof(int)*(M+2)*(N+2) );
	memset( hOSV->tag,0x0,sizeof(int)*(M+2)*(N+2) );
	ldu=hOSV->ldu;			
	BLK_shift=M*N;
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			hOSV->f[pos] = val[pos];
			hOSV->u[pos] = val[pos];
			if( c==0 )		BIT_SET( hOSV->tag[pos],GIP_BD_LEFT );
			if( c==N-1 )	BIT_SET( hOSV->tag[pos],GIP_BD_RIGHT );
			if( r==0 )		BIT_SET( hOSV->tag[pos],GIP_BD_DOWN );
			if( r==M-1 )	BIT_SET( hOSV->tag[pos],GIP_BD_UP );
			fx = FD_R1_0( c==N-1,hOSV->fdm,val,pos,1 );
			fy = FD_R1_0( r==M-1,hOSV->fdm,val,pos,ldu );
//			fx=(c<N-1) ? val[pos+1]-val[pos] : val[pos-1]-val[pos];
//			fy=(r<M-1) ? val[pos+ldu]-val[pos] : val[pos-ldu]-val[pos];
			hOSV->w[pos] = fx/sqrt(fx*fx+fy*fy+beta);					//w1			
			hOSV->w[pos+MN] = fy/sqrt(fx*fx+fy*fy+beta);				//w2		
		}
	}
	for( r = 0; r < M; r++ )	{					//k=div(w1,w2)
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
//			fx = FD_R1_0( c==N-1,hOSV->fdm,hOSV->w,pos,1 );
//			fy = FD_R1_0( r==M-1,hOSV->fdm,hOSV->w+MN,pos,ldu );
			fx = FD_R1_0( c==0,fdm_w,hOSV->w,pos,1 );
			fy = FD_R1_0( r==0,fdm_w,hOSV->w+MN,pos,ldu );
			hOSV->k[pos]=fx+fy;
		}
	}

	hOSV->A_dim=MN*4;			
	hOSV->m_type=0;		//unsym
	A_nz=hOSV->A_dim*6;			
	hOSV->rhs=(double*)GIP_alloc( sizeof(double)*(hOSV->A_dim) );
	hOSV->A_ptr=(int*)GIP_alloc( sizeof(int)*(hOSV->A_dim+1) );
	hOSV->A_ind=(int*)GIP_alloc( sizeof(int)*(A_nz) );
	hOSV->A_val=(double*)GIP_alloc( sizeof(double)*(A_nz) );
	hOSV->A_ptr[0]=0;			A_nz=0;
	ASSERT( hOSV->m_type!=11 && hOSV->m_type!=12 );
	scheme=GIP_FD_SCHEME(hOSV->fdm);			ASSERT( scheme==GIP_SCHEME_FORWARD );
//	|U|				
	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			hOSV->A_ind[A_nz++]=pos;
			if( r>0 )	hOSV->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu )+BLK_shift;
			if( c>0 )	hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu )+BLK_shift;
						hOSV->A_ind[A_nz++]=pos+BLK_shift;
			if( c<N-1 )	hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu )+BLK_shift;
			if( r<M-1 )	hOSV->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu )+BLK_shift;
 			hOSV->A_ptr[pos+1]=A_nz;
		} 
	}
//	| V |			
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			hOSV->A_ind[A_nz++]=pos+BLK_shift;
/*			if( c==N-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu )+BLK_shift*2;
			hOSV->A_ind[A_nz++]=pos+BLK_shift*2;
			if( c<N-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu )+BLK_shift*2;
			if( r==M-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu )+BLK_shift*3;
			hOSV->A_ind[A_nz++]=pos+BLK_shift*3;
			if( r<M-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu )+BLK_shift*3;
			hOSV->A_ptr[pos+BLK_shift+1]=A_nz;*/
			if( c>0 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu )+BLK_shift*2;
			hOSV->A_ind[A_nz++]=pos+BLK_shift*2;
			if( c==0 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu )+BLK_shift*2;
			if( r>0 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu )+BLK_shift*3;
			hOSV->A_ind[A_nz++]=pos+BLK_shift*3;
			if( r==0 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu )+BLK_shift*3;
			hOSV->A_ptr[pos+BLK_shift+1]=A_nz;
		}
	}
//	| W1 |			
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			if( r==M-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu );
			if( c==N-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu );
			hOSV->A_ind[A_nz++]=pos;
			if( c<N-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu );
			if( r<M-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu );
			hOSV->A_ind[A_nz++]=pos+BLK_shift*2;
			hOSV->A_ptr[pos+BLK_shift*2+1]=A_nz;
		}
	}
//	| W2 |			
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			if( r==M-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu );
			if( c==N-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu );
			hOSV->A_ind[A_nz++]=pos;
			if( c<N-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu );
			if( r<M-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu );
			hOSV->A_ind[A_nz++]=pos+BLK_shift*3;
			hOSV->A_ptr[pos+BLK_shift*3+1]=A_nz;
		}
	}
	hOSV->A_nz = A_nz;			ASSERT( A_nz<=hOSV->A_dim*6 );

	OSV_update_equations( hOSV );	
	hOSV->g_norm=OSV_get_g( hOSV,hOSV->g,0x0 );

	hOSV->hSS=GIP_NULL;
	if( BIT_TEST( hOSV->alg,OSV_SS_DIRECT) )	{
		int *colptr,*rowind;
		GIP_FLOATING *data;
		GIP_CRS2CCS( hOSV->A_dim,hOSV->A_dim,GIP_NULL,hOSV->A_ptr,hOSV->A_ind,&data,&colptr,&rowind,hOSV->temp );
		for( i = 0; i < 32; i++ )	setting[i]=0.0;
		if( (ret=GSS_init_id(hOSV->A_dim,hOSV->A_dim,colptr,rowind,GIP_NULL,hOSV->m_type,setting) ) != GRUS_OK )	
		{	GIP_ERROR( "TVM_init",ret,"\tERROR at init GSS solver. ERROR CODE:%d\r\n" );		}
		hOSV->hSS = GSS_symbol_id( hOSV->A_dim,hOSV->A_dim,colptr,rowind,hOSV->A_val );	
		if( hOSV->hSS == NULL )	
		{	GIP_ERROR( "TVM_init",-2,"\tERROR at SYMBOLIC ANALYSIS.\r\n" );		}
		GIP_free(colptr);		GIP_free(rowind);
	}else	{
		hOSV->hSS=GIP_Krylov_init_( hOSV->A_dim,15,0x0,hOSV->A_ptr,hOSV->A_ind,hOSV->A_val  );
		ASSERT( hOSV->hSS!=GIP_NULL );
	}

GIP_exit:	
	GIP_free( val );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void OSV_clear( OSV_MODEL *hOSV )	{
	if( hOSV == GIP_NULL )	return;

	hOSV->M=-1;		hOSV->N=-1;
	if( hOSV->A_ptr != GIP_NULL )	
	{	GIP_free( hOSV->A_ptr );	hOSV->A_ptr=GIP_NULL; }
	if( hOSV->A_ind != GIP_NULL )	
	{	GIP_free( hOSV->A_ind );	hOSV->A_ind=GIP_NULL; }
	if( hOSV->A_val != GIP_NULL )	
	{	GIP_free( hOSV->A_val );	hOSV->A_val=GIP_NULL; }
//	if( hOSV->w != GIP_NULL )	
//	{	GIP_free( hOSV->w );	hOSV->w=GIP_NULL; }
	if( hOSV->f != GIP_NULL )	
	{	GIP_free( hOSV->f );	hOSV->f=GIP_NULL; }
	if( hOSV->g != GIP_NULL )	
	{	GIP_free( hOSV->g );	hOSV->g=GIP_NULL; }
	if( hOSV->tag != GIP_NULL )	
	{	GIP_free( hOSV->tag );	hOSV->tag=GIP_NULL; }
	if( hOSV->temp != GIP_NULL )	
	{	GIP_free( hOSV->temp );	hOSV->temp=GIP_NULL; }
	if( hOSV->hSS!=GIP_NULL )	{
		if( BIT_TEST( hOSV->alg,OSV_SS_DIRECT) )
			GSS_clear_id( hOSV->hSS );
//		else
//			GIP_Krylov_clear_( (GIP_KRYLOV_SOLVER*)hOSV->hSS );
		hOSV->hSS=GIP_NULL;
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/18/2008	
*/
void OSV_step( OSV_MODEL *hOSV,double dT,double *temp )	{
	int M=hOSV->M,N=hOSV->N,MN=M*N,i,j,r,c,ldu=hOSV->ldu,*tag=hOSV->tag,ret=GIP_OK;
	int dim=hOSV->A_dim,*ptr=hOSV->A_ptr,*ind=hOSV->A_ind,isCheckRes=1;
	double *rhs=hOSV->rhs,u_max,u_min,w_max,w_min,rhs_norm=0.0,*val=hOSV->A_val,res,a,*data;
	int *colptr,*rowind;

//use GSS to solver A*u{n+1}=f(u{n})
	GIP_CRS2CCS( hOSV->A_dim,hOSV->A_dim,val,ptr,ind,&data,&colptr,&rowind,hOSV->temp );
//	if( osv_step%5==0 )	{
		ret = GSS_numeric_id( hOSV->A_dim,hOSV->A_dim,colptr,rowind,data,hOSV->hSS );
		if( ret != GRUS_OK )	{
			hOSV->hSS=NULL;		//±ØÐëÉèÖÃÎªNULL,GSSÒÑ×Ô¶¯ÊÍ·ÅÄÚ´æ
			GIP_ERROR( "TVM_step",ret,"\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n" );
		}
//	}

	memcpy( temp,rhs,sizeof(double)*dim );
	if( isCheckRes==1 )		{
		for( i = 0; i < dim; i++ )	
		{	rhs_norm += (temp[i]*temp[i]);		}
		rhs_norm = sqrt(rhs_norm);	
	}
	GSS_solve_id( hOSV->hSS,hOSV->A_dim,hOSV->A_dim,colptr,rowind,data,rhs );
	GIP_free(colptr);		GIP_free(rowind);
	if( isCheckRes==1 )	{
		double *x=rhs;
		res = 0.0;
		u_max=0.0,		u_min=DBL_MAX;			
		w_max=0.0,		w_min=DBL_MAX;			
		for( i = 0; i < dim; i++ )	{
			a = fabs(x[i]);
			if( i>= MN*2 )	{
				u_max=MAX( u_max,a );			u_min=MIN( u_min,a );
			}else	{
				w_max=MAX( w_max,a );			w_min=MIN( w_min,a );
			}
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				r = ind[j];
				temp[i] -= val[j]*x[r];
			}		
		}
		for( i = 0; i < dim; i++ )	
		{	res += (temp[i]*temp[i]);		}
		res = sqrt(res);	
		if( res/rhs_norm > 1.0e-8 )
			G_PRINTF( "\r\n\t\t|Ax-b|=%g,|b|=%g,|du|__=%g,|dw|__=%g\r\n",res,rhs_norm,u_max,w_max );	
//		G_PRINTF( "\t\t|Ax-b|=%g,|b|=%g,u_max=%g,u_min=%g\r\n",res,rhs_norm,u_max,u_min );	
	}

GIP_exit:
	;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/25/2008	
*/
void OSV_step_iter( OSV_MODEL *hOSV )	{
	int M=hOSV->M,N=hOSV->N,MN=M*N,i,j,ldu=hOSV->ldu,*tag=hOSV->tag,ret=GIP_OK,t;
	int dim=hOSV->A_dim,*ptr=hOSV->A_ptr,*ind=hOSV->A_ind,isCheckRes=1,r,c,pos;
	double *rhs=hOSV->rhs,*x_0=hOSV->w,u_max,u_min,rhs_norm=0.0,*val=hOSV->A_val,a,d,*x,*temp=hOSV->temp;
	GIP_KRYLOV_SOLVER* hKS;

	x=temp+dim;
//	memcpy( x,hOSV->w,sizeof(double)*dim );
	memset( x,0x0,sizeof(double)*dim );
	if( isCheckRes==1 )		{
		for( i = 0; i < dim; i++ )	
		{	rhs_norm += (rhs[i]*rhs[i]);		}
		rhs_norm = sqrt(rhs_norm);	
	}

//	hKS=GIP_Krylov_init_( dim,10,0x0,ptr,ind,val  );
	hKS=(GIP_KRYLOV_SOLVER*)(hOSV->hSS);			//GIP_Krylov_init_( dim,15,0x0,ptr,ind,val  );

	GIP_Krylov_core_( hKS,x,rhs,GIP_ROW_MAJOR  );

	if( isCheckRes==1 )	{
		double res = 0.0;
		u_max=DBL_MIN,		u_min=DBL_MAX;
		memcpy( temp,rhs,sizeof(double)*dim );
		for( i = 0; i < dim; i++ )	{
			a = fabs(x[i]);
			u_max=MAX( u_max,a );			u_min=MIN( u_min,a );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				c = ind[j];
				temp[i] -= val[j]*x[c];
			}		
		}
		for( i = 0; i < dim; i++ )	
		{	res += (temp[i]*temp[i]);		}
		res = sqrt(res);	
//		G_PRINTF( "\t\t|Ax-b|=%g,|b|=%g,u_max=%g,u_min=%g\r\n",res,rhs_norm,u_max,u_min );	
	}

	memcpy( rhs,x,sizeof(double)*dim );

}

/*	
	采用近似估计:	|du|+|dk|/2
	[REF]:	Image decomposition and restoration using total variation minimization and the H-1 norm

	ISSUE-NEEDED
	1 需要f-u在H-1空间内的norm
	is exactly the seminorm of H−1(Ω), the dual to H1(Ω) (L.C. Evans, Partial Differential Equations)

	v0.1	cys
		9/5/2008
*/
double  OSV_Energy_( OSV_MODEL *hOSV,double *u_cur,int flag )	{
	int r,c,ldu=hOSV->ldu,M=hOSV->M,N=hOSV->N,pos,pos_D,BLK_shift=M*N;
	double res=0.0,energy=0.0,a,fx,fy,x,y,a_max;
	double *f=hOSV->f,*k=hOSV->k,*w1_,*w2_,*k_,lenda=hOSV->lenda,beta=hOSV->beta;
	GIP_FD_METHOD fdm_w=GIP_SCHEME_BACKWARD;

	ASSERT( BIT_TEST(hOSV->fdm,GIP_SCHEME_FORWARD) );
	w1_=hOSV->temp;		w2_=w1_+BLK_shift;		k_=w2_+BLK_shift;
	if( 1 )	{
		for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
			for( c = 0; c < N; c++ )	{
				pos = GIP_RC2POS( r,c,ldu );
				fx = FD_R1_0( c==N-1,hOSV->fdm,u_cur,pos,1 );
				fy = FD_R1_0( r==M-1,hOSV->fdm,u_cur,pos,ldu ); 
				a = sqrt( fx*fx+fy*fy+beta );
				w1_[pos]=fx/a;		w2_[pos]=fy/a;
			}
		}
		for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
			for( c = 0; c < N; c++ )	{
				pos = GIP_RC2POS( r,c,ldu );
				fx = FD_R1_0( c==0,fdm_w,w1_,pos,1 );
				fy = FD_R1_0( r==0,fdm_w,w2_,pos,ldu );
				k_[pos] = fx+fy;
			}
		}
	}

	a_max = 0.0;
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos_D = GIP_RC2POS( r,c,ldu );
			fx = FD_R1_0( c==N-1,hOSV->fdm,u_cur,pos_D,1 );
			fy = FD_R1_0( r==M-1,hOSV->fdm,u_cur,pos_D,ldu );
			energy += sqrt(fx*fx+fy*fy+beta);
			if( 1 )	{
				fx = FD_R1_0( c==N-1,hOSV->fdm,k_,pos_D,1 );
				fy = FD_R1_0( r==M-1,hOSV->fdm,k_,pos_D,ldu );
				energy += 0.25*(fx*fx+fy*fy)/lenda ;
			}
		}
	}	

	return energy;
}


/*
	k=div(g1,g2)

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/11/2008	
*/
void OSV_get_v( OSV_MODEL *hOSV,double *k )	{
	int ldu=hOSV->ldu,M=hOSV->M,N=hOSV->N,r,c,pos,poo=-1;
	double *f=hOSV->f,*u=hOSV->u;

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			k[pos] = f[pos]-u[pos];
		}
	}	;
}

/*
	[REF] A nonlinear primal-dual method for total variation-based image restoration (1999)

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/19/2008	
*/
void OSV_newton_step( OSV_MODEL *hOSV,double *temp )	{
	int i,M=hOSV->M,N=hOSV->N,MN,pos,r,c,ldu=hOSV->ldu,no;
	double *rhs=hOSV->rhs,*u=hOSV->u,*k=hOSV->k,*w=hOSV->w,*g=hOSV->g,sp=1.0,sd=1.0,rou=0.9,alpha=0.0,s,p_2,q,a1,a2;
	double en_0,en,*u_,*u_0,dot,ksi=0.5,off,off_max,u_max,u_min,k_sum;

//	OSV_update_La( hOSV,temp );

	if( BIT_TEST( hOSV->alg,OSV_SS_DIRECT ) )
		OSV_step( hOSV,0.0,temp );
	else if( BIT_TEST( hOSV->alg,OSV_SS_ITER ) )
		OSV_step_iter( hOSV );

	MN=M*N;			ASSERT( MN*4==hOSV->A_dim );
/*
	u_0 = GIP_alloc(sizeof(double)*3*MN );
	u_ = u_0+MN;	
	for( i=0; i < MN; i++ )	u_0[i]=u[i];
	en_0 = OSV_Energy_( hOSV,u_0,0x0 );
	sp=2.0;
	do	{
		sp /= 2;
		for( i=0; i < MN; i++ )	u_[i]=u_0[i]+rhs[i]*sp;
		en = OSV_Energy_( hOSV,u_,0x0 );
		for( i=0; i < MN; i++ )	dot=rhs[i]*g[i];
		en -= ksi*sp*dot;
	}while( en >= en_0 );
*/

	for( i = 0; i < MN; i++ )	{
		a1=rhs[2*MN+i];		a2=rhs[3*MN+i];
		if( a1==0.0 && a2==0.0 )
			continue;

		s=a1*a1+a2*a2;
		p_2=(w[i]*a1+w[MN+i]*a2)/s;
		q=(w[i]*w[i]+w[MN+i]*w[MN+i]-1.0)/s;		
//		q=fabs(q)<FLT_EPSILON ? 0.0:q;			ASSERT( q<=0.0 );	
		q=MIN( q,0.0 );
		s=sqrt(-q+p_2*p_2)-p_2;			ASSERT( s>=0);	//s=sqrt(-q+(p*p/4))-p/2
//		if( s>0.0 )
			sd = MIN( sd,s );				//sd=sup{s:|w+sdw|<1}
	}
	sd *= rou;
	k_sum=0;
	for( i = 0; i < MN; i++ )	{
		k[i] += rhs[MN+i]*sd;		
		k_sum+=k[i];
		w[i]+=rhs[2*MN+i]*sd;			
		w[i+MN]+=rhs[3*MN+i]*sd;	
		s=sqrt( w[i]*w[i]+w[i+MN]*w[i+MN] );		ASSERT( s <= 1.0 );
	}
	u_max=DBL_MIN,			u_min=DBL_MAX;
	for( i = 0; i < MN; i++ )	{
		u[i]+=rhs[i]*sp;	
		u_max=MAX(u_max,u[i]);
		u_min=MIN(u_min,u[i]);
	}

	OSV_update_equations( hOSV );
	hOSV->g_norm=OSV_get_g( hOSV,g,0x0 );
	off_max=0.0;
	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			off = fabs(g[pos]-hOSV->rhs[pos]);
			if( off > off_max )	{
				off_max=off;		no=pos;
			}
		}
	}
//	GIP_free( u_0 );
	G_PRINTF( "\r\n  ****  sp=%.3g,sd=%.3g,K_sum=%.3g,off=%.2g(%d),contrast=%d(%d,%d)\r\n",sp,sd,k_sum,
		off_max,no,(int)(u_max-u_min),(int)u_max,(int)u_min );		

}

/*
	A testing code 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/5/2008	
*/
void OSV_MODEL_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName )	{
	double *temp,energy,res,tol;
	OSV_MODEL cgm,*hOSV=&cgm;
	int M,N,isTrace=0,nKey,t_step=1;
	char sTitle[80];

	OSV_init_( hIMG,hIMG_mask,hOSV );
	M=hOSV->M,		N=hOSV->N;
	temp=hOSV->temp;
	res=hOSV->nn_res;
	tol = res*1.0e-4;
	G_PRINTF( "OSV:\tM=%d,N=%d,lenda=%g,beta=%g,alg=%x,res=%g\r\n",M,N,hOSV->lenda,hOSV->beta,hOSV->alg,res );		

//	res = OSV_update_res( hOSV,temp );
//	GIP_show_IplImage_d( M,N,sWndName,hOSV->z,hOSV->ldu,0x0 );
	while( res>tol || BIT_TEST(hOSV->flag,GIP_SOLVER_CONVERGE) )	{
			_stprintf( sTitle,".\\trace\\u_%d.bmp",osv_step );		
			GIP_save_IplImage_d( M,N,sTitle,(hOSV->u+hOSV->shift),hOSV->ldu,0x0 );
			OSV_get_v( hOSV,temp );
			_stprintf( sTitle,".\\trace\\v_%d.bmp",osv_step );		
			GIP_save_IplImage_d( M,N,sTitle,temp,hOSV->ldu,0x0 );
//			GIP_show_IplImage_d( M,N,sWndName,temp,hOSV->ldu,0x0 );
		if( isTrace==1 && osv_step%t_step==0 )	{
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
//			GIP_save_IplImage_d( M,N,sTitle,(temp+hOSV->shift),hOSV->ldu,0x0 );
		}
		OSV_newton_step( hOSV,temp );
		res = hOSV->nn_res;
		G_PRINTF( "OSV_%d: energy=%g,res=%g,|g|=%g\r\n",osv_step,
			OSV_Energy_(hOSV,hOSV->u,0x0),hOSV->nn_res,hOSV->g_norm );

//		res = OSV_update_res( hOSV,temp );
		if( BIT_TEST(hOSV->flag,GIP_SOLVER_CONVERGE) )
			break;
		osv_step++;
	}
	G_PRINTF( "OSV iteration finish\r\n"	);
GIP_exit:
	OSV_clear( hOSV );
}
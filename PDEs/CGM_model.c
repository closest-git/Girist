/*
	[REF] Modeling Textures with Total Variation Minimization and Oscillating Patterns in Image Processing
	[REF] Image denoising and decomposition with total variation. minimization and oscillatory functions

	ISSUE-NEEDED:
	1、CGM_A_AT,CGM_A_BOUNDARY的区别
		理论上似乎根据boundary condition修正A更合理，但数值求解有困难。
	2、A'对应于前向差分，而A对应于负的后向差分


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
#include "CGM_model.h"
#include "GIP_fmm.h"
#include "../../lib/gss_6_demo.h"

static int cgm_step=0;

/*
	The matrix A is symmetric 

	v0.1	cys
		8/13/2008	
*/
void CGM_verify_UU( CGM_MODEL *hCGM,double *temp,double alpha,double *F )	{
	int dim=hCGM->A_dim,*ptr=hCGM->A_ptr,*ind=hCGM->A_ind,isCheckRes=1,M=hCGM->M,N=hCGM->N,MN;
	int isFind,i,j,k,r;
	double *val=hCGM->A_val,a,sum,off;

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
	[REF] Modeling Textures with Total Variation Minimization and Oscillating Patterns in Image Processing

	v0.1	cys
		8/7/2008	
*/
double CGM_get_g( CGM_MODEL *hCGM,double *g_arr,int flag )		{
	int ldu=hCGM->ldu,M=hCGM->M,N=hCGM->N,i,j,r,c,pos,BLK_shift,no;
	double *g,g_norm=0,fx,fy,a,*f=hCGM->f,*u=hCGM->u,*w1_,*w2_,*t_Ag,alpha=hCGM->alpha,beta=hCGM->beta,off,off_max;
	GIP_FD_METHOD fdm_w;

	fdm_w=GIP_SCHEME_BACKWARD;	ASSERT( BIT_TEST(hCGM->fdm,GIP_SCHEME_FORWARD) );
	BLK_shift=M*N;
	w1_=hCGM->temp;		w2_=w1_+BLK_shift;		t_Ag=w2_+BLK_shift;
	if( g_arr==GIP_NULL )
		g=w2_+BLK_shift;	
	else
		g=g_arr;

	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			fx = FD_R1_0( c==N-1,hCGM->fdm,u,pos,1 );
			fy = FD_R1_0( r==M-1,hCGM->fdm,u,pos,ldu );
			a = sqrt(fx*fx+fy*fy+beta);
			w1_[pos]=fx/a;		w2_[pos]=fy/a;			
		}
	}
	if( 0 )	{
		int dim=hCGM->A_dim,*ptr=hCGM->A_ptr,*ind=hCGM->A_ind;
		double *val=hCGM->A_val;
		for( i = 0; i < BLK_shift; i++ )	t_Ag[i]=0.0;
		for( i = 0; i < BLK_shift; i++ )	{
			for( j = ptr[2*i]; j < ptr[2*i+1]; j++ )	{
				r = ind[j];
				if( r < BLK_shift*2 )	continue;
				t_Ag[r-BLK_shift*2] += val[j]*w1_[i];
			}
			for( j = ptr[2*i+1]; j < ptr[2*i+2]; j++ )	{
				r = ind[j];
				if( r < BLK_shift*2 )	continue;
				t_Ag[r-BLK_shift*2] += val[j]*w2_[i];
			}
		}
	}
	off_max=0.0;
	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
//			fx = FD_R1_0( c==N-1,hCGM->fdm,w1_,pos,1 );
//			fy = FD_R1_0( r==M-1,hCGM->fdm,w2_,pos,ldu );
			fx = FD_R1_0( c==0,fdm_w,w1_,pos,1 );
			fy = FD_R1_0( r==0,fdm_w,w2_,pos,ldu );
			g[pos]=-alpha*(fx+fy)+(u[pos]-f[pos]);
			g_norm += g[pos]*g[pos];
			off = fabs((-alpha*(fx+fy)-t_Ag[pos])/fabs(t_Ag[pos]));
			if( off > off_max )	{
				off_max=off;		no=pos;
			}
		}
	}
	if( BIT_TEST( hCGM->alg,CGM_A_BOUNDARY) )
		ASSERT( off_max < FLT_EPSILON );
//	G_PRINTF( "\t  ****  off_max=%g,no=%d\r\n",off_max,no );
	g_norm = sqrt(g_norm);

	return g_norm;
}

/*
	[REF] Modeling Textures with Total Variation Minimization and Oscillating Patterns in Image Processing

	注意：
		1 A与scheme要匹配

	v0.1	cys
		8/7/2008	
*/
void CGM_update_equations( CGM_MODEL *hCGM )	{
	int dim=hCGM->A_dim,*ptr=hCGM->A_ptr,*ind=hCGM->A_ind,ldu=hCGM->ldu,M=hCGM->M,N=hCGM->N;
	int i,j,r,c,pos,r_D,A_nz,pos_D,alg=hCGM->alg,BLK_shift,nz_D,scheme;
	double *val=hCGM->A_val,*u=hCGM->u,*w=hCGM->w,*f=hCGM->f,*rhs=hCGM->rhs,*temp=hCGM->temp,alpha=hCGM->alpha,beta=hCGM->beta;
	double s,a,b,a_max,a_min,*ksi,fx,fy,res,s1,s2,u_1,u_N;
	
	scheme=GIP_FD_SCHEME(hCGM->fdm);			ASSERT( scheme==GIP_SCHEME_FORWARD );
	BLK_shift=M*N;
	ksi=temp+dim;		
//	F_=temp+BLK_shift;
	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			fx = FD_R1_0( c==N-1,hCGM->fdm,u,pos,1 );
			fy = FD_R1_0( r==M-1,hCGM->fdm,u,pos,ldu );
//			fx=(c<N-1) ? u[pos+1]-u[pos] : u[pos-1]-u[pos];
//			fy=(r<M-1) ? u[pos+ldu]-u[pos] : u[pos-ldu]-u[pos];
			ksi[pos] = sqrt( fx*fx+fy*fy+beta );
		}
	}
	A_nz=0;
//	|E|		(2MN:2MN)		
//	|A|		(MN:2MN)		
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos_D=GIP_RC2POS( r,c,ldu );
			if( BIT_TEST(hCGM->alg,CGM_A_AT) )	{
				val[A_nz++]=ksi[pos_D];			//[:2*pos_D]
				if( c==N-1 )		//boundary condition
				{	val[A_nz++]=alpha;		val[A_nz++]=-alpha;	}
				else
				{	val[A_nz++]=-alpha;		val[A_nz++]=alpha;	}
				ASSERT( A_nz==hCGM->A_ptr[2*pos_D+1] );
		
				val[A_nz++]=ksi[pos_D];			//[:2*pos_D+1]
				if( r==M-1 )		//boundary condition
				{	val[A_nz++]=alpha;		val[A_nz++]=-alpha;	}
				else
				{	val[A_nz++]=-alpha;		val[A_nz++]=alpha;	}
				ASSERT( A_nz==hCGM->A_ptr[2*pos_D+2] );
			}else if( BIT_TEST(hCGM->alg,CGM_A_BOUNDARY) )	{
				val[A_nz++]=ksi[pos_D];			//[:2*pos_D]
				if( c==1 )		val[A_nz++]=alpha;
								val[A_nz++]=-alpha;
				if( c<N-1 )		val[A_nz++]=alpha;
				ASSERT( A_nz==hCGM->A_ptr[2*pos_D+1] );
		
				val[A_nz++]=ksi[pos_D];			//[:2*pos_D+1]
				if( r==1 )		val[A_nz++]=alpha;
								val[A_nz++]=-alpha;
				if( r<M-1 )		val[A_nz++]=alpha;
				ASSERT( A_nz==hCGM->A_ptr[2*pos_D+2] );
			}
/*			val[A_nz++]=ksi[pos_D];			//[:2*pos_D]
			if( c>0 )		val[A_nz++]=alpha;
			val[A_nz++]=-alpha;	
			if( c==N-2 )	val[A_nz++]=alpha;		
			ASSERT( A_nz==hCGM->A_ptr[2*pos_D+1] );
	
			val[A_nz++]=ksi[pos_D];			//[:2*pos_D+1]
			if( r>0 )		val[A_nz++]=alpha;		
			val[A_nz++]=-alpha;	
			if( r==M-2 )	val[A_nz++]=alpha;	
			ASSERT( A_nz==hCGM->A_ptr[2*pos_D+2] );*/
		}
	}
//	| B |		(2MN:MN)		
//	| I |		( MN:MN)		
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos_D = GIP_RC2POS( r,c,ldu );
			if( r > 0 )	{
				s1=w[2*(pos_D-ldu)]/ksi[pos_D-ldu];			s2=w[2*(pos_D-ldu)+1]/ksi[pos_D-ldu];			
				val[A_nz++]=s1*(u[pos_D]-u[pos_D-ldu]);
				val[A_nz++]=s2*(u[pos_D]-u[pos_D-ldu])-1;
			}
			if( c>0 )	{
				s1=w[2*(pos_D-1)]/ksi[pos_D-1];				s2=w[2*(pos_D-1)+1]/ksi[pos_D-1];			
				val[A_nz++]=s1*(u[pos_D]-u[pos_D-1])-1;
				val[A_nz++]=s2*(u[pos_D]-u[pos_D-1]);
			}
			s1=w[2*pos_D]/ksi[pos_D];					s2=w[2*pos_D+1]/ksi[pos_D];			
			u_1= c<N-1 ? u[pos_D+1]:u[pos_D-1];				u_N= r<M-1 ? u[pos_D+ldu]:u[pos_D-ldu];
			val[A_nz++]=s1*(2*u[pos_D]-u_1-u_N)+1;
			val[A_nz++]=s2*(2*u[pos_D]-u_1-u_N)+1;
			if( c==N-2 )	{
 				s1=w[2*(pos_D+1)]/ksi[pos_D+1];				s2=w[2*(pos_D+1)+1]/ksi[pos_D+1];			
				val[A_nz++]=s1*(u[pos_D]-u[pos_D+1])-1;
 				val[A_nz++]=s2*(u[pos_D]-u[pos_D+1]);
			}
			if( r==M-2 )	{
				s1=w[2*(pos_D+ldu)]/ksi[pos_D+ldu];			s2=w[2*(pos_D+ldu)+1]/ksi[pos_D+ldu];			
				val[A_nz++]=s1*(u[pos_D]-u[pos_D+ldu]);
				val[A_nz++]=s2*(u[pos_D]-u[pos_D+ldu])-1;
			}
			val[A_nz++]=1.0;
			ASSERT( hCGM->A_ptr[pos_D+2*BLK_shift+1]==A_nz );
		}
	}

//	CGM_verify_UU( hCGM,temp,alpha,F_ );

	a_max=DBL_MIN;		a_min=DBL_MAX;
	for( i = 0; i < A_nz; i++ )		{
		a = fabs(val[i]); 
		a_max=MAX(a_max,a);		a_min=MIN(a_min,a);	
	}
	{
		for( r = 0; r < M; r++ )	{		
			for( c = 0; c < N; c++ )	{
				pos_D = GIP_RC2POS( r,c,ldu );
				if( pos_D==9361 )
				{	a = 0.0;}
				u_1= c<N-1 ? u[pos_D+1]:u[pos_D-1];				u_N= r<M-1 ? u[pos_D+ldu]:u[pos_D-ldu];
				rhs[2*pos_D]=-( ksi[pos_D]*w[2*pos_D]-(u_1-u[pos_D]) );			
				rhs[2*pos_D+1]=-( ksi[pos_D]*w[2*pos_D+1]-(u_N-u[pos_D]) );
				rhs[pos_D+2*BLK_shift]=-(u[pos_D]-f[pos_D]);				
			}
		}
		for( i = 0; i < 2*BLK_shift; i++ )	{				
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				r = ind[j];
				if( r <2*BLK_shift )		continue;
				rhs[r] += -(val[j]*w[i]);					//u1=-(alpha*A*w+u-u0)
			}
		}
	}
	res = 0.0;
	for( i = 0; i < dim; i++ )	res += rhs[i]*rhs[i];
	hCGM->nn_res = sqrt(res);

	if( 0 && cgm_step==0 )	{
		char sPath[80];
		_stprintf( sPath,"H:\\GIPack\\data\\CGM_%d_0.mtx",dim );
		Equation_output( M,N,dim,ptr,ind,val,rhs,sPath );
	}

}

/*
	[REF] A nonlinear primal-dual method for total variation-based image restoration (1999)

	注意：
		1 hCGM->ldu=N
		2 A与scheme要匹配

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/17/2008	
*/
void CGM_init_( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,CGM_MODEL *hCGM )	{
	int M,N,MN,r,c,pos,A_nz,i,setting[32],ret=GIP_OK,ldu,BLK_shift,t_len,scheme;
	GIP_FLOATING *val=GIP_NULL;
	double a,s,*f,fx,fy,beta,*temp;

	cgm_step=0;

	memset( hCGM,0x0,sizeof(CGM_MODEL) );
	GIP_import_IplImage( hIMG,&M,&N,&val,0x0 );
	GIP_conmpact_index( M,N,N,val );

	hCGM->alg = CGM_SS_DIRECT | CGM_A_AT | CGM_METHOD_MEDIAN;				//CGM_SS_ITER;	CGM_METHOD_DUAL
	hCGM->fdm = GIP_SCHEME_FORWARD;						//GIP_SCHEME_FORWARD;
	hCGM->lenda=1.0;
	hCGM->alpha = 30;			hCGM->beta = beta = 0.0001;
	hCGM->M=M;				hCGM->N=N;
	MN=M*N;
	hCGM->ldu=N;			hCGM->shift=0;
	t_len=MAX((M+2)*(N+2)*6,N);
	hCGM->temp=GIP_alloc( sizeof(double)*t_len );			temp=hCGM->temp;
	hCGM->f=(double*)GIP_alloc( sizeof(double)*(MN) );
	hCGM->g=(double*)GIP_alloc( sizeof(double)*(MN) );
//	hCGM->u=(double*)GIP_alloc( sizeof(double)*MN );
//	hCGM->w=(double*)GIP_alloc( sizeof(double)*MN*2 );
	hCGM->w=(double*)GIP_alloc( sizeof(double)*MN*3 );
	hCGM->u=hCGM->w+MN*2;
	hCGM->tag=(int*)GIP_alloc( sizeof(int)*(M+2)*(N+2) );
	memset( hCGM->tag,0x0,sizeof(int)*(M+2)*(N+2) );
	ldu=hCGM->ldu;			
	BLK_shift=M*N;
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			hCGM->f[pos] = val[pos];
			hCGM->u[pos] = val[pos];
			if( c==0 )		BIT_SET( hCGM->tag[pos],GIP_BD_LEFT );
			if( c==N-1 )	BIT_SET( hCGM->tag[pos],GIP_BD_RIGHT );
			if( r==0 )		BIT_SET( hCGM->tag[pos],GIP_BD_DOWN );
			if( r==M-1 )	BIT_SET( hCGM->tag[pos],GIP_BD_UP );
//			if( !BIT_TEST(hCGM->alg,CGM_METHOD_MEDIAN) )	{
				fx = FD_R1_0( c==N-1,hCGM->fdm,val,pos,1 );
				fy = FD_R1_0( r==M-1,hCGM->fdm,val,pos,ldu );
	//			fx=(c<N-1) ? val[pos+1]-val[pos] : val[pos-1]-val[pos];
	//			fy=(r<M-1) ? val[pos+ldu]-val[pos] : val[pos-ldu]-val[pos];
				hCGM->w[2*pos] = fx/sqrt(fx*fx+fy*fy+beta);				
				hCGM->w[2*pos+1] = fy/sqrt(fx*fx+fy*fy+beta);	
//			}
		}
	}

	if( BIT_TEST( hCGM->alg,CGM_METHOD_MEDIAN) )	{
		GIP_FD_EQUATION mat;
		FDM_init_equations_m_( M,N,ldu,hCGM->fdm,0x0,-hCGM->alpha,1.0,&mat );
		hCGM->A_dim=mat.dim;	hCGM->A_nz=mat.nz;
		GIP_CRS2CCS( hCGM->A_dim,hCGM->A_dim,mat.val,mat.ptr,mat.ind,
			&(hCGM->A_val),&(hCGM->A_ptr),&(hCGM->A_ind),hCGM->temp );
		hCGM->rhs=(double*)GIP_alloc( sizeof(double)*(hCGM->A_dim) );			ASSERT(hCGM->A_dim==M*N);
		a=0.0;
		for( r = 0; r < M; r++ )	{		
			for( c = 0; c < N; c++ )	{
				pos = GIP_RC2POS( r,c,ldu );
				hCGM->rhs[pos] = val[pos];
				a += val[pos]*val[pos];
			}
		}
		hCGM->nn_res = sqrt(a);	
	}else	{
		hCGM->A_dim=MN*3;			
		hCGM->m_type=0;		//unsym
		A_nz=hCGM->m_type==11 || hCGM->m_type==12 ? hCGM->A_dim*6 : hCGM->A_dim*6;			//È¡ÉÏÏÞ
		hCGM->rhs=(double*)GIP_alloc( sizeof(double)*(hCGM->A_dim) );
		hCGM->A_ptr=(int*)GIP_alloc( sizeof(int)*(hCGM->A_dim+1) );
		hCGM->A_ind=(int*)GIP_alloc( sizeof(int)*(A_nz) );
		hCGM->A_val=(double*)GIP_alloc( sizeof(double)*(A_nz) );
		hCGM->A_ptr[0]=0;			A_nz=0;
		ASSERT( hCGM->m_type!=11 && hCGM->m_type!=12 );
		scheme=GIP_FD_SCHEME(hCGM->fdm);			ASSERT( scheme==GIP_SCHEME_FORWARD );
	//	|E|		(2MN:2MN)		
	//	|A|		(MN:2MN)		
		for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
			for( c = 0; c < N; c++ )	{
				pos = GIP_RC2POS( r,c,ldu );
				if( BIT_TEST(hCGM->alg,CGM_A_AT) )	{
					hCGM->A_ind[A_nz++]=2*pos;
					if( c==N-1 )	hCGM->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu )+BLK_shift*2;
									hCGM->A_ind[A_nz++]=pos+BLK_shift*2;
					if( c<N-1 )		hCGM->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu )+BLK_shift*2;
 					hCGM->A_ptr[2*pos+1]=A_nz;
					hCGM->A_ind[A_nz++]=2*pos+1;
					if( r==M-1 )	hCGM->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu )+BLK_shift*2;
									hCGM->A_ind[A_nz++]=pos+BLK_shift*2;
					if( r<M-1 )		hCGM->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu )+BLK_shift*2;
 					hCGM->A_ptr[2*pos+2]=A_nz;
				}else if( BIT_TEST(hCGM->alg,CGM_A_BOUNDARY) )	{		//A is modified for boundary prixels
					hCGM->A_ind[A_nz++]=2*pos;
					if( c==1 )		hCGM->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu )+BLK_shift*2;
									hCGM->A_ind[A_nz++]=pos+BLK_shift*2;
					if( c<N-1 )		hCGM->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu )+BLK_shift*2;
 					hCGM->A_ptr[2*pos+1]=A_nz;
					hCGM->A_ind[A_nz++]=2*pos+1;
					if( r==1 )		hCGM->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu )+BLK_shift*2;
									hCGM->A_ind[A_nz++]=pos+BLK_shift*2;
					if( r<M-1 )		hCGM->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu )+BLK_shift*2;
 					hCGM->A_ptr[2*pos+2]=A_nz;
				}
			} 
		}
	//	| B |		(2MN:MN)		
	//	| I |		( MN:MN)		
		for( r = 0; r < M; r++ )	{		
			for( c = 0; c < N; c++ )	{
				pos = GIP_RC2POS( r,c,ldu );
				if( r>0 )		hCGM->A_ind[A_nz++]=2*GIP_RC2POS( r-1,c,ldu );
				if( r>0 )		hCGM->A_ind[A_nz++]=2*GIP_RC2POS( r-1,c,ldu )+1;
				if( c>0 )		hCGM->A_ind[A_nz++]=2*GIP_RC2POS( r,c-1,ldu );
				if( c>0 )		hCGM->A_ind[A_nz++]=2*GIP_RC2POS( r,c-1,ldu )+1;
								hCGM->A_ind[A_nz++]=2*pos;
								hCGM->A_ind[A_nz++]=2*pos+1;
				if( c==N-2 )	hCGM->A_ind[A_nz++]=2*GIP_RC2POS( r,c+1,ldu );
				if( c==N-2 )	hCGM->A_ind[A_nz++]=2*GIP_RC2POS( r,c+1,ldu )+1;
				if( r==M-2 )	hCGM->A_ind[A_nz++]=2*GIP_RC2POS( r+1,c,ldu );
				if( r==M-2 )	hCGM->A_ind[A_nz++]=2*GIP_RC2POS( r+1,c,ldu )+1;
				hCGM->A_ind[A_nz++]=pos+2*BLK_shift;
				hCGM->A_ptr[pos+2*BLK_shift+1]=A_nz;
			}
		}
		hCGM->A_nz = A_nz;			ASSERT( A_nz<=hCGM->A_dim*6 );
		CGM_update_equations( hCGM,temp );	
	}
	hCGM->g_norm=CGM_get_g( hCGM,hCGM->g,0x0 );

	hCGM->hSS=GIP_NULL;
	if( BIT_TEST( hCGM->alg,CGM_SS_DIRECT) )	{
		for( i = 0; i < 32; i++ )	setting[i]=0.0;
		if( (ret=GSS_init_id(hCGM->A_dim,hCGM->A_dim,hCGM->A_ptr,hCGM->A_ind,GIP_NULL,hCGM->m_type,setting) ) != GRUS_OK )	
		{	GIP_ERROR( "TVM_init",ret,"\tERROR at init GSS solver. ERROR CODE:%d\r\n" );		}
		hCGM->hSS = GSS_symbol_id( hCGM->A_dim,hCGM->A_dim,hCGM->A_ptr,hCGM->A_ind,hCGM->A_val );	
		if( hCGM->hSS == NULL )	
		{	GIP_ERROR( "TVM_init",-2,"\tERROR at SYMBOLIC ANALYSIS.\r\n" );		}
	}

GIP_exit:	
	GIP_free( val );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void CGM_clear( CGM_MODEL *hCGM )	{
	if( hCGM == GIP_NULL )	return;

	hCGM->M=-1;		hCGM->N=-1;
	if( hCGM->A_ptr != GIP_NULL )	
	{	GIP_free( hCGM->A_ptr );	hCGM->A_ptr=GIP_NULL; }
	if( hCGM->A_ind != GIP_NULL )	
	{	GIP_free( hCGM->A_ind );	hCGM->A_ind=GIP_NULL; }
	if( hCGM->A_val != GIP_NULL )	
	{	GIP_free( hCGM->A_val );	hCGM->A_val=GIP_NULL; }
//	if( hCGM->w != GIP_NULL )	
//	{	GIP_free( hCGM->w );	hCGM->w=GIP_NULL; }
	if( hCGM->f != GIP_NULL )	
	{	GIP_free( hCGM->f );	hCGM->f=GIP_NULL; }
	if( hCGM->g != GIP_NULL )	
	{	GIP_free( hCGM->g );	hCGM->g=GIP_NULL; }
	if( hCGM->tag != GIP_NULL )	
	{	GIP_free( hCGM->tag );	hCGM->tag=GIP_NULL; }
	if( hCGM->temp != GIP_NULL )	
	{	GIP_free( hCGM->temp );	hCGM->temp=GIP_NULL; }
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/18/2008	
*/
void CGM_solve_direct( CGM_MODEL *hCGM,double dT,double *temp )	{
	int M=hCGM->M,N=hCGM->N,MN=M*N,i,j,r,c,ldu=hCGM->ldu,*tag=hCGM->tag,ret=GIP_OK;
	int dim=hCGM->A_dim,*ptr=hCGM->A_ptr,*ind=hCGM->A_ind,isCheckRes=1;
	double *rhs=hCGM->rhs,u_max,u_min,w_max,w_min,rhs_norm=0.0,*val=hCGM->A_val,res,a;

//use GSS to solver A*u{n+1}=f(u{n})
	ret = GSS_numeric_id( hCGM->A_dim,hCGM->A_dim,hCGM->A_ptr,hCGM->A_ind,hCGM->A_val,hCGM->hSS );
	if( ret != GRUS_OK )	{
		hCGM->hSS=NULL;		//±ØÐëÉèÖÃÎªNULL,GSSÒÑ×Ô¶¯ÊÍ·ÅÄÚ´æ
		GIP_ERROR( "TVM_step",ret,"\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n" );
	}

	memcpy( temp,rhs,sizeof(double)*dim );
	if( isCheckRes==1 )		{
		for( i = 0; i < dim; i++ )	
		{	rhs_norm += (temp[i]*temp[i]);		}
		rhs_norm = sqrt(rhs_norm);	
	}
	GSS_solve_id( hCGM->hSS,hCGM->A_dim,hCGM->A_dim,hCGM->A_ptr,hCGM->A_ind,hCGM->A_val,rhs );
	if( isCheckRes==1 )	{
		double *x=rhs;
		res = 0.0;
		u_max=0.0,		u_min=DBL_MAX;			
		for( i = 0; i < dim; i++ )	{
			a = fabs(x[i]);
			if( i>= MN*2 )	{
				u_max=MAX( u_max,a );			u_min=MIN( u_min,a );
			}else	{
				w_max=MAX( u_max,a );			w_min=MIN( u_min,a );
			}
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				r = ind[j];
				temp[r] -= val[j]*x[i];
			}		
		}
		for( i = 0; i < dim; i++ )	
		{	res += (temp[i]*temp[i]);		}
		res = sqrt(res);	
		if( res/rhs_norm > 1.0e-8 )
			G_PRINTF( "\t  ****  |Ax-b|=%g,|b|=%g,|du|__=%g,|dw|__=%g\r\n",res,rhs_norm,u_max,w_max );	
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
void CGM_solve_iter( CGM_MODEL *hCGM,int flag )	{
	int M=hCGM->M,N=hCGM->N,MN=M*N,i,j,ldu=hCGM->ldu,*tag=hCGM->tag,ret=GIP_OK,t;
	int dim=hCGM->A_dim,*ptr=hCGM->A_ptr,*ind=hCGM->A_ind,isCheckRes=1,r,c,pos;
	double *rhs=hCGM->rhs,*x_0=hCGM->w,u_max,u_min,rhs_norm=0.0,*val=hCGM->A_val,a,d,*x,*temp=hCGM->temp;
	GIP_KRYLOV_SOLVER* hKS;

	x=temp+dim;
//	memcpy( x,hCGM->w,sizeof(double)*dim );
	memset( x,0x0,sizeof(double)*dim );
	if( isCheckRes==1 )		{
		for( i = 0; i < dim; i++ )	
		{	rhs_norm += (rhs[i]*rhs[i]);		}
		rhs_norm = sqrt(rhs_norm);	
	}

	hKS=GIP_Krylov_init_( dim,10,0x0,ptr,ind,val  );
	GIP_Krylov_core_( hKS,x,rhs,0x0  );

	if( isCheckRes==1 )	{
		double res = 0.0;
		u_max=DBL_MIN,		u_min=DBL_MAX;
		memcpy( temp,rhs,sizeof(double)*dim );
		for( i = 0; i < dim; i++ )	{
			a = fabs(x[i]);
			u_max=MAX( u_max,a );			u_min=MIN( u_min,a );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				r = ind[j];
				temp[r] -= val[j]*x[i];
			}		
		}
		for( i = 0; i < dim; i++ )	
		{	res += (temp[i]*temp[i]);		}
		res = sqrt(res);	
		G_PRINTF( "\t\t|Ax-b|=%g,|b|=%g,u_max=%g,u_min=%g\r\n",res,rhs_norm,u_max,u_min );	
	}

	memcpy( rhs,x,sizeof(double)*dim );

}

/*	
	energy
	注意：
		u_cur不一定等于hCGM->u

	v0.1	cys
		8/29/2008
*/
double  CGM_Energy_( CGM_MODEL *hCGM,double *u_cur,int flag )	{
	int r,c,ldu=hCGM->ldu,M=hCGM->M,N=hCGM->N,pos,pos_D;
	double res=0.0,energy=0.0,a,fx,fy,x,y,a_max;
	double *f=hCGM->f,alpha=hCGM->alpha,beta=hCGM->beta;

	a_max = 0.0;
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos_D = GIP_RC2POS( r,c,ldu );
			fx = FD_R1_0( c==N-1,hCGM->fdm,u_cur,pos_D,1 );
			fy = FD_R1_0( r==M-1,hCGM->fdm,u_cur,pos_D,ldu );
			energy += sqrt(fx*fx+fy*fy+beta)*alpha;
			a = u_cur[pos_D]-f[pos_D];
			a_max = MAX( a_max,fabs(a) );
			energy += 0.5*a*a ;
		}
	}	

	return energy;
}

/*
	[REF] A nonlinear primal-dual method for total variation-based image restoration (1999)

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/19/2008	
*/
void newton_step( CGM_MODEL *hCGM,double *temp )	{
	int i,M=hCGM->M,N=hCGM->N,MN,pos,r,c,ldu=hCGM->ldu,no;
	double *rhs=hCGM->rhs,*u=hCGM->u,*w=hCGM->w,*g=hCGM->g,sp=1.0,sd=1.0,rou=0.9,alpha=0.0,s,p_2,q;
	double en_0,en,*u_,*u_0,dot,ksi=0.5,off,off_max;
;

//	CGM_update_La( hCGM,temp );
	if( BIT_TEST( hCGM->alg,CGM_SS_DIRECT ) )
		CGM_solve_direct( hCGM,0.0,temp );
	else if( BIT_TEST( hCGM->alg,CGM_SS_ITER ) )
		CGM_solve_iter( hCGM,0x0 );

	if( BIT_TEST(hCGM->alg,CGM_METHOD_MEDIAN) )	{	
		for( i = 0; i < M*N; i++ )	{
			u[i]+=rhs[i];			
		}
		return;
	}


	MN=M*N;			ASSERT( MN*3==hCGM->A_dim );

	u_0 = GIP_alloc(sizeof(double)*3*MN );
	u_ = u_0+MN;	
	for( i=0; i < MN; i++ )	u_0[i]=u[i];
	en_0 = CGM_Energy_( hCGM,u_0,0x0 );
	sp=2.0;
	do	{
		sp /= 2;
		for( i=0; i < MN; i++ )	u_[i]=u_0[i]+rhs[2*MN+i]*sp;
		en = CGM_Energy_( hCGM,u_,0x0 );
		for( i=0; i < MN; i++ )	dot=rhs[2*MN+i]*g[i];
		en -= ksi*sp*dot;
	}while( en >= en_0 );

	for( i = 0; i < MN; i++ )	{
		pos=2*i;
		s=rhs[pos]*rhs[pos]+rhs[pos+1]*rhs[pos+1];
		p_2=(w[pos]*rhs[pos]+w[pos+1]*rhs[pos+1])/s;
		q=(w[pos]*w[pos]+w[pos+1]*w[pos+1]-1.0)/s;		
//		q=fabs(q)<FLT_EPSILON ? 0.0:q;			ASSERT( q<=0.0 );	
		q=MIN( q,0.0 );
		s=sqrt(-q+p_2*p_2)-p_2;			ASSERT( s>=0);	//s=sqrt(-q+(p*p/4))-p/2
//		if( s>0.0 )
			sd = MIN( sd,s );				//sd=sup{s:|w+sdw|<1}
	}
	sd *= rou;
	for( i = 0; i < MN; i++ )	{
		pos=2*i;	
		w[pos]+=rhs[pos]*sd;			
		w[pos+1]+=rhs[pos+1]*sd;	
		s=sqrt( w[pos]*w[pos]+w[pos+1]*w[pos+1] );		ASSERT( s <= 1.0 );
	}
	for( i = 0; i < MN; i++ )	{
		u[i]+=rhs[2*MN+i]*sp;			
	}

	CGM_update_equations( hCGM );
	hCGM->g_norm=CGM_get_g( hCGM,g,0x0 );
	off_max=0.0;
	for( r = 0; r < M; r++ )	{			//forward difference with Newmann boundary condition for dU	
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			off = fabs(g[pos]-hCGM->rhs[pos+2*MN]);
			if( off > off_max )	{
				off_max=off;		no=pos;
			}
		}
	}
	GIP_free( u_0 );
	G_PRINTF( "\t  ****  sp=%g,sd=%g,off_max=%g,no=%d\r\n",sp,sd,off_max,no );		

}

/*
	v=div(g1,g2)

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/11/2008	
*/
void CGM_get_v( CGM_MODEL *hCGM,double *v )	{
	int ldu=hCGM->ldu,M=hCGM->M,N=hCGM->N,r,c,pos,poo=-1;
	double *f=hCGM->f,*u=hCGM->u;

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			v[pos] = f[pos]-u[pos];
		}
	}	;
}

/*
	A testing code 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/5/2008	
*/
void CGM_MODEL_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName )	{
	double *temp,res,tol;
	CGM_MODEL cgm,*hCGM=&cgm;
	int M,N,isTrace=1,nKey,t_step=1;
	char sTitle[80];

	CGM_init_( hIMG,hIMG_mask,hCGM );
	M=hCGM->M,		N=hCGM->N;
	temp=hCGM->temp;
	res=hCGM->nn_res;
	tol = res*1.0e-4;
	G_PRINTF( "START:\tM=%d,N=%d,alpha=%g,beta=%g,alg=%x,res=%g\r\n",M,N,hCGM->alpha,hCGM->beta,hCGM->alg,res );		

//	res = CGM_update_res( hCGM,temp );
//	GIP_show_IplImage_d( M,N,sWndName,hCGM->z,hCGM->ldu,0x0 );
	while( res>tol || BIT_TEST(hCGM->flag,GIP_SOLVER_CONVERGE) )	{
		if( isTrace==1 && cgm_step%t_step==0 )	{
			_stprintf( sTitle,".\\trace\\u_%d.bmp",cgm_step );		
			GIP_save_IplImage_d( M,N,sTitle,(hCGM->u+hCGM->shift),hCGM->ldu,0x0 );
			CGM_get_v( hCGM,temp );
			_stprintf( sTitle,".\\trace\\v_%d.bmp",cgm_step );		
			GIP_save_IplImage_d( M,N,sTitle,temp,hCGM->ldu,0x0 );
//			GIP_show_IplImage_d( M,N,sWndName,temp,hCGM->ldu,0x0 );
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
//			GIP_save_IplImage_d( M,N,sTitle,(temp+hCGM->shift),hCGM->ldu,0x0 );
		}
		newton_step( hCGM,temp );
		res = hCGM->nn_res;
		G_PRINTF( "CGM_%d: energy=%g,res=%g,|g|=%g\r\n",cgm_step,
			CGM_Energy_(hCGM,hCGM->u,0x0),hCGM->nn_res,hCGM->g_norm	);

//		res = CGM_update_res( hCGM,temp );
		if( BIT_TEST(hCGM->flag,GIP_SOLVER_CONVERGE) )
			break;
		cgm_step++;
	}
	G_PRINTF( "GGM iteration finish\r\n"	);
GIP_exit:
	nKey = cvWaitKey( );
	CGM_clear( hCGM );
}
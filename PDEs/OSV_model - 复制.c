/*
	[REF] Image decomposition and restoration using total variation minimization and the H-1 norm
	by Stanley Osher, Andrs Sol, Luminita Vese

	Fixed point iteration with implicit format
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
#include "OSV_model.h"
#include "GIP_fmm.h"
#include "../../lib/gss_6_demo.h"

static int osv_step=0;

/*
	[REF] Modeling Textures with Total Variation Minimization and Oscillating Patterns in Image Processing

	注意：
		1、hOSV->ldu=N

	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/5/2008	
*/
void OSV_init_( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,OSV_MODEL *hOSV )	{
	int M,N,MN,r,c,pos,A_nz,i,setting[32],ret=GIP_OK,ldu,BLK_shift;
	GIP_FLOATING *val=GIP_NULL;
	double *g,a,s,*f,fx,fy;

	osv_step=0;

	memset( hOSV,0x0,sizeof(OSV_MODEL) );
	GIP_import_IplImage( hIMG,&M,&N,&val,0x0 );

	hOSV->alg = OSV_SS_DIRECT;			//OSV_SS_ITER;			//
	hOSV->lenda=1.0;
	hOSV->mui =0.01;			hOSV->beta = 1.0e-12;
	hOSV->M=M;				hOSV->N=N;
	MN=M*N;
	hOSV->ldu=N;			hOSV->shift=0;
	hOSV->f=(double*)GIP_alloc( sizeof(double)*(MN) );
	hOSV->u=(double*)GIP_alloc( sizeof(double)*MN*2 );
	hOSV->w=hOSV->u+MN;		
	hOSV->tag=(int*)GIP_alloc( sizeof(int)*(M+2)*(N+2) );
	memset( hOSV->tag,0x0,sizeof(int)*(M+2)*(N+2) );
	ldu=hOSV->ldu;			
	s = 1.0/(2*hOSV->lenda);
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			hOSV->f[pos] = val[pos];
			hOSV->u[pos] = val[pos];
			if( c==0 )		BIT_SET( hOSV->tag[pos],GIP_BD_LEFT );
			if( c==N-1 )	BIT_SET( hOSV->tag[pos],GIP_BD_RIGHT );
			if( r==0 )		BIT_SET( hOSV->tag[pos],GIP_BD_DOWN );
			if( r==M-1 )	BIT_SET( hOSV->tag[pos],GIP_BD_UP );
			if( c==0 || c==N-1 || r==0 || r==M-1 )	{
				hOSV->w[pos] = 0.0;		
			}else	{
				fx=(val[GIP_RC2POS(r,c+1,ldu)]-val[GIP_RC2POS(r,c-1,ldu)])/2.0;
				fy=(val[GIP_RC2POS(r+1,c,ldu)]-val[GIP_RC2POS(r-1,c,ldu)])/2.0;
				hOSV->w[pos] = fx+fy;				
			}			
		}
	}

	hOSV->A_dim=MN*2;			
	hOSV->m_type=0;		//unsym
	A_nz=hOSV->m_type==11 || hOSV->m_type==12 ? hOSV->A_dim*6 : hOSV->A_dim*6;			//取上限
	hOSV->rhs=(int*)GIP_alloc( sizeof(double)*(hOSV->A_dim) );
	hOSV->A_ptr=(int*)GIP_alloc( sizeof(int)*(hOSV->A_dim+1) );
	hOSV->A_ind=(int*)GIP_alloc( sizeof(int)*(A_nz) );
	hOSV->A_val=(double*)GIP_alloc( sizeof(double)*(A_nz) );
	hOSV->A_ptr[0]=0;			A_nz=0;
	BLK_shift=M*N;
	ASSERT( hOSV->m_type!=11 && hOSV->m_type!=12 );
//	U Block		
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
		//[U:U]
			hOSV->A_ind[A_nz++]=pos;
		//[W:U]
			if( r>0 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu )+BLK_shift;
			if( c>0 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu )+BLK_shift;
			hOSV->A_ind[A_nz++]=pos+BLK_shift;
			if( c<N-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu )+BLK_shift;
			if( r<M-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu )+BLK_shift;
 			hOSV->A_ptr[pos+1]=A_nz;
		}
	}
//	W Block		
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			if( r>0 )	hOSV->A_ind[A_nz++]=GIP_RC2POS( r-1,c,ldu );
			if( c>0 )	hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c-1,ldu );
			hOSV->A_ind[A_nz++]=pos;
			if( c<N-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r,c+1,ldu );
			if( r<M-1 )		hOSV->A_ind[A_nz++]=GIP_RC2POS( r+1,c,ldu );
			hOSV->A_ind[A_nz++]=pos+BLK_shift;
			hOSV->A_ptr[pos+BLK_shift+1]=A_nz;
		}
	}
	hOSV->A_nz = A_nz;			ASSERT( A_nz<=hOSV->A_dim*9 );

	f=hOSV->f;
	memset( hOSV->rhs,0x0,sizeof(double)*hOSV->A_dim );
	s=2.0*hOSV->lenda;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
		hOSV->rhs[pos]=s*f[pos];
		hOSV->rhs[pos+MN]=0.0;
	}
	}

	hOSV->hSS=GIP_NULL;
	if( BIT_TEST( hOSV->alg,OSV_SS_DIRECT) )	{
		for( i = 0; i < 32; i++ )	setting[i]=0.0;
		if( (ret=GSS_init_id(hOSV->A_dim,hOSV->A_dim,hOSV->A_ptr,hOSV->A_ind,GIP_NULL,hOSV->m_type,setting) ) != GRUS_OK )	
		{	GIP_ERROR( "TVM_init",ret,"\tERROR at init GSS solver. ERROR CODE:%d\r\n" );		}
		hOSV->hSS = GSS_symbol_id( hOSV->A_dim,hOSV->A_dim,hOSV->A_ptr,hOSV->A_ind,hOSV->A_val );	
		if( hOSV->hSS == NULL )	
		{	GIP_ERROR( "TVM_init",-2,"\tERROR at SYMBOLIC ANALYSIS.\r\n" );		}
	}

/*
	OSV_Init_u0( hOSV,hIMG_mask,0x0 );
	GIP_FMM_smooth_u_( hOSV->M,hOSV->N,hOSV->u,hOSV->ldu,hOSV->tag,0x0 );
*/
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
	if( hOSV->u != GIP_NULL )	
	{	GIP_free( hOSV->u );	hOSV->u=GIP_NULL; }
	if( hOSV->f != GIP_NULL )	
	{	GIP_free( hOSV->f );	hOSV->f=GIP_NULL; }
	if( hOSV->tag != GIP_NULL )	
	{	GIP_free( hOSV->tag );	hOSV->tag=GIP_NULL; }
}

/*

	v0.1	cys
		5/30/2008	
*/
void OSV_Du( OSV_MODEL *hOSV,double *Va,double *Du )	{
	int M=hOSV->M,N=hOSV->N,ldu=hOSV->ldu,r,c,pos;
	double ux,uy,g,beta=hOSV->beta;

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
void OSV_verify_UU( OSV_MODEL *hOSV,double *temp )	{
	int dim=hOSV->A_dim,*ptr=hOSV->A_ptr,*ind=hOSV->A_ind,isCheckRes=1,M=hOSV->M,N=hOSV->N,MN;
	int isFind,i,j,k,r;
	double *val=hOSV->A_val,a,sum;

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
void OSV_update_La( OSV_MODEL *hOSV,double *Va,double *temp )	{
	int dim=hOSV->A_dim,*ptr=hOSV->A_ptr,ldu=hOSV->ldu,M=hOSV->M,N=hOSV->N;
	int i,r,c,pos,r_D,c_D,A_nz,pos_D,alg=hOSV->alg,BLK_shift,nz_D;
	double *val=hOSV->A_val,*u,*w=hOSV->w,mui=hOSV->mui,s,a,b,a_max,a_min,lenda=hOSV->lenda,x,y,beta=hOSV->beta;
	
//	OSV_Du( hOSV,Va,temp );

	BLK_shift=M*N;
	u = hOSV->u;
	A_nz=0;
	s = 1.0;
//	U column		
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos_D=GIP_RC2POS( r,c,ldu );
//[U:U]
			val[A_nz++]=2*lenda;
//[W:U]
			a=0.0;
			if( r>0 )	{
				pos = GIP_RC2POS( r-1,c,ldu );
				x=c>0?u[pos]-u[pos-1]:0;		y=u[pos_D]-u[pos];
				b = 1.0/sqrt(x*x+y*y+beta);
				a += s*b;			val[A_nz++]=-s*b;
			}
			if( c>0 )	{
				pos = GIP_RC2POS( r,c-1,ldu );
				x=u[pos_D]-u[pos];				y=r>0?u[pos]-u[pos-ldu]:0;		
				b = 1.0/sqrt(x*x+y*y+beta);
				a += s*b;			val[A_nz++]=-s*b;
			}
			nz_D=A_nz++;
			if( c<N-1 )		{
				pos = GIP_RC2POS( r,c+1,ldu );
				x=u[pos_D]-u[pos];				y=r>0?u[pos_D]-u[pos_D-ldu]:0;
				b = 1.0/sqrt(x*x+y*y+beta);
				a += s*b;			val[A_nz++]=-s*b;
			}
			if( r<M-1 )		{
				pos = GIP_RC2POS( r+1,c,ldu );
				x=c>0?u[pos_D]-u[pos_D-1]:0;	y=u[pos_D]-u[pos];
				b = 1.0/sqrt(x*x+y*y+beta);
				a += s*b;			val[A_nz++]=-s*b;
			}
			val[nz_D]=a;
			ASSERT( A_nz==hOSV->A_ptr[pos_D+1] );
		}
	}
//	W column		
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos_D=GIP_RC2POS( r,c,ldu )+BLK_shift;
//		U:W
			if( r>0 )		val[A_nz++]=1.0;
			if( c>0 )		val[A_nz++]=1.0;
			val[A_nz++]=-4.0;
			if( c<N-1 )		val[A_nz++]=1.0;
			if( r<M-1 )		val[A_nz++]=1.0;
//		W:W
			val[A_nz++]=1.0;
			
			ASSERT( A_nz==hOSV->A_ptr[pos_D+1] );
		}
	}

	a_max=DBL_MIN;		a_min=DBL_MAX;
	for( i = 0; i < A_nz; i++ )		{
		a = fabs(val[i]); 
		a_max=MAX(a_max,a);		a_min=MIN(a_min,a);	
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/8/2008	
*/
void OSV_step( OSV_MODEL *hOSV,double dT,double *temp )	{
	int M=hOSV->M,N=hOSV->N,MN=M*N,i,j,r,c,g_pos,pos,ldu=hOSV->ldu,*tag=hOSV->tag,ret=GIP_OK;
	int dim=hOSV->A_dim,*ptr=hOSV->A_ptr,*ind=hOSV->A_ind,isCheckRes=1;
	double *rhs=hOSV->rhs,*u=hOSV->u,u_max,u_min,rhs_norm=0.0,*val=hOSV->A_val,res,a;

	OSV_update_La( hOSV,u,temp );
//use GSS to solver A*u{n+1}=f(u{n})
//	if( osv_step%5==0 )	{
		ret = GSS_numeric_id( hOSV->A_dim,hOSV->A_dim,hOSV->A_ptr,hOSV->A_ind,hOSV->A_val,hOSV->hSS );
		if( ret != GRUS_OK )	{
			hOSV->hSS=NULL;		//必须设置为NULL,GSS已自动释放内存
			GIP_ERROR( "TVM_step",ret,"\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n" );
		}
//	}

//	OSV_update_rhs( hOSV,rhs,temp+M*N*3 );
	memcpy( temp,rhs,sizeof(double)*dim );
	if( isCheckRes==1 )		{
		for( i = 0; i < dim; i++ )	
		{	rhs_norm += (temp[i]*temp[i]);		}
		rhs_norm = sqrt(rhs_norm);	
	}
	GSS_solve_id( hOSV->hSS,hOSV->A_dim,hOSV->A_dim,hOSV->A_ptr,hOSV->A_ind,hOSV->A_val,temp );
	memcpy( hOSV->u,temp,sizeof(double)*dim );
	if( isCheckRes==1 )	{
		res = 0.0;
		u_max=DBL_MIN,		u_min=DBL_MAX;
		memcpy( temp,rhs,sizeof(double)*dim );
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
		printf( "\t\t|Ax-b|=%g,|b|=%g,u_max=%g,u_min=%g\r\n",res,rhs_norm,u_max,u_min );	
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
void OSV_step_iter( OSV_MODEL *hOSV,double *temp )	{
	int M=hOSV->M,N=hOSV->N,MN=M*N,i,j,ldu=hOSV->ldu,*tag=hOSV->tag,ret=GIP_OK,t;
	int dim=hOSV->A_dim,*ptr=hOSV->A_ptr,*ind=hOSV->A_ind,isCheckRes=1,r,c,pos;
	double *rhs=hOSV->rhs,*u=hOSV->u,u_max,u_min,rhs_norm=0.0,*val=hOSV->A_val,a,d,*x;

	OSV_update_La( hOSV,u,temp );
	
	x=temp+dim;
	memcpy( x,rhs,sizeof(double)*dim );

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
	memcpy( u,x,sizeof(double)*dim );

	if( isCheckRes==1 )	{
		double res = 0.0;
		u_max=DBL_MIN,		u_min=DBL_MAX;
		memcpy( temp,rhs,sizeof(double)*dim );
		for( i = 0; i < dim; i++ )	{
			a = fabs(rhs[i]);
			u_max=MAX( u_max,a );			u_min=MIN( u_min,a );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				r = ind[j];
				temp[r] -= val[j]*u[i];
			}		
		}
		for( i = 0; i < dim; i++ )	
		{	res += (temp[i]*temp[i]);		}
		res = sqrt(res);	
		printf( "\t\t|Ax-b|=%g,u_max=%g,u_min=%g\r\n",res,u_max,u_min );	
	}

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/8/2008	
*/
double OSV_energy_( OSV_MODEL *hOSV,double *temp )	{
	int ldu=hOSV->ldu,M=hOSV->M,N=hOSV->N,r,c,alg=hOSV->alg,pos,poo=-1;
	double energy=0.0,beta=hOSV->beta,mui=hOSV->mui,lenda=hOSV->lenda,a,a_max,a_min,ux,uy,g1x,g2y;
	double *u=hOSV->u,*f=hOSV->f,*w=hOSV->w,u_oo=0,g1_oo=0,g2_oo=0;
/*
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

	printf( "OSV_%d: energy=%g,|f-u|=%g,|g1|=%g(%d),|g2|=%g\r\n",osv_step,energy,u_oo,g1_oo,poo,g2_oo	);*/
	printf( "OSV_%d: \r\n",osv_step	);
	return energy;
}

/*	
	off=Vb-L(Va)*Va
	Normalized residual.[REF]:	ACCELERATION METHODS FOR TOTAL VARIATION-BASED IMAGE DENOISING

	v0.1	cys
		5/23/2008

double  OSV_NormalRes( OSV_MODEL *hOSV,double *temp )	{
	int i,j,M=hTVM->M,N=hTVM->N,dim=hTVM->A_dim,*ptr=hTVM->A_ptr,*ind=hTVM->A_ind;
	double off=0.0,a,*val=hTVM->A_val,d;

	OSV_update_La( hTVM,temp );

	memcpy( temp,hOSV->u,MN_size );
	memcpy( temp+MN,hOSV->g1,MN_size );
	memcpy( temp+MN*2,hOSV->g2,MN_size );
	rhs=temp+dim;
	memset( rhs,0x0,sizeof(double)*dim );
	memcpy( rhs,hOSV->z,MN_size );

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
void OSV_get_v( OSV_MODEL *hOSV,double *v )	{
	int ldu=hOSV->ldu,M=hOSV->M,N=hOSV->N,r,c,pos,poo=-1;
	double g1x,g2y,*f=hOSV->f,*u=hOSV->u;

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
void OSV_MODEL_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName )	{
	double *temp,energy;
	OSV_MODEL meyer,*hOSV=&meyer;
	int M,N,t_len,isTrace=1,nKey,t_step=1;
	char sTitle[80];

	OSV_init_( hIMG,hIMG_mask,hOSV );
	printf( "START:\tmui=%g,lenda=%g,alg=%x\r\n",hOSV->mui,hOSV->lenda,hOSV->alg );		

	M=meyer.M,		N=meyer.N;
	t_len=MAX((M+2)*(N+2)*6,N);
	temp=GIP_alloc( sizeof(double)*t_len );
//	tol = TVM_NormalRes( hOSV,temp );
//	tol *= 1.0e-4;
//	GIP_show_IplImage_d( M,N,sWndName,meyer.z,meyer.ldu,0x0 );
	while( !BIT_TEST(meyer.flag,GIP_SOLVER_CONVERGE) )	{
		if( isTrace==1 && osv_step%t_step==0 )	{
			sprintf( sTitle,".\\trace\\u_%d.bmp",osv_step );		
			GIP_save_IplImage_d( M,N,sTitle,(meyer.u+meyer.shift),meyer.ldu,0x0 );
			OSV_get_v( hOSV,temp );
			sprintf( sTitle,".\\trace\\v_%d.bmp",osv_step );		
			GIP_save_IplImage_d( M,N,sTitle,temp,meyer.ldu,0x0 );
//			GIP_show_IplImage_d( M,N,sWndName,temp,meyer.ldu,0x0 );
			nKey = cvWaitKey( );
			switch( nKey )	{
				case 27:
					printf( "IMAGE decomposition iteration break out\r\n"	);
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
 		if( BIT_TEST( hOSV->alg,OSV_SS_DIRECT ) )
			OSV_step( &meyer,0.0,temp );
		else if( BIT_TEST( hOSV->alg,OSV_SS_ITER ) )
			OSV_step_iter( &meyer,temp );
		energy=OSV_energy_( &meyer,temp ); 
		if( BIT_TEST(meyer.flag,GIP_SOLVER_CONVERGE) )
			break;
		osv_step++;
		//if( average > 2 )	{
		//	OSV_ReInit_u0( &meyer,0x0 );
		//	OSV_smooth_u_( &meyer );
		//}
	}
	printf( "GAC iteration finish\r\n"	);
	GIP_free( temp );
GIP_exit:
	OSV_clear( &meyer );
}
#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "../PCH/GIP_def.h"
#include "GIP_ss.h"
#include "gip_util.h"
#include "Compact_DataStruc.h"
#include "../../lib/gss_6_demo.h"

static char BLAS_N[1]={'N'},BLAS_T[1]={'T'},BLAS_C[1]={'C'},BLAS_L[1]={'L'},BLAS_R[1]={'R'},BLAS_U[1]={'U'};
static int inc_1=1; 
static GIP_FLOATING one_=1.0,fuyi_=-1.0,zero_=0.0;		 

/*
	b=Ax, A 缺省采用压缩列格式存储

	v0.1	cys
		8/25/2008
*/
void GIP_AV( int dim,int *ptr,int *ind,GIP_FLOATING *val,GIP_FLOATING *x,GIP_FLOATING *b,int flag )	{
	int i,k,r,c;
	double s;
	for( i = 0; i < dim; i++ )	{
		b[i] = 0.0;
	}
	if( BIT_TEST(flag,GIP_ROW_MAJOR) )	{		//压缩行格式存储
		for( i = 0; i < dim; i++ )	{
			for( k = ptr[i];k < ptr[i+1]; k++ )	{
				c = ind[k];
				b[i] += val[k]*x[c];
			}
		};
	}else	{									//缺省采用压缩列格式存储
		for( i = 0; i < dim; i++ )	{
			s = x[i];
			for( k = ptr[i];k < ptr[i+1]; k++ )	{
				r = ind[k];
				b[r] += val[k]*s;
			}
		};
	}
}

#define AMG_CF_C	0x1
#define AMG_CF_F	0x2
#define AMG_CF_U	0x4
/*
	[REF]	An Accelerated Algebraic Multigrid Algorithm for Total-Variation Denoising

	C/F-Splitting Algorithm
	temp[5*dim]

	v0.1	cys
		8/24/2008
*/
void AMG_CF_Split( GIP_SPARSE_SOLVER* hSS,int flag )	{
	int i,cur,r,j,l,dim=hSS->dim,*ptr=hSS->ptr,*ind=hSS->ind,*S_ptr=0x0,*S_ind=0x0,*ST_ptr=0x0,*ST_ind=0x0,nnz;
	int *lenda=hSS->I_temp,*split=lenda+dim,*connect=split+dim,pact=0;
	PACT_order_list POL;

	for( i = 0; i < dim; i++ )	{
		connect[i]=0;
		lenda[i]=0;
		pact=MAX( pact,ptr[i+1]-ptr[i] );
		split[i]=AMG_CF_U;
	}
	PACT_ol_init( &POL,pact,dim,0x0 );
	lenda=POL.list;
	nnz=S_ind[dim];
	S_ptr=(int*)GIP_alloc( sizeof(int)*(dim+1) );
	S_ind=(int*)GIP_alloc( sizeof(int)*(nnz) );
	for( i = 0; i < dim; i++ )	{
		lenda[i]=S_ptr[i+1]-S_ptr[i];
		PACT_ol_insert( &POL,i,lenda[i] );		
	}

	while( (cur=PACT_ol_maxpop(&POL) ) != -1 )	{
		split[cur]=AMG_CF_C;
		for( i = ST_ptr[cur]; i < ST_ptr[cur+1]; i++ )	{
			r = ST_ind[i];
			if( split[r]!=AMG_CF_U )	continue;
			split[r]=AMG_CF_F;
			for( j=S_ptr[r]; j < S_ptr[r+1]; j++ )	{
				l = S_ind[j];
				if( split[l]!=AMG_CF_U )	continue;
				PACT_ol_update( &POL,l,lenda[l]+1 );		
			}
		}
		for( j=S_ptr[cur]; j < S_ptr[cur+1]; j++ )	{
			l = S_ind[j];
			if( split[l]!=AMG_CF_U )	continue;
			PACT_ol_update( &POL,l,lenda[l]-1 );		
		}
	}

}

/*
	z=A*(v1,v2,...vj)*H(:,j)		其中z=V(:,j+1)=A*vj
	1、orthogonalization
	2、update H

	v0.1	cys
		8/25/2007
*/
void GIP_Krylov_update_( GIP_KRYLOV_SOLVER* hSolver,int j,int flag  )	{
	GIP_FLOATING *h,*z,*c,*u,s;
	int dim=hSolver->dim,ldh=hSolver->ldh,n=j+1,i,t;
	double norm,norm_0,norm_1;
	clock_t start=clock( );

	h=hSolver->H+j*ldh;
	z=hSolver->V+(j+1)*dim;
	norm_0 = BLAS_NRM2_0( dim,z,inc_1 );
//	orthogonal:		z = z-V*(V,z)
	BLAS_GEMV_0( BLAS_T,dim,n,one_,hSolver->V,dim,z,inc_1,zero_,h,inc_1 );	//h=V'z
	BLAS_GEMV_0( BLAS_N,dim,n,fuyi_,hSolver->V,dim,h,inc_1,one_,z,inc_1 );	
//DGKS correction if needed
	norm = BLAS_NRM2_0( dim,z,inc_1 );	//normalize z
/*	if( norm < hSolver->ortho_thresh*norm_0 )	{
		clock_t s_0=clock( );
		norm_0 = norm;
		c = hSolver->work;
		BLAS_GEMV( BLAS_C,dim,n,one_,hSolver->V,dim,z,inc_1,zero_,c,inc_1 );	//h=-V'z
		BLAS_GEMV( BLAS_N,dim,n,fuyi_,hSolver->V,dim,c,inc_1,one_,z,inc_1 );	
		i=j+1;			BLAS_AXPY( i,one_,c,inc_1,h,inc_1 );	
		norm = BLAS_NRM2( dim,z,inc_1 );	//normalize z
		if( norm < 0.707*norm_0 )	{	//z lies in the span of V_{j} numerically
			G_PRINTF("\r\n!!!z lies in the span of V_{j} numerically!!!\r\n");
//			for( i = 0; i < dim; i++ )	CLEAR( z[i] );
		}
		GRUS_SE_INFO[GRUS_SE_REORTHO_TIME] += clock( )-s_0;	
	}*/
	h[j+1]=norm;
	s=1.0/norm;
	BLAS_SCAL_0( dim,s,z,inc_1 );
}

/*
	z=A*(v1,v2,...vj)*H(:,j)		其中z=V(:,j+1)=A*vj
	1、orthogonalization
	2、update H

	v0.1	cys
		9/6/2007
*/
void GIP_Krylov_expand_( GIP_KRYLOV_SOLVER* hKS,int p,int m  )	{
	int i,j,dim=hKS->dim,ldh=hKS->ldh;
	double norm;
	GIP_FLOATING s,*V=hKS->V;

//	hKS->V_m = m;
	for( j =p; j < m; j++ )	{
/*		if( hSolver->balance_mode != GSE_BALANCE_NO )
		{	for( i = 0; i< dim; i++ )	SCALE_RECIP( V[j*dim+i],D[i] );	}
		GRUS_driver_switch( V+j*dim,V+(j+1)*dim,GRUS_DRIVER_OPX,0 );
		if( hSolver->balance_mode != GSE_BALANCE_NO )	{
			for( i = 0; i < dim; i++ )	SCALE_DIV( V[j*dim+i],D[i] );
			for( i = 0; i < dim; i++ )	SCALE_DIV( V[(j+1)*dim+i],D[i] );
		}
		if( hSolver->AV != GRUS_NULL )	
			BLAS_COPY( dim,V+(j+1)*dim,inc_1,hSolver->AV+j*dim,inc_1 );*/
//		GIP_AV( dim,hKS->ptr,hKS->ind,hKS->val,V+j*dim,V+(j+1)*dim );
		GIP_Krylov_update_( hKS,j,0x0 );		
		norm=fabs(hKS->H[j*ldh+j+1]);	
		if(  norm < hKS->tole )
			break;
	}
}

/*
	注意：
	1 	KRYLOV SPARSE SOLVER原则上采用矩阵的隐式表达，暂保留(ptr,ind,val）结构

	v0.1	cys
		8/25/2007
*/
GIP_KRYLOV_SOLVER* GIP_Krylov_init_( int dim,int m,int flag,int *ptr,int*ind,double *val  )	{
	int dimQH;
	GIP_KRYLOV_SOLVER* hKS=GIP_alloc(sizeof(GIP_KRYLOV_SOLVER) );

	memset( hKS,0x0,sizeof(GIP_KRYLOV_SOLVER) );
	hKS->dim=dim;	
	hKS->ptr=ptr;	hKS->ind=ind;	hKS->val=val;	
	hKS->V_m=m;
	hKS->MAX_RESTART=1;

	hKS->tole = 0.1;

	dimQH=hKS->ldh=hKS->V_m+1;
	hKS->d_temp=GIP_alloc( sizeof(double)*MAX(dim,dimQH*3) );	//see LAPACK_HSEIN,LAPACK_HEEVX
	hKS->i_temp=GIP_alloc( sizeof(double)*MAX(dim,dimQH*2) );
	hKS->V=GIP_alloc( sizeof(GIP_FLOATING)*dim*dimQH );	
	hKS->H=GIP_alloc( sizeof(GIP_FLOATING)*dimQH*dimQH );
	memset( hKS->H,0x0,sizeof(GIP_FLOATING)*dimQH*dimQH );
	hKS->nRestart = 0;
/*	
	hKS->a_norm = 0.0;
	mch_eps = dlamch( cmach );
	if( hKS->stop_mode==GSE_STOP_DIRECT )
		hKS->tole = 100*MAX(param[GSE_CONVERGE_TOL],mch_eps);		//The convergence tolerance
	else
		hKS->tole = param[GSE_CONVERGE_TOL];
	hKS->ortho_thresh = param[GSE_REORTHO_THRESH];	//sin(param[GSE_REORTHO_THRESH]/180.0*PI);	
*/
	return hKS;
}

/*
	copy from http://www.netlib.org/templates/cpp//gmres.h
	v0.1	cys
		8/26/2007
*/
void GeneratePlaneRotation( GIP_FLOATING dx,GIP_FLOATING dy,GIP_FLOATING *cs,GIP_FLOATING *sn )		{
  if (dy == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (fabs(dy) > fabs(dx)) {
    GIP_FLOATING temp = dx / dy;
    *sn = 1.0 / sqrt( 1.0 + temp*temp );
    *cs = temp * (*sn);
  } else {
    GIP_FLOATING temp = dy / dx;
    *cs = 1.0 / sqrt( 1.0 + temp*temp );
    *sn = temp * (*cs);
  }
}

/*
	copy from http://www.netlib.org/templates/cpp//gmres.h
	v0.1	cys
		8/26/2007
*/
void ApplyPlaneRotation(GIP_FLOATING *dx, GIP_FLOATING *dy, GIP_FLOATING cs, GIP_FLOATING sn){
	GIP_FLOATING temp  =  cs * (*dx) + sn * (*dy);
	*dy = -sn * (*dx) + cs * (*dy);
	*dx = temp;
}

/*
	copy from http://www.netlib.org/templates/cpp//gmres.h
	v0.1	cys
		8/26/2007
*/
void GIP_GMRES_Update( int m, GIP_KRYLOV_SOLVER* hKS, GIP_FLOATING *s, GIP_FLOATING *x )	{
	int dim=hKS->dim,k,j,ldh=hKS->ldh;
	GIP_FLOATING *V=hKS->V,*vH;

	for ( k = m-1; k >= 0; k--) {
		vH = hKS->H+ldh*k;
		s[k] /= vH[k];
		for ( j = k - 1; j >= 0; j--)
		  s[j] -= vH[j] * s[k];
	}
	//			for (int j = 0; j <= k; j++)
	//				x += v[j] * y(j);
	BLAS_GEMV_0( BLAS_N,dim,m,one_,V,dim,s,inc_1,one_,x,inc_1 );	//x+=Vy
}

/*
	v0.1	cys
		8/25/2007
*/
void GIP_Krylov_core_( GIP_KRYLOV_SOLVER* hKS,GIP_FLOATING *x_0,GIP_FLOATING *rhs,int flag  )	{
	int m=hKS->V_m,dim=hKS->dim,i,k,ldh=hKS->ldh,isConverge=0;
	GIP_FLOATING *V,*vH,normb,beta,*s,*cs,*sn,a;
	double res;

	normb = BLAS_NRM2_0( dim,rhs,inc_1 );
	V=hKS->V;
	GIP_AV( dim,hKS->ptr,hKS->ind,hKS->val,x_0,V,flag );
	for( i = 0; i < dim; i++ )	V[i]=rhs[i]-V[i];
	beta = BLAS_NRM2_0( dim,V,inc_1 );
	a=1.0/beta;
	BLAS_SCAL_0( dim,a,V,inc_1 );

	cs=hKS->d_temp;		sn=cs+m+1;		s=sn+m+1;
	for( i=0; i < m+1; i++ )	s[i]=0.0;
	s[0]=beta;
	for( i=0; i < m; i++ )	{
		GIP_AV( dim,hKS->ptr,hKS->ind,hKS->val,V+i*dim,V+(i+1)*dim,flag );
		GIP_Krylov_update_( hKS,i,0x0 );

		vH = hKS->H+ldh*i;
		for (k = 0; k < i; k++)	{
//			ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));  
			ApplyPlaneRotation( vH+k, vH+k+1, cs[k], sn[k] );      
		}
//		GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
		GeneratePlaneRotation( vH[i],vH[i+1], cs+i, sn+i);
		ApplyPlaneRotation( vH+i, vH+i+1, cs[i], sn[i]);
		ApplyPlaneRotation( s+i, s+i+1, cs[i], sn[i]);
		res = fabs(s[i+1]);

		if ( res < hKS->tole*normb ) {
			GIP_GMRES_Update( i+1, hKS,s,x_0 );	
			isConverge=1;
			break;
		}
	}
	if( isConverge==0 )
		GIP_GMRES_Update( m, hKS,s,x_0 );

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		9/27/2008		
*/
void *GIP_SS_init_( int dim,int nz_max,int *_ptr,int *_ind,double *_val,int m_type,int flag )		{
	int i,ret;
	GIP_SPARSE_SOLVER *hSS=GIP_NULL;

	hSS=GIP_alloc( sizeof(GIP_SPARSE_SOLVER) );
	memset( hSS,0x0,sizeof(GIP_SPARSE_SOLVER) );

	hSS->alg=GIP_SS_DIRECT;

	hSS->dim=dim;		
	hSS->nz=nz_max;
	hSS->m_type=m_type;
	hSS->ptr = GIP_alloc( sizeof(int)*(dim+1) );
	hSS->ind = GIP_alloc( sizeof(int)*nz_max );
	hSS->val = GIP_alloc( sizeof(double)*nz_max );
	hSS->rhs = GIP_alloc( sizeof(double)*dim );
	hSS->temp = GIP_alloc( sizeof(double)*dim );
	memcpy( hSS->ptr,_ptr,sizeof(int)*(dim+1) );
	memcpy( hSS->ind,_ind,sizeof(int)*nz_max );

	if( BIT_TEST( hSS->alg,GIP_SS_DIRECT) )	{
		double setting[32];
		for( i = 0; i < 32; i++ )	setting[i]=0.0;
	//	setting[5] = 65791;
		if( (ret=GSS_init_id(hSS->dim,hSS->dim,hSS->ptr,hSS->ind,GIP_NULL,hSS->m_type,setting) ) != GRUS_OK )	
		{	GIP_ERROR( "TVM_init",ret,"\tERROR at init GSS solver. ERROR CODE:%d\r\n" );		}
		hSS->hLU = GSS_symbol_id( hSS->dim,hSS->dim,hSS->ptr,hSS->ind,hSS->val );	
		if( hSS->hLU == NULL )	
		{	GIP_ERROR( "TVM_init",-2,"\tERROR at SYMBOLIC ANALYSIS.\r\n" );		}
	}
GIP_exit:
	return hSS;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/6/2008		
*/
void GIP_SS_solve_( GIP_SPARSE_SOLVER *hSS,int flag )		{
	int i,j,ret,isCheckRes=1,r;
	int dim=hSS->dim,*ptr=hSS->ptr,*ind=hSS->ind;
	double *val=hSS->val,*rhs=hSS->rhs,*temp=hSS->temp;
	double rhs_norm=0.0,res,u_max=0.0,u_min,a;

	if( 0 )	{
		char sPath[80];
		_stprintf( sPath,"H:\\GIPack\\data\\SS_%d_0.mtx",dim );
		Equation_output( 0,0,dim,ptr,ind,val,rhs,sPath );
	}

	if( BIT_TEST( hSS->alg,GIP_SS_DIRECT) )	{
		ret = GSS_numeric_id( hSS->dim,hSS->dim,hSS->ptr,hSS->ind,hSS->val,hSS->hLU );
		if( ret != GRUS_OK )	{
			hSS->hLU=NULL;		
			GIP_ERROR( "TVM_step",ret,"\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n" );
		}

		memcpy( temp,rhs,sizeof(double)*dim );
		if( isCheckRes==1 )		{
			for( i = 0; i < dim; i++ )	
			{	rhs_norm += (temp[i]*temp[i]);		}
			rhs_norm = sqrt(rhs_norm);	
		}
		GSS_solve_id( hSS->hLU,hSS->dim,hSS->dim,hSS->ptr,hSS->ind,hSS->val,rhs );
	}
	if( isCheckRes==1 )	{
		double *x=rhs;
		res = 0.0;
		u_max=0.0,		u_min=DBL_MAX;			
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
		if( res/rhs_norm > 1.0e-8 )
			G_PRINTF( "\r\n\t  ****  |Ax-b|=%g,|b|=%g,|du|__=%g\r\n",res,rhs_norm,u_max );	
	}	
GIP_exit:
	;
}	

/*
	修改系统方程:
		type[i]==FIX =>	x[i]=a[i]

	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/7/2008		
*/
void GIP_SS_update_( GIP_SPARSE_SOLVER *hSS,GIP_BOUNDARY_CONDITION *type,double *a )	{
	int i,j,r,nUpdate=0;
	int dim=hSS->dim,*ptr=hSS->ptr,*ind=hSS->ind;
	double *val=hSS->val,*rhs=hSS->rhs;

	for( i = 0; i < dim; i++ )	{
		if( type[i]==GIP_BC_FIX )
			nUpdate++;
	}
	if( nUpdate==0 )
		return;
	else
		G_PRINTF( "\t  ****  GIP_SS_update_:nUpdate=%d\r\n",nUpdate );	

	for( i = 0; i < dim; i++ )	{
		if( type[i]==GIP_BC_FIX )	{
			rhs[i]=a[i];
			for( j=ptr[i]; j < ptr[i+1]; j++ )	{
				r = ind[j];
				if( r==i )
					val[j]=1;
				else 	{
					if( type[r]!=GIP_BC_FIX )
						rhs[r] -= val[j]*a[i];
					val[j]=0.0;
				}
			}		
		}else	{
			for( j=ptr[i]; j < ptr[i+1]; j++ )	{
				r = ind[j];
				if( type[r]==GIP_BC_FIX )		
					val[j]=0.0;
			}
		}
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/8/2008		
*/
void GIP_SS_clear_( GIP_SPARSE_SOLVER *hSS )	{
	if( BIT_TEST( hSS->alg,GIP_SS_DIRECT) && hSS->hLU!=GIP_NULL )	{
		GSS_clear_id( hSS->hLU );
		hSS->hLU=GIP_NULL;
	}
	if( hSS->ptr!=GIP_NULL )
	{	GIP_free(hSS->ptr);			hSS->ptr=GIP_NULL;	}
	if( hSS->ind!=GIP_NULL )
	{	GIP_free(hSS->ind);			hSS->ind=GIP_NULL;	}
	if( hSS->val!=GIP_NULL )
	{	GIP_free(hSS->val);			hSS->val=GIP_NULL;	}
	if( hSS->x!=GIP_NULL )
	{	GIP_free(hSS->x);			hSS->x=GIP_NULL;	}
	if( hSS->rhs!=GIP_NULL )
	{	GIP_free(hSS->rhs);			hSS->rhs=GIP_NULL;	}
	if( hSS->temp!=GIP_NULL )
	{	GIP_free(hSS->temp);		hSS->temp=GIP_NULL;	}
	if( hSS->I_temp!=GIP_NULL )
	{	GIP_free(hSS->I_temp);		hSS->I_temp=GIP_NULL;	}
}
/*
	total variation model for 2-D Images denoise

	1 finite difference scheme
	2 2-D Images(M rows,N columns)	[32] P.13

	ISSUE-NEEDED:
	1、Spectral Total Variation:	["How to discretize the Total Variation of an image"]
	2、Iterative refine 
	3、based on analytical solutions for 4-pixel images ["A Four-Pixel Scheme for Singular Differential Equations"]
	4、high-precision ux,uy(Hamilton-Jacobi WENO)
	5、active set methods
	6、alpha的自动取值 可能the model prefers to remove disks of radius less than alpha
		[REF] "Edge-preserving and Scale-dependent Properties of Total Variation Regularization."
	7、Algorithmic(automatic) differentiation 
	8、Extensions to Total Variation Denoising (1997)
	9、MONTE: The Method of Nonflat Time Evolution in PDE-based Image Restoration
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
#include "../../lib/gss_6_demo.h"
#include "TV_Model.h"
#include "GIP_fmm.h"

#define TVM_POS( r,c,ld )	((r)*(ld)+(c))
#define TVM_R( pos,ld )		((pos)/(ld))
#define TVM_C( pos,ld )		((pos)%(ld))

static int tvm_step=0;

/*
	u_buffer定义为[0:M+1,0:N+1],其中u0为矩形{r=1,r=M,c=1,c=N}

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/2/2008	

void TVM_Init_u0( TV_MODEL *hTVM,GIP_IMAGE *hMask,int type )	{
	int M=hTVM->M,N=hTVM->N,r,c,pos,isOut,ldu=hTVM->ldu;
	int r0=1,r1=M,c0=1,c1=N;		//
	double *u=hTVM->u_buffer,ar,ac,dist;

	if( hMask==0x0 )	{
		for( r = 0; r < M+2; r++ )	{
			for( c = 0; c < N+2; c++ )	{
				pos = TVM_POS( r,c,ldu );
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
					BIT_SET( hTVM->tag[pos],GIP_INTERFACE );	
				}
			}
		}
	}else	{		//based on mask imag
		int height=hMask->height,width=hMask->width,step=hMask->widthStep/sizeof(uchar);
		uchar *i_data = (uchar *)hMask->imageData,a;
		ASSERT( height==M && width==N && hMask->depth==IPL_DEPTH_8U && hMask->nChannels==1 );
		for( r = 0; r < M+2; r++ )	{		//暂定点全部在边界之外
			for( c = 0; c < N+2; c++ )	{
				pos = TVM_POS( r,c,ldu );
				u_buffer[pos]=1.0;				
			}
		}
		for( r = 0; r < height; r++ )	{
			for( c = 0; c < width; c++ )	{
				a = i_data[r*step+c];
				pos = TVM_POS( r+1,c+1,ldu );	
				if( a == 0 )	{
					u_buffer[pos]=0.0;		
					BIT_SET( hTVM->tag[pos],GIP_INTERFACE );
				}else if( a != UCHAR_MAX )	{		//边界内部的点
					u_buffer[pos]=-1.0;
				}
			}
		}
			
	}

}
*/

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		5/22/2008	
*/
int TVM_init( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,TV_MODEL *hTVM )	{
	int M,N,MN,r,c,pos,A_nz,i,setting[32],ret=GIP_OK;
	GIP_FLOATING *val=GIP_NULL;

	tvm_step=0;

	memset( hTVM,0x0,sizeof(TV_MODEL) );
	GIP_import_IplImage( hIMG,&M,&N,&val,0x0 );

	hTVM->M=M;		hTVM->N=N;
	MN=M*N;

	hTVM->ldu=N;			hTVM->shift=0;
	hTVM->u=(double*)GIP_alloc( sizeof(double)*(M)*(N) );
	hTVM->p=(double*)GIP_alloc( sizeof(double)*(M)*(N) );
	hTVM->z=(double*)GIP_alloc( sizeof(double)*(M)*(N) );
	hTVM->tag=(int*)GIP_alloc( sizeof(int)*(M)*(N) );
	memset( hTVM->tag,0x0,sizeof(int)*(M)*(N) );
	hTVM->alpha=5;		hTVM->beta=1.0e-4;
//init the symbol structure of A
//网格行优先编号，对称矩阵A每列储存下列元素{(r,c),(r+1,c),(r,c)}
	hTVM->A_dim=MN;			
	hTVM->m_type=0;
	A_nz=hTVM->m_type==11 || hTVM->m_type==12 ? MN*3 : MN*5;			//取上限
	hTVM->A_ptr=(int*)GIP_alloc( sizeof(int)*(hTVM->A_dim+1) );
	hTVM->A_ind=(int*)GIP_alloc( sizeof(int)*(A_nz) );
	hTVM->A_val=(double*)GIP_alloc( sizeof(double)*(A_nz) );
	hTVM->A_ptr[0]=0;			A_nz=0;
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = TVM_POS( r,c,hTVM->ldu );
			hTVM->z[pos] = val[TVM_POS(r,c,N)];
			if( hTVM->m_type!=11 && hTVM->m_type!=12 )	{
				if( r>0 )	
					hTVM->A_ind[A_nz++]=TVM_POS( r-1,c,hTVM->ldu );
				if( c>0 )
					hTVM->A_ind[A_nz++]=TVM_POS( r,c-1,hTVM->ldu );
			}
			hTVM->A_ind[A_nz++]=pos;
			if( c<N-1 )
				hTVM->A_ind[A_nz++]=TVM_POS( r,c+1,hTVM->ldu );
			if( r<M-1 )	
				hTVM->A_ind[A_nz++]=TVM_POS( r+1,c,hTVM->ldu );
			hTVM->A_ptr[pos+1]=A_nz;

			if( c==0 )		BIT_SET( hTVM->tag[pos],GIP_BD_LEFT );
			if( c==N-1 )	BIT_SET( hTVM->tag[pos],GIP_BD_RIGHT );
			if( r==0 )		BIT_SET( hTVM->tag[pos],GIP_BD_DOWN );
			if( r==M-1 )	BIT_SET( hTVM->tag[pos],GIP_BD_UP );
		}
	}
	hTVM->A_nz = A_nz;

	memcpy( hTVM->u,hTVM->z,sizeof(double)*M*N );
//	memset( hTVM->u,0x0,sizeof(double)*M*N );

	hTVM->hSS=GIP_NULL;
	for( i = 0; i < 32; i++ )	setting[i]=0.0;
	if( (ret=GSS_init_id(hTVM->A_dim,hTVM->A_dim,hTVM->A_ptr,hTVM->A_ind,GIP_NULL,hTVM->m_type,setting) ) != GRUS_OK )	
	{	GIP_ERROR( "TVM_init",ret,"\tERROR at init GSS solver. ERROR CODE:%d\r\n" );		}
	hTVM->hSS = GSS_symbol_id( hTVM->A_dim,hTVM->A_dim,hTVM->A_ptr,hTVM->A_ind,hTVM->A_val );	
	if( hTVM->hSS == NULL )	
	{	GIP_ERROR( "TVM_init",-2,"\tERROR at SYMBOLIC ANALYSIS.\r\n" );		}

GIP_exit:	
	GIP_free( val );

	return ret;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void TVM_clear( TV_MODEL *hTVM )	{
	if( hTVM == GIP_NULL )	return;

	hTVM->M=-1;		hTVM->N=-1;
	if( hTVM->A_ptr != GIP_NULL )	
	{	GIP_free( hTVM->A_ptr );	hTVM->A_ptr=GIP_NULL; }
	if( hTVM->A_ind != GIP_NULL )	
	{	GIP_free( hTVM->A_ind );	hTVM->A_ind=GIP_NULL; }
	if( hTVM->A_val != GIP_NULL )	
	{	GIP_free( hTVM->A_val );	hTVM->A_val=GIP_NULL; }
	if( hTVM->u != GIP_NULL )	
	{	GIP_free( hTVM->u );		hTVM->u=GIP_NULL; }
	if( hTVM->p != GIP_NULL )	
	{	GIP_free( hTVM->p );		hTVM->p=GIP_NULL; }
	if( hTVM->z != GIP_NULL )	
	{	GIP_free( hTVM->z );		hTVM->z=GIP_NULL; }
	if( hTVM->tag != GIP_NULL )	
	{	GIP_free( hTVM->tag );	hTVM->tag=GIP_NULL; }
}

/*	
	adaptive p

	v0.1	cys
		5/29/2008	
*/
void TVM_update_p( TV_MODEL *hTVM,double *Va )	{
	int M=hTVM->M,N=hTVM->N,ldu=hTVM->ldu,r,c,pos;
	double s=0.9,g_max=0.0,g2,ux,uy,a,b,*p=hTVM->p,g;

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = TVM_POS( r,c,ldu );
			ux = c==N-1 ? 0.0 : Va[pos+1]-Va[pos];		//delta+
			uy = r==M-1 ? 0.0 : Va[pos+ldu]-Va[pos];		//delta+	
			g=sqrt(ux*ux+uy*uy );
			g_max = MAX( g,g_max );
		}
	}

	g2 = s*g_max;
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = TVM_POS( r,c,ldu );
			ux = c==N-1 ? 0.0 : Va[pos+1]-Va[pos];		//delta+
			uy = r==M-1 ? 0.0 : Va[pos+ldu]-Va[pos];		//delta+	
			g=sqrt(ux*ux+uy*uy );
			a = g/g2;		b=(g-g2)/g2;
			p[pos] = 1.5*(1+2*a)*b*b+(1-2*b)*a*a;			
			ASSERT( p[pos]>=1.0 && p[pos] <= 2.0 );
			if( g==0.0 )
 				ASSERT( p[pos]==1.5 );
			if( g >= g2 )
 				 p[pos]=1.0;
		}
	}
}


/*
	Dg的公式参见[53],Va定义于[1:M,1:N] 

	v0.1	cys
		5/28/2008	
*/
void TVM_Dg( TV_MODEL *hTVM,double *Va,double *Dg )	{
	int M=hTVM->M,N=hTVM->N,ldu=hTVM->ldu,r,c,pos,isAdaptP=BIT_TEST(hTVM->flag,GIP_TVM_ADAPTIVE_P);
	double *val=hTVM->A_val,*p=hTVM->p,ux,uy,g,alpha=hTVM->alpha,beta=hTVM->beta;

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = TVM_POS( r,c,ldu );
			ux = c==N-1 ? 0.0 : Va[pos+1]-Va[pos];		//delta+
			uy = r==M-1 ? 0.0 : Va[pos+ldu]-Va[pos];		//delta+
			g=sqrt(ux*ux+uy*uy+beta );
			if( isAdaptP )	{
 				Dg[pos]=pow(g,p[pos]-2.0);	
			}else	{
				Dg[pos]=1.0/g;
			}
		}
	}
}

/*
	Vb=Nh(Va)=Va-


	v0.1	cys
		5/23/2008	

void TVM_Nh( TV_MODEL *hTVM,double *Va,double *Vb,double *temp )	{
	int ldu=hTVM->ldu,M=hTVM->M,N=hTVM->N,r,c,pos,isBound,pos_l,pos_r,pos_u,pos_d;
	double alpha=hTVM->alpha,beta=hTVM->beta;
	double a,ux,uy;

	TVM_Dg( hTVM,temp );

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = TVM_POS( r,c,ldu );
			a = Va[pos];
			isBound = c==0 || c==N-1 || r==0 || r==M-1;
			if( !(c==0 || c==N-1) )	{
				pos_l=TVM_POS(r,c-1,ldu);		pos_r=TVM_POS(r,c+1,ldu);
				ux = temp[pos]*(Va[pos_r]-a)-temp[pos_l]*(a-Va[pos_l]);	
			}else
				ux=0.0;
			if( !(r==0 || r==M-1) )	{
				pos_d=TVM_POS(r-1,c,ldu);		pos_u=TVM_POS(r+1,c,ldu);
				uy = temp[pos]*(Va[pos_u]-a)-temp[pos_d]*(a-Va[pos_d]);		
			}else	{
				uy=0.0;
			}
			Vb[pos] = Va[pos]-alpha*(ux+uy);
		}
	}
}
*/

/*
	to update the linear operator	L(u)
	

	v0.1	cys
		5/22/2008	
*/
void TVM_update_La( TV_MODEL *hTVM,double *Va,double *temp )	{
	int dim=hTVM->A_dim,*ptr=hTVM->A_ptr,*ind=hTVM->A_ind,ldu=hTVM->ldu,M=hTVM->M,N=hTVM->N;
	int i,r,c,pos,r_D,c_D,A_nz,pos_D;
	double *val=hTVM->A_val,ux,uy,g,p,a;
	double alpha=hTVM->alpha,beta=hTVM->beta,a_max=0.0,a_min=DBL_MAX;
	
	TVM_Dg( hTVM,Va,temp );

	A_nz=0;
	for( i = 0; i < dim; i++ )	{
		r_D=TVM_R(i,ldu);	c_D=TVM_C(i,ldu);
		a = 0.0;
		if( r_D > 0 )	{	//at (r_D-1,c_D)
			pos = TVM_POS( r_D-1,c_D,ldu );
			a += temp[pos];
			val[A_nz++]=-alpha*temp[pos];
		}
		if( c_D > 0 )	{	//at (r_D,c_D-1)
			pos = TVM_POS( r_D,c_D-1,ldu );
			a += temp[pos];
			val[A_nz++]=-alpha*temp[pos];
		}
		pos_D = A_nz++;			//its the pos of diagonal
		if( c_D < N-1 )	{		//at (r_D,c_D+1)
			pos = TVM_POS( r_D,c_D,ldu );		a += temp[pos];		
			val[A_nz++]=-alpha*temp[pos];
		}
		if( r_D < M-1 )	{		//at (r_D+1,c_D)
			pos = TVM_POS( r_D,c_D,ldu );		a += temp[pos];		
			val[A_nz++]=-alpha*temp[pos];
		}
		val[pos_D]=1.0+a*alpha;
		ASSERT( A_nz==ptr[i+1] );
	}
	for( i = 0; i < A_nz; i++ )	{
		a_max=MAX(a_max,val[i]);		a_min=MIN(a_min,val[i]);
	}

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		5/22/2008	
*/
void TVM_step( TV_MODEL *hTVM,double *temp )	{
	int i,j,r,ret;
	int dim=hTVM->A_dim,*ptr=hTVM->A_ptr,*ind=hTVM->A_ind,ldu=hTVM->ldu,isCheckRes=1;
	double *val=hTVM->A_val,res=0.0,rhs_norm=0.0,*u,u_max,u_min;

//	ASSERT( hTVM->A_nz==M*N );
	u=hTVM->u;
//	TVM_update_Lu( hTVM,temp );
	TVM_update_La( hTVM,u,temp );
//use GSS to solver A*u{n+1}=f(u{n})
//	if( tvm_step%5==0 )	{
		ret = GSS_numeric_id( hTVM->A_dim,hTVM->A_dim,hTVM->A_ptr,hTVM->A_ind,hTVM->A_val,hTVM->hSS );
		if( ret != GRUS_OK )	{
			hTVM->hSS=NULL;		//必须设置为NULL,GSS已自动释放内存
			GIP_ERROR( "TVM_step",ret,"\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n" );
		}
//	}
	memcpy( u,hTVM->z,sizeof(double)*dim );
	if( isCheckRes==1 )		{
		memcpy( temp,u,sizeof(double)*dim );
		for( i = 0; i < dim; i++ )	
		{	rhs_norm += (temp[i]*temp[i]);		}
		rhs_norm = sqrt(rhs_norm);	
	}
	GSS_solve_id( hTVM->hSS,hTVM->A_dim,hTVM->A_dim,hTVM->A_ptr,hTVM->A_ind,hTVM->A_val,u );
	if( isCheckRes==1 )	{
		u_max=DBL_MIN,		u_min=DBL_MAX;
		for( i = 0; i < dim; i++ )	{
			u_max=MAX( u_max,u[i] );
			u_min=MIN( u_min,u[i] );
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				r = ind[j];
				temp[r] -= val[j]*u[i];
//				if( r != i )
//					temp[i] -= val[j]*u[r];
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
	off=Vb-L(Va)*Va
	Normalized residual.[REF]:	ACCELERATION METHODS FOR TOTAL VARIATION-BASED IMAGE DENOISING

	v0.1	cys
		5/23/2008
*/
double  TVM_NormalRes( TV_MODEL *hTVM,double *Va,double *Vb,double *temp )	{
	int i,j,M=hTVM->M,N=hTVM->N,dim=hTVM->A_dim,*ptr=hTVM->A_ptr,*ind=hTVM->A_ind;
	double off=0.0,a,*val=hTVM->A_val,d;

	TVM_update_La( hTVM,Va,temp );

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

/*
	A testing code 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		5/12/2008	
*/
void TVM_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_0,char* sWndName )	{
	double res=DBL_MAX,*temp,tol=0.0,psnr,delta,delta_0;
	TV_MODEL tvm,*hTVM=&tvm;
	int M,N,t_len,isTrace=1,t_step=1,pos; 

	delta=GIP_IMG_delta( hIMG );
	delta_0=GIP_IMG_delta( hIMG_0 );
	TVM_init( hIMG,GIP_NULL,&tvm );
	M=tvm.M,		N=tvm.N;
	t_len=MAX((M+2)*(N+2)*2,N);
	temp=GIP_alloc( sizeof(double)*t_len );
//	TVM_Nh( hTVM,hTVM->z,temp,temp+M*N );
//	tol = GIP_DOFF2( M*N,hTVM->z,temp );
//	pos = GIP_MAXOFF_POS( M*N,hTVM->z,temp,0x0 );
	memcpy( temp,hTVM->z,sizeof(double)*M*N );
	tol = TVM_NormalRes( hTVM,temp,hTVM->z,temp+M*N );
	tol *= 1.0e-4;
 	psnr = GIP_IMG_PSNR( M,N,hTVM->z,hTVM->ldu,hIMG_0,temp );

	G_PRINTF( "\tdelta=%g,delta_0=%g\r\n",delta,delta_0 );		
	G_PRINTF( "START:\t|z-Nh(z)|=%g,psnr=%g,alpha=%g,beta=%g\r\n",tol*1.0e4,psnr,hTVM->alpha,hTVM->beta );		
//total variation loop
	while( res>tol || BIT_TEST(tvm.flag,GIP_SOLVER_CONVERGE) )	{
		if( isTrace==1 && tvm_step%t_step==0 )	
		{	GIP_interactive( M,N,tvm.u+tvm.shift,tvm.ldu,tvm_step,&t_step,0x0 );	}
		TVM_step( &tvm,temp );
// 		TVM_Nh( hTVM,hTVM->u,temp,temp+M*N );
//		res=GIP_DOFF2( M*N,hTVM->z,temp );
		res = TVM_NormalRes( hTVM,hTVM->u,hTVM->z,temp );
		psnr = GIP_IMG_PSNR( M,N,hTVM->u,hTVM->ldu,hIMG_0,temp );
		G_PRINTF( "%d:\t|z-Nh(u)|=%g,psnr=%g\r\n",tvm_step,res,psnr );		
		tvm_step++;
	}
/*adaptive total variation loop
	BIT_SET( tvm.flag,GIP_TVM_ADAPTIVE_P );
	TVM_update_p( hTVM,hTVM->u );
	memcpy( temp,hTVM->z,sizeof(double)*M*N );
	tol = TVM_NormalRes( hTVM,temp,hTVM->z,temp+M*N );
	tol *= 1.0e-5;
	G_PRINTF( "START:\t|z-Nh(z)|=%g,psnr=%g,alpha=%g,beta=%g\r\n",tol*1.0e4,psnr,hTVM->alpha,hTVM->beta );		
	while( res>tol || BIT_TEST(tvm.flag,GIP_SOLVER_CONVERGE) )	{
		if( isTrace==1 && tvm_step%t_step==0 )	
		{	GIP_interactive( M,N,tvm.u+tvm.shift,tvm.ldu,tvm_step,&t_step,0x0 );	}
		TVM_step( &tvm,temp );
// 		TVM_Nh( hTVM,hTVM->u,temp,temp+M*N );
//		res=GIP_DOFF2( M*N,hTVM->z,temp );
		res = TVM_NormalRes( hTVM,hTVM->u,hTVM->z,temp );
		psnr = GIP_IMG_PSNR( M,N,hTVM->u,hTVM->ldu,hIMG_0,temp );
		G_PRINTF( "%d:\t|z-Nh(u)|=%g,psnr=%g\r\n",tvm_step,res,psnr );		
		tvm_step++;
	}
*/
	G_PRINTF( "GAC iteration finish\r\n"	);
	GIP_free( temp );
GIP_exit:
	if( isTrace==1 )	{
		char sTitle[80];
		_stprintf( sTitle,".\\trace\\u_%d.bmp",tvm_step );		
		GIP_save_IplImage_d( M,N,sTitle,tvm.u+tvm.shift,tvm.ldu,0x0 );
	}

	TVM_clear( &tvm );
}
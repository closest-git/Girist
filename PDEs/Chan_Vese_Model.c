/*
	chan-vese model for 2-D Images

	1 finite difference scheme
	2 2-D Images(M rows,N columns)	[32] P.13

	ISSUE-NEEDED:
	1、reinitialize u to the signed distance function. prevent interior contour from growing.
	2、mui，the scale(length) parameter
	3、use other infomation from initial imag(curvature or orientation...)
	4、used the k-means segmentation as an initial guess;	Canny and the SUSAN edge detectors 
	5、Gibou-Fedkiw algorithm.	"A fast level set based algorithm for segmentation"	2005
	6、Song-Chan algorithm.		"Fast algorithm for level set based optimization"	2002
	7、Some Recent Developments in Variational Image Segmentation
	8、Image Segmentation by Piecewise Constant Mumford-Shah Model without Estimating the Constants
	9、Median Scheme "MODERN THEORY OF NUMERICAL METHODS FOR MOTION BY MEAN CURVATURE"
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
#include "Chan_Vese_Model.h"
#include "GIP_fmm.h"
#include "../../lib/gss_6_demo.h"

static int cvm_step=0;

/*
	u_buffer定义为[0:M+1,0:N+1],其中u0为矩形{r=1,r=M,c=1,c=N}

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/2/2008	
*/
void CVM_Init_u0( CV_MODEL *hCVM,GIP_IMAGE *hMask,int type )	{
	int M=hCVM->M,N=hCVM->N,r,c,pos,isOut,ldu=hCVM->ldu;
	int r0=10,r1=M-10,c0=20,c1=N-20;		//
	double *u_buffer=hCVM->u_buffer,ar,ac,dist;

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
				u_buffer[pos] = (isOut==1) ? dist : -dist;
				if( u_buffer[pos]==0.0 )			{	
					BIT_SET( hCVM->tag[pos],GIP_INTERFACE );	
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
				u_buffer[pos]=1.0;				
			}
		}
		for( r = 0; r < height; r++ )	{
			for( c = 0; c < width; c++ )	{
				a = i_data[r*step+c];
				pos = GIP_RC2POS( r+1,c+1,ldu );	
				if( a == 0 )	{
					u_buffer[pos]=0.0;		
					BIT_SET( hCVM->tag[pos],GIP_INTERFACE );
				}else if( a != UCHAR_MAX )	{		//边界内部的点
					u_buffer[pos]=-1.0;
				}
			}
		}
			
	}

}

/*
	computes the reinitialized distance function in a band 
	that is one-layer thick on each side of the interface.
	[REF]: The Fast Construction of Extension Velocities in Level Set Methods, D. Adalsteinsson and J. A. Sethian, J. Comput. Phys. 148, 2 (1999).
	[REF]: A Study of Numerical Methods for the Level Set Approach
	
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/21/2008	
	v0.2	cys
		5/6/2008	

void CVM_ReInit_u0( CV_MODEL *hCVM,int type )	{
	int M=hCVM->M,N=hCVM->N,r,c,pos,isOut,ldu=hCVM->ldu,*tag=hCVM->tag,fix;
	int neibor[]={-1,1,-ldu,ldu},nCross,i,sdf_alg=GIP_SDF_INTERSECTION;
	double *u_buffer=hCVM->u_buffer,a,b,cross[4],s,t;

	for( r = 0; r < M+2; r++ )	{
	for( c = 0; c < N+2; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
		tag[pos] = tag[pos]&TAG_INIT_FIX;
	}
	}

	for( r = 1; r < M+1; r++ )	{
		for( c = 1; c < N+1; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			a = u_buffer[pos];
			if( a==0.0 )	{
				BIT_SET( tag[pos],GIP_CVM_INTERFACE );	
			}else {
				nCross = 0;			s=DBL_MAX;		t=DBL_MAX;
				for( i = GIP_NEIBOR_LEFT; i <= GIP_NEIBOR_UP; i++ )	{
					b = u_buffer[pos+neibor[i]];
					if( a*b<0.0 )			
					{	cross[i]=fabs(a)/( fabs(a)+fabs(b) );	nCross++;	}
					else
						cross[i]=DBL_MAX;
				}
				if( nCross == 0 )	
					continue;
				BIT_SET( tag[pos],GIP_CVM_INTERFACE );	
				switch( sdf_alg )	{
				case GIP_SDF_INTERSECTION:
					s = MIN(cross[0],cross[1]);			t = MIN(cross[2],cross[3]);		//X,Y	
					if( s>1.0 || t>1.0 )	{	
						ASSERT( nCross==1 || nCross==2 );		//Figure 4b 4c 4e
						u_buffer[pos] =	MIN(s,t);		//Figure 4a 4d
					}	else	{
						ASSERT( s<1.0 && t < 1.0 );		//Figure 4b 4c 4e
						u_buffer[pos] =	s*t/sqrt(s*s+t*t);
					}
				break;
				case GIP_SDF_PROJECT_1:	{
					double *u,ux,uy,uxx,uyy,uxy,b,c,alpha;
					u = u_buffer+pos;
					ux=(u[1]-u[-1])/2;				uy=(u[ldu]-u[-ldu])/2;
					uxx=(u[1]-2*u[0]+u[-1]);		uyy=(u[ldu]-2*u[0]+u[-ldu]);
					uxy=((u[1+ldu]-u[1-ldu])/2-(u[-1+ldu]-u[-1-ldu])/2)/2;	
					b = ux*ux+uy*uy;
					c = (ux*ux*uxx+2*ux*uy*uxy+uy*uy*uyy)/2.0;
					ASSERT( (b*b-4*a*c)>=0.0 );
					alpha = (-b+sqrt(b*b-4*a*c))/2.0/c;
					u_buffer[pos] = fabs(alpha*sqrt(b));										}
				break;
				}
				}
				u_buffer[pos] = a > 0.0 ? u_buffer[pos] : -u_buffer[pos];	
//				ASSERT( -1.0 <= u_buffer[pos] && u_buffer[pos]<= 1.0 );
		}
	}
}
*/

/*
	ISSUE-NEEDED
	1、heavyside function

	Copyright 2008-present, Grusoft.

	v0.1	cys
		5/30/2008	
*/
void CVM_Update_c( CV_MODEL *hCVM )	{
	int M=hCVM->M,N=hCVM->N,MN=M*N,r,c,pos,ldu=hCVM->ldu,n1=0,n2=0;
	double c1=0.0,c2=0.0,*u=hCVM->u_buffer,*z=hCVM->z,h,h1,h2,eps=hCVM->eps;
	double thresh=hCVM->alg==CVM_GLOBAL_CONVEX ? 0.5 : 0.0;

	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{		
		pos = GIP_RC2POS( r,c,ldu );
		if( u[pos]>=thresh )	{
			c1 += z[pos];
			n1++;
		}else	{
			c2 += z[pos];
			n2++;
		}
	}
	}
	c1/=n1;		c2/=n2;
/*	h1=0.0,		h2=0.0;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{		
		pos = GIP_RC2POS( r,c,ldu );
		h = 0.5*(1+2.0/GIP_PI*atan(u[pos]/eps) );
		c1 += z[pos]*h;
		c2 += z[pos]*(1-h);
		h1 += h;		h2+=1-h;
	}
	}
	c1/=h1;		c2/=h2;*/
	hCVM->c1=c1;		hCVM->c2=c2;

}

/*
	Init edge dector function from Intensity Image
	
	1	edge dector function    from [10] P.153
	2	边界处理：暂时空出

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/19/2008	
*/
void CVM_init( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,CV_MODEL *hCVM )	{
	int M,N,MN,r,c,pos,A_nz,i,setting[32],ret=GIP_OK;
	GIP_FLOATING *val=GIP_NULL;
	double *g,a;

	cvm_step=0;

	memset( hCVM,0x0,sizeof(CV_MODEL) );
	GIP_import_IplImage( hIMG,&M,&N,&val,0x0 );

	hCVM->alg = CVM_ALG_0;			//CVM_GLOBAL_CONVEX,CVM_ALG_0
	hCVM->lenda1=hCVM->lenda2=1.0;
	hCVM->eps = 1.0;
	hCVM->dt = 0.1*MIN(M,N);		//0.01;	
	hCVM->mui = 0.001*255*255;
	hCVM->M=M;		hCVM->N=N;
	MN=M*N;
	hCVM->ldu=N;			hCVM->shift=0;
	hCVM->u_buffer=(double*)GIP_alloc( sizeof(double)*(M+2)*(N+2) );
	hCVM->z=(double*)GIP_alloc( sizeof(double)*(M+2)*(N+2) );
	hCVM->tag=(int*)GIP_alloc( sizeof(int)*(M+2)*(N+2) );
	memset( hCVM->tag,0x0,sizeof(int)*(M+2)*(N+2) );

	hCVM->A_dim=MN;			
	hCVM->m_type=0;
	A_nz=hCVM->m_type==11 || hCVM->m_type==12 ? MN*3 : MN*5;			//取上限
	hCVM->A_ptr=(int*)GIP_alloc( sizeof(int)*(hCVM->A_dim+1) );
	hCVM->A_ind=(int*)GIP_alloc( sizeof(int)*(A_nz) );
	hCVM->A_val=(double*)GIP_alloc( sizeof(double)*(A_nz) );
	hCVM->A_ptr[0]=0;			A_nz=0;
	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,hCVM->ldu );
			hCVM->z[pos] = val[GIP_RC2POS(r,c,N)];
			if( hCVM->m_type!=11 && hCVM->m_type!=12 )	{
				if( r>0 )	
					hCVM->A_ind[A_nz++]=GIP_RC2POS( r-1,c,hCVM->ldu );
				if( c>0 )
					hCVM->A_ind[A_nz++]=GIP_RC2POS( r,c-1,hCVM->ldu );
			}
			hCVM->A_ind[A_nz++]=pos;
			if( c<N-1 )
				hCVM->A_ind[A_nz++]=GIP_RC2POS( r,c+1,hCVM->ldu );
			if( r<M-1 )	
				hCVM->A_ind[A_nz++]=GIP_RC2POS( r+1,c,hCVM->ldu );
			hCVM->A_ptr[pos+1]=A_nz;

			if( c==0 )		BIT_SET( hCVM->tag[pos],GIP_BD_LEFT );
			if( c==N-1 )	BIT_SET( hCVM->tag[pos],GIP_BD_RIGHT );
			if( r==0 )		BIT_SET( hCVM->tag[pos],GIP_BD_DOWN );
			if( r==M-1 )	BIT_SET( hCVM->tag[pos],GIP_BD_UP );
		}
	}
	hCVM->A_nz = A_nz;
	hCVM->hSS=GIP_NULL;
	for( i = 0; i < 32; i++ )	setting[i]=0.0;
	if( (ret=GSS_init_id(hCVM->A_dim,hCVM->A_dim,hCVM->A_ptr,hCVM->A_ind,GIP_NULL,hCVM->m_type,setting) ) != GRUS_OK )	
	{	GIP_ERROR( "TVM_init",ret,"\tERROR at init GSS solver. ERROR CODE:%d\r\n" );		}
	hCVM->hSS = GSS_symbol_id( hCVM->A_dim,hCVM->A_dim,hCVM->A_ptr,hCVM->A_ind,hCVM->A_val );	
	if( hCVM->hSS == NULL )	
	{	GIP_ERROR( "TVM_init",-2,"\tERROR at SYMBOLIC ANALYSIS.\r\n" );		}


	CVM_Init_u0( hCVM,hIMG_mask,0x0 );
	GIP_FMM_smooth_u_( hCVM->M,hCVM->N,hCVM->u_buffer,hCVM->ldu,hCVM->tag,0x0 );
	CVM_Update_c( hCVM );
GIP_exit:	
	GIP_free( val );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	
*/
void CVM_clear( CV_MODEL *hCVM )	{
	if( hCVM == GIP_NULL )	return;

	hCVM->M=-1;		hCVM->N=-1;
	if( hCVM->A_ptr != GIP_NULL )	
	{	GIP_free( hCVM->A_ptr );	hCVM->A_ptr=GIP_NULL; }
	if( hCVM->A_ind != GIP_NULL )	
	{	GIP_free( hCVM->A_ind );	hCVM->A_ind=GIP_NULL; }
	if( hCVM->A_val != GIP_NULL )	
	{	GIP_free( hCVM->A_val );	hCVM->A_val=GIP_NULL; }
	if( hCVM->u_buffer != GIP_NULL )	
	{	GIP_free( hCVM->u_buffer );	hCVM->u_buffer=GIP_NULL; }
	if( hCVM->z != GIP_NULL )	
	{	GIP_free( hCVM->z );	hCVM->z=GIP_NULL; }
	if( hCVM->tag != GIP_NULL )	
	{	GIP_free( hCVM->tag );	hCVM->tag=GIP_NULL; }
}

/*

	v0.1	cys
		5/30/2008	
*/
void CVM_Dg( CV_MODEL *hCVM,double *Va,double *Dg )	{
	int M=hCVM->M,N=hCVM->N,ldu=hCVM->ldu,r,c,pos;
	double ux,uy,g,beta=1.0E-12;

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,ldu );
			ux = c==N-1 ? 0.0 : Va[pos+1]-Va[pos];		//delta+
			uy = r==M-1 ? 0.0 : Va[pos+ldu]-Va[pos];		//delta+
			g=sqrt(ux*ux+uy*uy+beta );
			Dg[pos]=1.0/g;
		}
	}
}

/*

	v0.1	cys
		5/30/2008	
*/
void CVM_update_La( CV_MODEL *hCVM,double *Va,double *temp )	{
	int dim=hCVM->A_dim,*ptr=hCVM->A_ptr,*ind=hCVM->A_ind,ldu=hCVM->ldu,M=hCVM->M,N=hCVM->N;
	int i,r,c,pos,r_D,c_D,A_nz,pos_D,alg=hCVM->alg;
	double *val=hCVM->A_val,*u,a,s,s_0,eps,delta,a_max,a_min;
	
	CVM_Dg( hCVM,Va,temp );

	u = hCVM->u_buffer;
	A_nz=0;
	eps=hCVM->eps;
	s_0 = hCVM->dt*hCVM->mui;
	for( i = 0; i < dim; i++ )	{
		r_D=GIP_POS2R(i,ldu);	c_D=GIP_POS2C(i,ldu);
		if( alg == CVM_GLOBAL_CONVEX )
			delta = 1.0;
		else
			delta = eps/( (eps*eps+u[i]*u[i]) * GIP_PI );
		s = s_0*delta;
		a = 0.0;
		if( r_D > 0 )	{	//at (r_D-1,c_D)
			pos = GIP_RC2POS( r_D-1,c_D,ldu );
			a += temp[pos];
			val[A_nz++]=-s*temp[pos];
		}
		if( c_D > 0 )	{	//at (r_D,c_D-1)
			pos = GIP_RC2POS( r_D,c_D-1,ldu );
			a += temp[pos];
			val[A_nz++]=-s*temp[pos];
		}
		pos_D = A_nz++;			//its the pos of diagonal
		if( c_D < N-1 )	{		//at (r_D,c_D+1)
			pos = GIP_RC2POS( r_D,c_D,ldu );		a += temp[pos];		
			val[A_nz++]=-s*temp[pos];
		}
		if( r_D < M-1 )	{		//at (r_D+1,c_D)
			pos = GIP_RC2POS( r_D,c_D,ldu );		a += temp[pos];		
			val[A_nz++]=-s*temp[pos];
		}
		val[pos_D]=1.0+a*s;
		ASSERT( A_nz==ptr[i+1] );
	}

	a_max=DBL_MIN;		a_min=DBL_MAX;
	for( i = 0; i < A_nz; i++ )	{
		a_max=MAX(a_max,val[i]);		a_min=MIN(a_min,val[i]);
	}

}

/*
	rhs=u+dt*delta*mu(lenda2-lenda1)		[17]

	v0.1	cys
		5/30/2008	
*/
void CVM_update_rhs( CV_MODEL *hCVM,double *rhs,double *temp )	{
	int M=hCVM->M,N=hCVM->N,MN=M*N,r,c,pos,ldu=hCVM->ldu,alg=hCVM->alg;
	double *u=hCVM->u_buffer,*z=hCVM->z,lenda1,lenda2,s,a,b,d,c1,c2,eps,delta;
	double eps_3=1.0e-2,eps_4=eps_3/sqrt(2.0),silen,alpha,s_max=0.0;

	lenda1=hCVM->lenda1,		lenda2=hCVM->lenda2;
	c1=hCVM->c1,				c2=hCVM->c2;
	eps=hCVM->eps;

	s = hCVM->dt;
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
		if( alg == CVM_GLOBAL_CONVEX )	{
			delta = 1.0;
/*			if( a < -eps_4 )	silen=-1.0;
			else if( -eps_4<=a && a<eps_3 ) 
				silen=-1.0;
			else if( eps_3<=a && a<1-eps_3 )	silen=0.0;
			else if( 1-eps_3<=a && a<1+eps_4 ) 
				silen=1.0;
			else if( 1.0+eps_4<= a )	silen=1.0;
			d -= alpha*silen;*/
		}
		else
			delta = eps/( (eps*eps+a*a) * GIP_PI );
		d *= delta*s;
		rhs[pos] = a+d*s;
	}
	}
	if( alg == CVM_GLOBAL_CONVEX )	{
		alpha = MAX(lenda1,lenda2)*s_max*0.6;
		for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{	
			pos = GIP_RC2POS( r,c,ldu );
			a = u[pos];
			if( a < -eps_4 )	silen=-1.0;
			else if( -eps_4<=a && a<eps_3 ) 
				silen=-1.0;
			else if( eps_3<=a && a<1-eps_3 )	silen=0.0;
			else if( 1-eps_3<=a && a<1+eps_4 ) 
				silen=1.0;
			else if( 1.0+eps_4<= a )	silen=1.0;
			rhs[pos] -= alpha*silen*s;
		}
		}
	}

	return;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		5/13/2008	
*/
void CVM_step( CV_MODEL *hCVM,double dT,double *temp )	{
	int M=hCVM->M,N=hCVM->N,MN=M*N,i,j,r,c,g_pos,pos,ldu=hCVM->ldu,*tag=hCVM->tag,ret=GIP_OK;
	int dim=hCVM->A_dim,*ptr=hCVM->A_ptr,*ind=hCVM->A_ind,isCheckRes=1;
	double *rhs,*u,u_max,u_min,rhs_norm=0.0,*val=hCVM->A_val,res;

	u=hCVM->u_buffer;
	CVM_update_La( hCVM,u,temp );
//use GSS to solver A*u{n+1}=f(u{n})
//	if( cvm_step%5==0 )	{
		ret = GSS_numeric_id( hCVM->A_dim,hCVM->A_dim,hCVM->A_ptr,hCVM->A_ind,hCVM->A_val,hCVM->hSS );
		if( ret != GRUS_OK )	{
			hCVM->hSS=NULL;		//必须设置为NULL,GSS已自动释放内存
			GIP_ERROR( "TVM_step",ret,"\r\n\tERROR at NUMERIC FACTORIZATION. ERROR CODE:%d\r\n" );
		}
//	}
	rhs=temp;
	CVM_update_rhs( hCVM,rhs,temp+M*N );
	memcpy( u,rhs,sizeof(double)*dim );
	if( isCheckRes==1 )		{
		memcpy( temp,u,sizeof(double)*dim );
		for( i = 0; i < dim; i++ )	
		{	rhs_norm += (temp[i]*temp[i]);		}
		rhs_norm = sqrt(rhs_norm);	
	}
	GSS_solve_id( hCVM->hSS,hCVM->A_dim,hCVM->A_dim,hCVM->A_ptr,hCVM->A_ind,hCVM->A_val,u );
	if( isCheckRes==1 )	{
		res = 0.0;
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
	CVM_Update_c( hCVM );

GIP_exit:
	;
}

/*
	A testing code 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		5/12/2008	
*/
void CVM_test( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_mask,char* sWndName )	{
	double *temp;
	CV_MODEL cvm,*hCVM=&cvm;
	int M,N,t_len,isTrace=1,nKey,t_step=1;
	char sTitle[80];

	CVM_init( hIMG,hIMG_mask,&cvm );
	G_PRINTF( "START:\tmui=%g,dT=%g,eps=%g\r\n",hCVM->mui,hCVM->dt,hCVM->eps );		

	M=cvm.M,		N=cvm.N;
	t_len=MAX((M+2)*(N+2)*2,N);
	temp=GIP_alloc( sizeof(double)*t_len );
//	GIP_show_IplImage_d( M,N,sWndName,cvm.z,cvm.ldu,0x0 );
	while( !BIT_TEST(cvm.flag,GIP_SOLVER_CONVERGE) )	{
		if( isTrace==1 && cvm_step%t_step==0 )	{
			_stprintf( sTitle,".\\trace\\u_%d.bmp",cvm_step );		
			GIP_save_IplImage_d( M,N,sTitle,(cvm.u_buffer+cvm.shift),cvm.ldu,0x0 );
//			_stprintf( sTitle,"u_%d.bmp",cvm_step );		
//			GIP_show_IplImage_d( M,N,sTitle,(cvm.u_buffer+cvm.shift),cvm.ldu,GIP_OPENCV_CREATEWND );
//			_stprintf( sTitle,".\\trace\\interface_%d.bmp",step );		
			GIP_CONTOUR_1( M,N,cvm.u_buffer+cvm.shift,temp,cvm.ldu,0x0 );
			GIP_show_IplImage_d( M,N,sWndName,temp,cvm.ldu,0x0 );
			nKey = cvWaitKey( );
			switch( nKey )	{
				case 27:
					G_PRINTF( "GAC iteration break out\r\n"	);
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
//			GIP_save_IplImage_d( M,N,sTitle,(temp+cvm.shift),cvm.ldu,0x0 );
		}
		G_PRINTF( "%d:\tc1=%g,c2=%g\r\n",cvm_step,hCVM->c1,hCVM->c2 );
//		average = CVM_norm_Hu_( M,N,hCVM->u_buffer+hCVM->shift,hCVM->ldu,temp+hCVM->shift );		
		CVM_step( &cvm,0.0,temp );
 		if( BIT_TEST(cvm.flag,GIP_SOLVER_CONVERGE) )
			break;
		cvm_step++;
		//if( average > 2 )	{
		//	CVM_ReInit_u0( &cvm,0x0 );
		//	CVM_smooth_u_( &cvm );
		//}
	}
	G_PRINTF( "GAC iteration finish\r\n"	);
	GIP_free( temp );
GIP_exit:
	CVM_clear( &cvm );
}
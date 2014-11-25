#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <FLOAT.h>
#include <limits.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "GIP_fem.h"
#include "GIP_mesh.h"
#include "../util/GIP_util.h"
#include "../util/GIP_ss.h"

#define _FEM_RANK_MAX 9
#define _DU_1_KSI		1.0e-2					//用于避免|Du|=0,理论上越大，越难收敛，还需具体分析。
//_DU_ALPHA尺度系数：越大，则细节越少。(|du|:edge-preserving)
#define _DU_ALPHA		3					
#define _DU2_BETA		0					//越大，渲染越明显

/*	
	节省内存
		K,J可以不用存储,通过装配生成. hTri->K可以采用索引存储
*/

/*
	N(形状函数)由其系数代表
	[22] P.49

	注意：
		1、和标准公式相比，N并没有/2A

	Copyright 2008-present, Grusoft.
	v0.1	9/26/2008		
*/
void _FEM_tri_init( GIP_FEM_TRI* hTri,GIP_FEM_*hFEM,int no_1,int no_2,int no_3 )	{
	int i,j;
	GIP_FLOATING s_C,x1,x2,x3,y1,y2,y3,D=0;
	GIP_FEM_NODE *nodes=hFEM->arrNode;

	hTri->nodes[0]=no_1;		hTri->nodes[1]=no_2;	hTri->nodes[2]=no_3;		
	hTri->flag[0]=0;			hTri->flag[1]=0;		hTri->flag[2]=0;		
	hTri->f[0]=0;				hTri->f[1]=0;			hTri->f[2]=0;
	x1=nodes[no_1].c;			x2=nodes[no_2].c;		x3=nodes[no_3].c;
	y1=nodes[no_1].r;			y2=nodes[no_2].r;		y3=nodes[no_3].r;
//	hTri->B[0]=y2-y3;			hTri->B[1]=y3-y1;		hTri->B[2]=y1-y2;
//	hTri->B[3]=x3-x2;			hTri->B[4]=x1-x3;		hTri->B[5]=x2-x1;
	hTri->N[0]=x2*y3-x3*y2;		hTri->N[1]=x3*y1-x1*y3;		hTri->N[2]=x1*y2-x2*y1;
	hTri->N[3]=y2-y3;			hTri->N[4]=y3-y1;			hTri->N[5]=y1-y2;
	hTri->N[6]=x3-x2;			hTri->N[7]=x1-x3;			hTri->N[8]=x2-x1;

	hTri->area=( (x1*y2-x2*y1)-(x1*y3-x3*y1)+(x2*y3-x3*y2) )/2;
	ASSERT( hTri->area==(hTri->N[0]+hTri->N[1]+hTri->N[2])/2 );

//	ASSERT( fabs(hTri->area-0.5)<DBL_EPSILON );
/*	b=hTri->B;		c=hTri->B+3;
	u=hFEM->u;
	switch( hFEM->model )	{
		case GIP_MODEL_130:
			D=200;
			break;
		case GIP_MODEL_CGM:
			Dx=(b[0]*u[no_1]+b[1]*u[no_2]+b[2]*u[no_3])/2/hTri->area;
			Dy=(c[0]*u[no_1]+c[1]*u[no_2]+c[2]*u[no_3])/2/hTri->area;
			D=1.0/sqrt( (Dx*Dx+Dy*Dy)+1.0e-4 );
			break;
		default:
			ASSERT( 0 );
			break;
	}*/
	s_C=hTri->area/12.0;
//	s_K=1.0/4/hTri->area*D;
	for( i = 0; i < 3; i++ )	{
	for( j = 0; j < 3; j++ )	{
//		hTri->K[i*3+j]=(b[i]*b[j]+c[i]*c[j])*s_K;
		hTri->C[i*3+j]=(i==j)?s_C*2:s_C;
	}
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/11/2008		
*/
double _FEM_DU( GIP_FEM_ *hFEM,int flag ){
	int e,rank=3,no_1,no_2,no_3;
	GIP_FEM_TRI* hEle;
	GIP_FLOATING *b,*c,*u;
	double dU=0.0,Ux,Uy;

	u=hFEM->u;
	for( e = 0; e < hFEM->nEle; e++ )	{
		hEle = (GIP_FEM_TRI*)(hFEM->arrEle)+e;	
		b=hEle->N+rank;				c=hEle->N+rank*2;
		no_1=hEle->nodes[0];		no_2=hEle->nodes[1];		no_3=hEle->nodes[2];		
		Ux=(b[0]*u[no_1]+b[1]*u[no_2]+b[2]*u[no_3])/2/hEle->area;
		Uy=(c[0]*u[no_1]+c[1]*u[no_2]+c[2]*u[no_3])/2/hEle->area;
		dU+=sqrt( (Ux*Ux+Uy*Uy)+_DU_1_KSI );
	}
//	dU = sqrt(dU);

	return dU;
}

/*
[**]
		雪	梦中偶得
	飞玉润寒蕊，滴水下苍苔。
	风吹柴门启，见友携酒来
[**]

	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/14/2008		
*/
void _FEM_form_K( GIP_FEM_ *hFEM,int flag )	{
	int e,i,j,rank=3,no_1,no_2,no_3;
	GIP_FEM_TRI* hEle;
	GIP_FLOATING *b,*c,*u;
	double s_K,Dx,Dy,D,nn,r_sum,u_scale=hFEM->u_scale;

	u=hFEM->u;
	for( e = 0; e < hFEM->nEle; e++ )	{
		hEle = (GIP_FEM_TRI*)(hFEM->arrEle)+e;	
//		b=hEle->B;					c=hEle->B+3;
		b=hEle->N+rank;				c=hEle->N+rank*2;
		no_1=hEle->nodes[0];		no_2=hEle->nodes[1];		no_3=hEle->nodes[2];		
		switch( hFEM->model )	{
			case GIP_MODEL_130:
				D=200;
				break;
			case GIP_SEG_MUMSHAH_2:
			case GIP_SEG_MUMSHAH_s:
			case GIP_MUMSHAH_1:
			case GIP_MUMSHAH_filter:
				Dx=(b[0]*u[no_1]+b[1]*u[no_2]+b[2]*u[no_3])/2/hEle->area;
				Dy=(c[0]*u[no_1]+c[1]*u[no_2]+c[2]*u[no_3])/2/hEle->area;
				D=_DU_ALPHA*1.0/sqrt( (Dx*Dx+Dy*Dy)+_DU_1_KSI )+_DU2_BETA;
				break;
			case GIP_MODEL_CGM:
//			case GIP_DECOMP_Hs:
			case GIP_SEG_CV:
				Dx=(b[0]*u[no_1]+b[1]*u[no_2]+b[2]*u[no_3])/2/hEle->area;
				Dy=(c[0]*u[no_1]+c[1]*u[no_2]+c[2]*u[no_3])/2/hEle->area;
				D=_DU_ALPHA*1.0/sqrt( (Dx*Dx+Dy*Dy)+_DU_1_KSI );
				break;
			default:
				ASSERT( 0 );
				break;
		}
//		s_K=1.0/4/hEle->area*(_DU_ALPHA*D+_DU_BETA);
		s_K=1.0/4/hEle->area*D;
		for( i = 0; i < rank; i++ )	{
			r_sum=0.0;
			for( j = 0; j < rank; j++ )	{
				hEle->K[i*rank+j]=(b[i]*b[j]+c[i]*c[j])*s_K;
				if( (hFEM->model==GIP_MODEL_CGM || hFEM->model==GIP_SEG_MUMSHAH_2) )	{
					nn = (i==j)?hEle->area/6.0:hEle->area/12.0;
					hEle->K[i*rank+j]+=nn;
				}else if( BIT_TEST(hFEM->method,_FEM_U_SCALE) )	{
					nn = (i==j)?hEle->area/6.0:hEle->area/12.0;
					hEle->K[i*rank+j]+=nn*u_scale;
				}
				r_sum += hEle->K[i*rank+j];
			}
//			ASSERT( r_sum==0.0 );		//	!!testing!!		
		}		
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/14/2008		
*/
void _FEM_form_J( GIP_FEM_ *hFEM,int e,double *J,int flag )	{
	int i,j,rank=3,no_1,no_2,no_3;
	GIP_FEM_TRI* hEle;
	GIP_FLOATING *b,*c,*u;
	double s,s_K,Ux,Uy,dU,k_u_;

	u=hFEM->u;
//	for( e = 0; e < hFEM->nEle; e++ )	{
		hEle = (GIP_FEM_TRI*)(hFEM->arrEle)+e;	
//		b=hEle->B;					c=hEle->B+3;
		b=hEle->N+rank;				c=hEle->N+rank*2;
		no_1=hEle->nodes[0];		no_2=hEle->nodes[1];		no_3=hEle->nodes[2];		
		switch( hFEM->model )	{
			case GIP_SEG_MUMSHAH_2:
			case GIP_SEG_MUMSHAH_s:
			case GIP_MUMSHAH_1:
			case GIP_MUMSHAH_filter:
			case GIP_MODEL_CGM:
//			case GIP_DECOMP_Hs:
			case  GIP_SEG_CV:
				Ux=(b[0]*u[no_1]+b[1]*u[no_2]+b[2]*u[no_3])/2/hEle->area;
				Uy=(c[0]*u[no_1]+c[1]*u[no_2]+c[2]*u[no_3])/2/hEle->area;
				dU=sqrt( (Ux*Ux+Uy*Uy)+_DU_1_KSI );
				break;
			default:
				ASSERT( 0 );
				break;
		}
		s_K=1.0/4/hEle->area;
		for( i = 0; i < rank; i++ )	{
			k_u_ = 0.0;
			for( j = 0; j < rank; j++ )	{
				k_u_ += (b[i]*b[j]+c[i]*c[j])*u[hEle->nodes[j]]*s_K;
			}
			for( j = 0; j < rank; j++ )	{
				s = -_DU_ALPHA*(Ux*b[j]+Uy*c[j])/dU/dU/dU;
				J[i*rank+j]=hEle->K[i*rank+j]+s*k_u_;
			}
		}		
//	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/7/2008		
*/
void  _FEM_form_CU( GIP_FEM_ *hFEM,double *cu,int flag )	{
	int dim=hFEM->K_dim,e,i,j,rank=3,n_i,n_j;
	GIP_FEM_TRI* hEle;
	GIP_FLOATING *u=hFEM->u;
	double nn,u_scale=hFEM->u_scale;

	for( i = 0; i < dim; i++ )	cu[i]=0.0;

	for( e = 0; e < hFEM->nEle; e++ )	{
		hEle = (GIP_FEM_TRI*)(hFEM->arrEle)+e;	
		rank = 3;
		for( i = 0; i < rank; i++ )	{
			n_i=hEle->nodes[i];
			for( j = 0; j < rank; j++ )	{
				n_j = hEle->nodes[j];
				nn = (i==j)?hEle->area/6.0:hEle->area/12.0;
				cu[n_j] += nn*u[n_i]*u_scale;
			}
		}	
	}
}

/*
	ISSUE-NEEDED:
		1、	为什么f[i]!=z[i]

	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/14/2008		
*/
void _FEM_form_f( GIP_FEM_ *hFEM,int flag )	{
	GIP_OBJECT *hObj=hFEM->hObj;
	int e,i,j,rank=3,no,dim=hFEM->K_dim,M=hObj->M,N=hObj->N;
	GIP_FEM_TRI* hEle;
	GIP_FLOATING *f=hFEM->f,*f_n=hFEM->D_temp,*u=hFEM->u;
	double nn,z_avg=0.0,*cu,u_scale=hFEM->u_scale,*data=hObj->data;

	memcpy( f_n,hFEM->u,sizeof(double)*M*N );
	switch( hFEM->model	)	{
	case GIP_SEG_MUMSHAH_s:
//	case GIP_DECOMP_Hs:
		_DECOM_Hs_force( M,N,N,f_n,hObj->data,0x0 );
		break;
	case GIP_MUMSHAH_filter:
		_DECOM_filer_force( M,N,N,f_n,hObj->data,0x0 );
		break;
	case GIP_MUMSHAH_1:
		for( i = 0; i < dim; i++ )	f_n[i]=data[i]==f_n[i] ? 1.0 : (data[i]-f_n[i])/fabs(data[i]-f_n[i]);
		break;
	case GIP_SEG_CV:
		_SEG_CV_force( M,N,N,f_n,hObj->data,0x0 );
		break;
	default:
		f_n = hObj->data;
		break;
	}
/*	if( BIT_TEST( hFEM->method,_FEM_U_SCALE ) )	{		//U_SPLIT似乎可取代time transient
		cu = hFEM->D_temp+dim;
		 _FEM_form_CU( hFEM,cu,0x0 );
		for( i = 0; i < dim; i++ )
			f_n[i] += cu[i];
	}*/
/*	
	for( i = 0; i < dim; i++ )	z_avg+=z[i];`	
	z_avg /= dim;
*/
	memset( f,0x0,sizeof(double)*dim );

	for( e = 0; e < hFEM->nEle; e++ )	{
		hEle = (GIP_FEM_TRI*)(hFEM->arrEle)+e;	
		for( i = 0; i < rank; i++ )	{
			no=hEle->nodes[i];
			for( j = 0; j < rank; j++ )	{
				nn = (i==j)?hEle->area/6.0:hEle->area/12.0;
//				f[no] += nn*z_avg;				//	!!testing!!			
				f[no] += nn*f_n[hEle->nodes[j]];		
				if( BIT_TEST( hFEM->method,_FEM_U_SCALE ) )
					f[no] += nn*u_scale*u[hEle->nodes[j]];
			}
		}		
	}
/*
	for( i = 0; i < dim; i++ )	{			//10/19/2008
		ASSERT( f[i]==z[i] );
	}
*/
}

/*
	copy from "Fundamentals of the Finite  Element Method for Heat and Fluid Flow" P.130
		r=0;r=M-1;c=0;c=N-1

	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/7/2008		
*/
void _FEM_boundary_condition_130( int r,int c,int M,int N,GIP_BOUNDARY_CONDITION *bc_type,double *bc_val )	{
	if( r==0 || r==M-1 || c==0 || c==N-1 )	{
		*bc_type=GIP_BC_FIX;
		*bc_val=r==M-1 ? 500 : 100.0;
	}
}

/*
	u_0有多个来源，可能来自hObj，也可能来自hMesh
	但长度始终为K_dim

	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/30/2008					
*/
void 	_FEM_init_u(  GIP_FEM_ *hFEM,GIP_FLOATING *u_0,int flag )	{
	int dim=hFEM->K_dim,i;
	GIP_FLOATING u_avg,*u=hFEM->u;
	
	switch( hFEM->model )	{
	case GIP_SEG_CV:
		u_avg=0.0;
		for( i = 0; i < dim; i++ )	u_avg +=u_0[i];
		u_avg /= dim;
		for( i = 0; i < dim; i++ )	
			u[i] = u_0[i]>u_avg ? 1.0 : -1.0;
		break;
	default:
		for( i = 0; i < dim; i++ )	u[i]=u_0[i];
		break;
	}
}

/*
	按照格点划分

	Copyright 2008-present, Grusoft.
	v0.1	9/26/2008		
*/
void _FEM_mesh_generate_0( GIP_FEM_*hFEM,int M,int N,int ldu,GIP_FEM_METHOD method,GIP_FLOATING *u )	{
	int r,c,pos,nTri=0,bc_kind=hFEM->bc_kind;
	GIP_FEM_NODE *hNode;
	GIP_FEM_TRI* hTri;

	hFEM->nNode=M*N;
	hFEM->K_dim=hFEM->nNode;
	hFEM->nEle=(M-1)*(N-1)*2;
	hFEM->arrNode=GIP_alloc( sizeof(GIP_FEM_NODE)*hFEM->nNode );
	hFEM->arrEle=GIP_alloc( sizeof(GIP_FEM_TRI)*hFEM->nEle );
	memset( hFEM->arrEle,0x0,sizeof(GIP_FEM_TRI)*hFEM->nEle );
	hFEM->I_temp=GIP_alloc( sizeof(int)*hFEM->nNode*2 );
	hFEM->D_temp=GIP_alloc( sizeof(double)*hFEM->K_dim*2 );
	hFEM->u=GIP_alloc( sizeof(GIP_FLOATING)*M*N );
	hFEM->f=GIP_alloc( sizeof(GIP_FLOATING)*M*N );
	memset( hFEM->f,0x0,sizeof(GIP_FLOATING)*M*N );
//boundary condition
	hFEM->bc_type=GIP_alloc( sizeof(GIP_BOUNDARY_CONDITION)*hFEM->K_dim );
	memset( hFEM->bc_type,0x0,sizeof(GIP_BOUNDARY_CONDITION)*hFEM->K_dim );
	hFEM->bc_val=GIP_alloc( sizeof(double)*hFEM->K_dim );

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			pos = GIP_RC2POS( r,c,N ); 
//			hFEM->u[pos] = u[GIP_RC2POS( r,c,ldu )];
			hNode = hFEM->arrNode+pos;
			hNode->r=r;		hNode->c=c;
			if( bc_kind==BOUNDARY_CONDITION_130 )
				_FEM_boundary_condition_130( r,c,M,N,hFEM->bc_type+pos,hFEM->bc_val+pos );
		}
	}
	for( r = 0; r < M-1; r++ )	{		
		for( c = 0; c < N-1; c++ )	{
			hTri = (GIP_FEM_TRI*)(hFEM->arrEle)+nTri;		nTri++;
			_FEM_tri_init( hTri,hFEM,GIP_RC2POS(r,c,N),GIP_RC2POS(r,c+1,N),GIP_RC2POS(r+1,c,N) );
			hTri = (GIP_FEM_TRI*)(hFEM->arrEle)+nTri;		nTri++;
			_FEM_tri_init( hTri,hFEM,GIP_RC2POS(r+1,c+1,N),GIP_RC2POS(r+1,c,N),GIP_RC2POS(r,c+1,N) );			
		}
	}
	ASSERT( hFEM->nEle == nTri);

	_FEM_init_u( hFEM,u,0x0 );

}

/*
	将nodal的值插到整个区域

	Copyright 2008-present, Grusoft

	v0.1	cys
		10/22/2008	
*/	
void _FEM_tri2domain( GIP_FEM_*hFEM,GIP_MESH *hMesh,int flag )	{
	int e,i,rank=3,no_1,no_2,no_3,r,c,r_0,r_1,c_0,c_1,cur;
	int M=hFEM->hObj->M,N=hFEM->hObj->N;
	GIP_FEM_TRI* hEle;
	GIP_FEM_NODE* nodes=hFEM->arrNode;
	GIP_FLOATING *Na,*Nb,*Nc,*u=hFEM->u,u1,u2,u3,*domain;
	double n1,n2,n3,a;

	domain=GIP_alloc( sizeof(GIP_FLOATING)*M*N );
	for( i = 0; i < M*N; i++ )		domain[i]=-1;

	u=hFEM->u;
	for( e = 0; e < hFEM->nEle; e++ )	{
		hEle = (GIP_FEM_TRI*)(hFEM->arrEle)+e;	
		r_0=INT_MAX,	r_1=INT_MIN;
		c_0=INT_MAX,	c_1=INT_MIN;
		for( i = 0; i < rank; i++ )	{
			cur=hEle->nodes[i];
			r=nodes[cur].r;			c=nodes[cur].c;
			r_0=MIN( r,r_0 );		r_1=MAX( r,r_1 );
			c_0=MIN( c,c_0 );		c_1=MAX( c,c_1 );
			domain[r*N+c] = u[cur];
		}
		if( hEle->area<=0.5 )		
			continue;
		if( r_1-r_0<=1 && c_1-c_0<=1 )		
			continue;
		no_1=hEle->nodes[0];	no_2=hEle->nodes[1];		no_3=hEle->nodes[2];	
		u1=u[no_1];				u2=u[no_2];					u3=u[no_3];	
		Na=hEle->N;			Nb=hEle->N+rank;			Nc=hEle->N+rank*2;
		a=hEle->area*2;
		for( r=r_0; r <= r_1; r++ )	{
		for( c=c_0; c <= c_1; c++ )	{
			n1=fabs(Na[0]+Nb[0]*c+Nc[0]*r);
			n2=fabs(Na[1]+Nb[1]*c+Nc[1]*r);
			n3=fabs(Na[2]+Nb[2]*c+Nc[2]*r);
			if( n1+n2+n3 != a )
				continue;
			if( n1==0 || n2==0 || n3== 0 )	{
				domain[r*N+c] = (n1*u1+n2*u2+n3*u3)/a;
			}else	{
				ASSERT( domain[r*N+c]==-1 );
				domain[r*N+c] = (n1*u1+n2*u2+n3*u3)/a;
			}
		}
		}
	}

	GIP_save_IplImage_d( M,N,GIP_path_(_T("mesh_domain"),hFEM->nNode),domain,N,0x0 );
	GIP_free( domain );
}

/*
	按照格点划分

	Copyright 2008-present, Grusoft.
	v0.1	10/22/2008		
*/
void _FEM_mesh_generate_( GIP_FEM_*hFEM,int M,int N,int ldu,GIP_FEM_METHOD method,GIP_FLOATING *u )	{
	int i,pos,nTri=0,bc_kind=hFEM->bc_kind,nNode,n_max;
	GIP_FEM_NODE *hNode;
	GIP_FEM_TRI* hTri;
	GIP_MESH mesh;

	n_max=MIN(M*N,M*N/100);
	GIP_MESH_gen( M,N,ldu,n_max,u,&mesh,0x0 );
	hFEM->nNode=mesh.nNode;				nNode=mesh.nNode;
	hFEM->K_dim=hFEM->nNode;
	hFEM->nEle=mesh.nEle;
	hFEM->arrNode=GIP_alloc( sizeof(GIP_FEM_NODE)*hFEM->nNode );
	hFEM->arrEle=GIP_alloc( sizeof(GIP_FEM_TRI)*hFEM->nEle );
	memset( hFEM->arrEle,0x0,sizeof(GIP_FEM_TRI)*hFEM->nEle );
	hFEM->I_temp=GIP_alloc( sizeof(int)*hFEM->nNode*2 );
	hFEM->D_temp=GIP_alloc( sizeof(double)*hFEM->K_dim*2 );
	hFEM->u=GIP_alloc( sizeof(GIP_FLOATING)*hFEM->K_dim );
	hFEM->f=GIP_alloc( sizeof(GIP_FLOATING)*hFEM->K_dim );
	memset( hFEM->f,0x0,sizeof(GIP_FLOATING)*hFEM->K_dim );
//boundary condition
	hFEM->bc_type=GIP_alloc( sizeof(GIP_BOUNDARY_CONDITION)*hFEM->K_dim );
	memset( hFEM->bc_type,0x0,sizeof(GIP_BOUNDARY_CONDITION)*hFEM->K_dim );
	hFEM->bc_val=GIP_alloc( sizeof(double)*hFEM->K_dim );

	for( i = 0; i < nNode; i++ )	{		
		hFEM->u[i] = mesh.n_val[i];
		hNode = hFEM->arrNode+i;
		hNode->r=mesh.n_cord[i];		hNode->c=mesh.n_cord[i+nNode];
	}
	for( i = 0; i < mesh.nEle; i++ )	{		
		hTri = (GIP_FEM_TRI*)(hFEM->arrEle)+nTri;		
		_FEM_tri_init( hTri,hFEM,mesh.tri_list[3*i],mesh.tri_list[3*i+1],mesh.tri_list[3*i+2] );
		nTri++;
	}	

//	_FEM_tri2domain( hFEM,&mesh,0x0 );

}

/*
	生成graph结构
	
	注意：
		1、总是hFEM->K_nz=hFEM->nNode

	Copyright 2008-present, Grusoft.
	v0.1	9/26/2008		
*/
void _FEM_graph_generate_( GIP_FEM_*hFEM,int flag )	{
	int *ptr,*ind,dim,nz_max,rank,*temp=hFEM->I_temp,*g_map;
	int i,j,k,e,n_i,n_j,isFind,nz;
	GIP_FEM_TRI* hEle;

	dim=hFEM->nNode;		nz_max=0;
	hFEM->K_nz=hFEM->nNode;
	for( i = 0; i < dim; i++ )		temp[i]=0;
	for( i = 0; i < hFEM->nEle; i++ )	{
		hEle = (GIP_FEM_TRI*)(hFEM->arrEle)+i;	
		rank = 3;
		for( j = 0; j < rank; j++ )	{
			temp[hEle->nodes[j]]+=rank;
			nz_max += rank;
		}
	}
	ptr=GIP_alloc( sizeof(int)*(dim+1) );		hFEM->K_ptr=ptr;
	ind=GIP_alloc( sizeof(int)*nz_max );		hFEM->K_ind=ind;
	ptr[dim] = nz_max;
	for( i = dim-1; i >= 0; i-- )	{
		ptr[i]=ptr[i+1]-temp[i];
		temp[i]=ptr[i];
	}
	ASSERT( ptr[0]==0 );

	for( e = 0; e < hFEM->nEle; e++ )	{
		hEle = (GIP_FEM_TRI*)(hFEM->arrEle)+e;	
		rank = 3;
		for( i = 0; i < rank; i++ )	{
			n_i=hEle->nodes[i];
			for( j = 0; j < rank; j++ )	{
				n_j = hEle->nodes[j];
				isFind = 0;
				for( k = ptr[n_i]; k < temp[n_i]; k++ )	{
					if( ind[k]==n_j )
					{	isFind=1;		break;	}
				}
				if( isFind==1 )			continue;
				ind[temp[n_i]]=n_j;		temp[n_i]++;	
			}
		}
	}
//整理
	nz = 0;
	for( i = 0; i < dim; i++ )	{
		ASSERT( temp[i]<=ptr[i+1] && nz<=ptr[i] );
		for( j = ptr[i]; j < temp[i]; j++ )	{
		for( k = j+1; k < temp[i]; k++ )	{
			if( ind[k] < ind[j] )	{
				e=ind[k];	ind[k]=ind[j];		ind[j]=e;
			}
		}
		}
		if( nz<ptr[i] )	{
			for( j = 0; j < temp[i]-ptr[i]; j++ )	{
				ind[nz+j]=ind[ptr[i]+j];
			}
		}
		nz += temp[i]-ptr[i];
		temp[i]=nz;
	}
	for( i = 0; i < dim; i++ )	{
		ptr[i+1] = temp[i];
	}
	ASSERT(	ptr[dim]==nz );
	hFEM->K_nz = nz;

	hFEM->K_val=GIP_alloc( sizeof(double)*nz );	

//get map
	for( e = 0; e < hFEM->nEle; e++ )	{
		hEle = (GIP_FEM_TRI*)(hFEM->arrEle)+e;	
		rank = 3;
		g_map = hEle->g_map;
		for( i = 0; i < rank; i++ )	{
			n_i=hEle->nodes[i];
			for( j = 0; j < rank; j++ )	{
				n_j = hEle->nodes[j];
				isFind = -1;
				for( k = ptr[n_i]; k < ptr[n_i+1]; k++ )	{
					if( ind[k]==n_j )
					{	isFind=k;		break;	}
				}
				if( isFind==-1 )
				{	isFind=-1;		}
				ASSERT( isFind!=-1 );			
				g_map[i*rank+j]=isFind;	
			}
		}
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	9/26/2008		
*/
void GIP_FEM_clear( GIP_FEM_ *hFEM )	{
	hFEM->nNode=-1;
	hFEM->nEle=-1;
	if( hFEM->hSS!=GIP_NULL )	{	
		GIP_SS_clear_( hFEM->hSS );	
		GIP_free( hFEM->hSS );
		hFEM->hSS=GIP_NULL;			
	}
	if( hFEM->arrNode!=GIP_NULL )
	{	GIP_free(hFEM->arrNode);		hFEM->arrNode=GIP_NULL;		}
	if( hFEM->arrEle!=GIP_NULL )
	{	GIP_free(hFEM->arrEle);			hFEM->arrEle=GIP_NULL;		}
	if( hFEM->I_temp!=GIP_NULL )
	{	GIP_free(hFEM->I_temp);			hFEM->I_temp=GIP_NULL;		}
	if( hFEM->D_temp!=GIP_NULL )
	{	GIP_free(hFEM->D_temp);			hFEM->D_temp=GIP_NULL;		}
	if( hFEM->K_ptr!=GIP_NULL )
	{	GIP_free(hFEM->K_ptr);			hFEM->K_ptr=GIP_NULL;		}
	if( hFEM->K_ind!=GIP_NULL )
	{	GIP_free(hFEM->K_ind);			hFEM->K_ind=GIP_NULL;		}
	if( hFEM->K_val!=GIP_NULL )
	{	GIP_free(hFEM->K_val);			hFEM->K_val=GIP_NULL;		}
	if( hFEM->u!=GIP_NULL )
	{	GIP_free(hFEM->u);				hFEM->u=GIP_NULL;			}
	if( hFEM->f!=GIP_NULL )
	{	GIP_free(hFEM->f);				hFEM->f=GIP_NULL;			}
	if( hFEM->bc_type!=GIP_NULL )
	{	GIP_free(hFEM->bc_type);		hFEM->bc_type=GIP_NULL;		}
	if( hFEM->bc_val!=GIP_NULL )
	{	GIP_free(hFEM->bc_val);			hFEM->bc_val=GIP_NULL;		}
}

/*
	The matrix A is symmetric(epsi精度校验)

	v0.1	cys
		10/8/2008	
*/
void _FEM_verify_MAT( int dim,int *ptr,int *ind,double *val,double epsi )	{
	int isFind,i,j,k,r;

	for( i = 0; i < dim; i++ )	{
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{			
			r = ind[j];
			if( j > ptr[i] )		ASSERT( r > ind[j-1] );
//			if( r > i )	{
				isFind = -1;
				for( k = ptr[r]; k < ptr[r+1];	k++ )	{
					if( ind[k]==i )	{
						isFind=k;				
						break;
					}
				}
				if( isFind==-1 )	
				{	i=0;	}
				ASSERT( isFind!=-1 );
				ASSERT( fabs(val[j]-val[k])<=epsi*fabs(val[j]) );
//			}
		}
	}
	
}


/*
[**]
	在梦里，我又回到和桥，温暖的气息漂浮在粼粼的波光里。。。	
	于是，我的灵魂得到安息而沉沉睡去。
[**]

	装配hFEM->K,hFEM->F，同时形成hSS->mat,hSS->rhs

	目标：
		1 linear(装配K),nonlinear(还需装配J)
		2 稳态,瞬态（C修正)
		3 Fix_Point,Newton
		4 单变量，多变量(mixed,vector)
	

	Copyright 2008-present, Grusoft.
	v0.1	cys
		9/27/2008		
*/
void _FEM_assemb( GIP_FEM_ *hFEM,int flag )	{
	GIP_SPARSE_SOLVER* hSS=hFEM->hSS;
	int i,j,rank,n_i,n_j,pos_e,pos,e,dim=hSS->dim,model=hFEM->model;
	int isTransient=0,isAssembJ=0/*,isUpdateRhs=0*/;
	double *val,*rhs,xita=hFEM->xita,*J,dT=hFEM->dT;
	GIP_FLOATING *C,*K,*u=hFEM->u,*f=hFEM->f,*K_val=hFEM->K_val;
	GIP_FEM_TRI* hEle;
	GIP_OBJECT *hObj=hFEM->hObj;

	isTransient=BIT_TEST( hFEM->method,GIP_FEM_TRANSIENT );
	isAssembJ=BIT_TEST( hFEM->method,_FEM_NONLINEAR_NEWTON );
//	isUpdateRhs=isAssembJ || model==GIP_SEG_CV || model==GIP_DECOMP_Hs;

	val = hSS->val;			rhs=hSS->rhs;
	memset( rhs,0x0,sizeof(double)*hSS->dim );
	memset( K_val,0x0,sizeof(double)*hSS->nz );
//	memset( rhs,0x0,sizeof(double)*hSS->dim );
	if( isAssembJ )	{
		memset( val,0x0,sizeof(double)*hSS->nz );
		J=GIP_alloc( sizeof(double)*_FEM_RANK_MAX*_FEM_RANK_MAX );
	}

	if( isTransient )	{
		for( i = 0; i < dim; i++ )		rhs[i]=dT*f[i];
	}else
		for( i = 0; i < dim; i++ )		rhs[i]=hFEM->f[i];

	for( e = 0; e < hFEM->nEle; e++ )	{
		hEle = (GIP_FEM_TRI*)(hFEM->arrEle)+e;	
		rank = 3;
		C=hEle->C;		K=hEle->K;		
		for( i = 0; i < rank; i++ )	{
			n_i=hEle->nodes[i];
			for( j = 0; j < rank; j++ )	{
				n_j = hEle->nodes[j];
				pos_e=i*rank+j;
				pos = hEle->g_map[pos_e];
				if( isTransient )	{
					rhs[n_j] += (C[pos_e]-(1-xita)*dT*K[pos_e])*u[n_i];
					K_val[pos] += C[pos_e]+xita*dT*K[pos_e];
				}else
					K_val[pos] += K[pos_e]; 
			}
		}	
		val = hSS->val;
		if( isAssembJ ) {
			_FEM_form_J( hFEM,e,J,0x0 );
			for( i = 0; i < rank; i++ )	{
				n_i=hEle->nodes[i];
				for( j = 0; j < rank; j++ )	{
					n_j = hEle->nodes[j];
					pos_e=i*rank+j;
					pos = hEle->g_map[pos_e];
					val[pos] += J[pos_e]; 
				}
			}	
		}
	}
	if( isAssembJ )	{
		for( i=0; i < dim; i++ )	{
			for( j = hSS->ptr[i]; j < hSS->ptr[i+1]; j++ )	{
				rhs[hSS->ind[j]] -= K_val[j]*u[i];
			}
		}
	}
	if( !isAssembJ )
		memcpy( val,K_val,sizeof(double)*hSS->nz );
	if( isAssembJ )			GIP_free( J );

//	_FEM_verify_MAT( hSS->dim,hSS->ptr,hSS->ind,hSS->val );
	if( hFEM->bc_kind !=BOUNDARY_CONDITION_NONE )		//boundary condition
		GIP_SS_update_( hSS,hFEM->bc_type,hFEM->bc_val );
	_FEM_verify_MAT( hSS->dim,hSS->ptr,hSS->ind,hSS->val,FLT_EPSILON );
	if( 0 )	{		//	!!testing!!		
		char sPath[80];
		_FEM_verify_MAT( hFEM->K_dim,hFEM->K_ptr,hFEM->K_ind,hFEM->K_val,DBL_EPSILON );
		_stprintf( sPath,"..\\data\\FEM_K_%d_%d_0.mtx",hFEM->K_dim,hFEM->loop );
		Equation_output( hObj->M,hObj->N, hFEM->K_dim,hFEM->K_ptr,hFEM->K_ind,hFEM->K_val,hFEM->f,sPath );
		_stprintf( sPath,"..\\data\\SS_K_%d_%d_0.mtx",hFEM->K_dim,hFEM->loop );
		Equation_output( hObj->M,hObj->N, hSS->dim,hSS->ptr,hSS->ind,hSS->val,hSS->rhs,sPath );
//		_stprintf( sPath,"..\\data\\FEM_J_%d_%d.mtx",hSS->dim,hFEM->loop );
//		Equation_output( hObj->M,hObj->N,hSS->dim,hSS->ptr,hSS->ind,hSS->val,hSS->rhs,sPath );
	}
}

/*
[*]	告别幼稚与幻想，走向成熟与现实。。。	

	Copyright 2008-present, Grusoft.
	v0.1	9/26/2008
		void GIP_FEM_preprocess( GIP_FEM_*hFEM,int M,int N,int ldu,GIP_FEM_METHOD method,GIP_FLOATING *u,int model )	
	v0.2	10/28/2008
		引用GIP_OBJECT
		
*/
void GIP_FEM_preprocess( GIP_FEM_*hFEM,GIP_OBJECT *hObj,int model,GIP_FEM_METHOD method )	{
	int dim,nz_max=0,m_type=0;
	double f_norm=0.0;

	memset( hFEM,0x0,sizeof(GIP_FEM_) );
	hFEM->hObj=hObj;
//	hFEM->M=M;			hFEM->N=N;			
	hFEM->xita=1.0;			//Fully implicit scheme; Backward Time difference 
	hFEM->model=model;
	hFEM->u_scale = 1.0;
	BIT_SET( hFEM->method,_FEM_ERR_ENERGY);
//Y.Meyer Lemma: If |f|* dose not exceed 0.5/lenda, then the optimal ROF decomposition is: u=0, v=f.
//	ASSERT(  f_norm > _DU_ALPHA/2.0 );			//from Y.Meyer Theorem.

	switch( model )	{
	case GIP_MODEL_130:
		hFEM->bc_kind=BOUNDARY_CONDITION_130;
		break;
	case GIP_SEG_MUMSHAH_s:
	case GIP_MUMSHAH_1:
	case GIP_MUMSHAH_filter:
		hFEM->bc_kind=BOUNDARY_CONDITION_NONE;
		BIT_SET( hFEM->method,GIP_FEM_NONLINEAR );
		BIT_SET( hFEM->method,_FEM_NONLINEAR_FIXPOINT );
//		BIT_SET( hFEM->method,_FEM_NONLINEAR_NEWTON );
		BIT_SET( hFEM->method,_FEM_U_SCALE );
		hFEM->u_scale = 1;
		break;
	case GIP_SEG_MUMSHAH_2:
	case GIP_MODEL_CGM:
		hFEM->bc_kind=BOUNDARY_CONDITION_NONE;
		BIT_RESET( hFEM->method,_FEM_ERR_ENERGY);		BIT_SET( hFEM->method,_FEM_ERR_RES);
//		BIT_SET( hFEM->method,_FEM_NONLINEAR_NEWTON );
		BIT_SET( hFEM->method,GIP_FEM_NONLINEAR );
		BIT_SET( hFEM->method,_FEM_NONLINEAR_FIXPOINT );
		hFEM->u_scale = 0.0;
		break;
	case GIP_SEG_CV:
		hFEM->bc_kind=BOUNDARY_CONDITION_NONE;
//		BIT_SET( hFEM->method,_FEM_NONLINEAR_NEWTON );
		BIT_SET( hFEM->method,_FEM_NONLINEAR_FIXPOINT );
		BIT_SET( hFEM->method,GIP_FEM_STEADY );
		BIT_SET( hFEM->method,GIP_FEM_NONLINEAR );
		break;
	default:
		break;
	}
	if( !BIT_TEST( hFEM->method,_FEM_U_SCALE ) )
		ASSERT( hFEM->u_scale==0.0);

	if( BIT_TEST( hFEM->method,GIP_FEM_TRANSIENT ) )
		hFEM->dT = 1.0;
	_FEM_mesh_generate_0( hFEM,hObj->M,hObj->N,hObj->N,method,hObj->data );
//	_FEM_mesh_generate_( hFEM,hObj->M,hObj->N,hObj->N,method,hObj->data );
	_FEM_form_f( hFEM,0x0 );

	_FEM_graph_generate_( hFEM,0x0 );
	dim=hFEM->K_dim;			nz_max=hFEM->K_nz;
	hFEM->hSS=GIP_SS_init_( dim,nz_max,hFEM->K_ptr,hFEM->K_ind,GIP_NULL,m_type,0x0  );

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/11/2008		
*/
void GIP_FEM_dump( GIP_FEM_*hFEM,int flag )	{
	int method=hFEM->method;
	char sModel[80];

	GIP_MODEL_info( sModel,hFEM->model,0x0 );

	G_PRINTF( "FEM:\tnode=%d,ele=%d,method=0x%x\r\n",hFEM->nNode,hFEM->nEle,method );	
	if( BIT_SET( method,GIP_FEM_NONLINEAR ) )		{	
		G_PRINTF( "\tNONLINEAR:\t%s\t%s",BIT_TEST( method,_FEM_NONLINEAR_FIXPOINT )?"FIXPOINT":
			BIT_TEST( method,_FEM_NONLINEAR_NEWTON )?"NEWTON":"!!!!",BIT_TEST( hFEM->method,_FEM_ERR_ENERGY)?
			"CHECK ENERGY":"CHECK RES" );
		G_PRINTF( "\r\n" );
	}	
	if( BIT_TEST( method,_FEM_U_SCALE ) )	
	{	G_PRINTF( "\t<scale=%g>",hFEM->u_scale );		}

	G_PRINTF( "\t<alpha=%g> <beta=%g> <ksi=%g> <bd=%d>\r\n",(double)_DU_ALPHA,(double)_DU2_BETA,(double)_DU_1_KSI,
		hFEM->bc_kind );
}

/*
	收敛指标：对比迭代过程的图像，似乎泛函的energy较好，而偏微分方程的余量有些问题。
	可能是FEM的线性模型与fourier周期模型不匹配，也可能偏微分方程更敏感（非线性更强）

	Copyright 2008-present, Grusoft.
	v0.1	cys
		10/19/2008		
*/
double _FEM_res( GIP_FEM_*hFEM,double *u_old,double resi_0,int flag )	{
	int dim=hFEM->K_dim,*ptr=hFEM->K_ptr,*ind=hFEM->K_ind,i,j,r,no_max;
	double *res=hFEM->D_temp,*f=hFEM->f,*K_val=hFEM->K_val,*u=hFEM->u;
	double a,r_norm=0.0,off=0.0,off_max,u_1,u_0,off_s;
	GIP_SPARSE_SOLVER* hSS=hFEM->hSS;
	GIP_OBJECT *hObj=hFEM->hObj;

	if( BIT_TEST( hFEM->method,_FEM_ERR_ENERGY) )	{
		u_1=_FEM_DU( hFEM,0x0 );
		if( hFEM->model==GIP_SEG_MUMSHAH_s || hFEM->model==GIP_MUMSHAH_filter )
			_DECOM_Hs_energy( hObj->M,hObj->N,hObj->N,u,hObj->data,&off_s,0x0 );
		else if( hFEM->model==GIP_SEG_MUMSHAH_2 )	{
			off_s = 0.0;
			for( i=0; i < dim; i++ )	off_s+=(hObj->data[i]-u[i])*(hObj->data[i]-u[i]);
		}	else if( hFEM->model==GIP_MUMSHAH_1 )	{
			off_s = 0.0;
			for( i=0; i < dim; i++ )	off_s+=fabs(hObj->data[i]-u[i]);
		}
		r_norm = _DU_ALPHA*u_1+off_s;
	}else	{
		for( i=0; i < dim; i++ )	res[i]=f[i];
		for( i=0; i < dim; i++ )	{
			for( j = ptr[i]; j < ptr[i+1]; j++ )	{
				r = ind[j];
				res[r] -= K_val[j]*u[i];
			}
		}
		for( i=0; i < dim; i++ )	{
	//		ASSERT( res[i]==hSS->rhs[i] );
			r_norm += res[i]*res[i];
		}
		r_norm=sqrt(r_norm);
	}
	r_norm/=dim;
	if( u_old != GIP_NULL )		{
		off = 0.0;			off_max=0.0;		no_max=-1;
		u_1=DBL_MIN;		u_0=DBL_MAX;
		for( i=0; i < dim; i++ )	{
			a = fabs( u_old[i]-u[i] );
			if( u[i] > u_1 )
			{	u_1=u[i]; }
			if( u[i] < u_0 )
			{	u_0=u[i]; }

			off += a*a;
			if( a > off_max )
			{	no_max=i;	off_max = a;	}
		}
		off=sqrt(off);
	}	

	if( hFEM->loop>0 )	{
		ASSERT( u_old != GIP_NULL );
		G_PRINTF( "\t*  %2d:res=%.3lf,|dU|=%.1lf(%.3g)\r\n",hFEM->loop,r_norm,r_norm/resi_0,off,off_max );
		if( 0 )
			G_PRINTF( "\t*  u_1=%g,u_0=%g\r\n",u_1,u_0 );
	}else
		G_PRINTF( "\t*  %2d:res=%.3lf(%g)\r\n",hFEM->loop,r_norm,r_norm*dim );

	return r_norm;
	
}

/*
	steady Analysis

	Copyright 2008-present, Grusoft.
	v0.1	cys
		9/26/2008		
*/
void GIP_FEM_core( GIP_FEM_*hFEM,int loop_max,int flag )	{
	int dim,i,isTransient;
	GIP_SPARSE_SOLVER* hSS=hFEM->hSS;
	GIP_FLOATING *u=hFEM->u;
	double *rhs=hSS->rhs,resi_0,res,*u_old;
	GIP_OBJECT *hObj=hFEM->hObj;
	
	isTransient=BIT_TEST( hFEM->method,GIP_FEM_TRANSIENT );
	//BIT_SET( hFEM->method,GIP_FEM_STEADY );

	dim=hFEM->nNode;
	u_old=GIP_alloc( sizeof(double)*dim*2 );
	_FEM_form_K( hFEM,0x0 );
	_FEM_assemb( hFEM,0x0 );
	_FEM_form_f( hFEM,0x0 );
	if( BIT_TEST(hFEM->method,GIP_FEM_NONLINEAR) || isTransient )	{
		hFEM->loop=0;
		resi_0 = _FEM_res( hFEM,GIP_NULL,0.0,0x0 );
		if( resi_0==0.0 )	{
			res=0.0;	goto QUIT;
		}
		do	{
			GIP_save_IplImage_d( hObj->M,hObj->N,GIP_path_(_T("u_fem"),hFEM->loop),hFEM->u,hObj->N,0x0 );
			hFEM->loop++;
			GIP_SS_solve_( hFEM->hSS,0x0 );
			for( i = 0; i < dim; i++ )	u_old[i]=u[i];
			if( BIT_TEST(hFEM->method,_FEM_NONLINEAR_NEWTON) )
				for( i = 0; i < dim; i++ )	u[i]+=rhs[i];
			else
				for( i = 0; i < dim; i++ )	u[i]=rhs[i];
			_FEM_form_K( hFEM,0x0 );
//			if( BIT_TEST( hFEM->method,_FEM_F_UPDATE )
			if( hFEM->model==GIP_SEG_CV || hFEM->model==GIP_MUMSHAH_filter  || hFEM->model==GIP_MUMSHAH_1 || hFEM->model==GIP_SEG_MUMSHAH_s )
				_FEM_form_f( hFEM,0x0 );
			_FEM_assemb( hFEM,0x0 );
			res=_FEM_res( hFEM,u_old,resi_0,0x0 );
		}while( res>resi_0*1.0e-4 && hFEM->loop<loop_max );
	}else	{
		GIP_SS_solve_( hFEM->hSS,0x0 );
		for( i = 0; i < dim; i++ )	u[i]=rhs[i];
	}	
QUIT:
	GIP_free( u_old );
	G_PRINTF( "\t******  res_0=%g,res=%g,loop=%d\r\n",resi_0,res,hFEM->loop );	

}

/*
	Transient Analysis

	Copyright 2008-present, Grusoft.
	v0.1	cys
		9/26/2008		
*/
void GIP_FEM_transient( GIP_FEM_*hFEM,double dT )	{
	int dim,i;
	GIP_SPARSE_SOLVER* hSS=hFEM->hSS;
	GIP_FLOATING *u=hFEM->u;
	double *rhs=hSS->rhs,dU;
	
	BIT_SET( hFEM->method,GIP_FEM_TRANSIENT );

	dim=hFEM->nNode;
	_FEM_form_K( hFEM,0x0 );
	_FEM_assemb( hFEM,0x0 );
	GIP_SS_solve_( hFEM->hSS,0x0 );

	dU=0.0;
	for( i = 0; i < dim; i++ )	{
		dU += (u[i]-rhs[i])*(u[i]-rhs[i]);
		u[i]=rhs[i];
	}
	hFEM->dU = sqrt(dU/dim);
 }


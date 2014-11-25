#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include <math.h>
#include "../PCH/GIP_def.h"
#include "../util/gip_heap.h"
#include "../util/gip_util.h"
#include "../util/Compact_DataStruc.h"
#include "../util/gip_convexhull.h"
#include "gip_mesh.h"

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		10/22/2008	
*/	
int _MESH_node_order( int M,int N,int ldu,int n_max,GIP_FLOATING *u,int *nodes )	{
	int r,c,nFix=0,nNode=M*N,no,reduce=0,pos,pact,nz=0,cur,r_step,c_step,isFix;
	double off,a,ratio=n_max*1.0/M/N;
	PACT_order_list POL;

	if( nNode<=n_max )	{
		for( no = 0; no < M*N; no++ )	nodes[no]=no;
		return nNode;
	}

	pact=255*4+1;
	PACT_ol_init( &POL,pact,MAX(nNode,M*ldu),0x0 );	
	r_step=1.0/ratio;			c_step=1.0/ratio;

	for( r=0; r<M; r++ )	{
	for( c=0; c<N; c++ )	{
		isFix=0;
		cur=GIP_RC2POS( r,c,ldu );
		if( r==0 || r==M-1 )	//boundary and other condition
		{	isFix = c==0||c==N-1||c%c_step==0;		}
		if( c==0 || c==N-1 )		
		{	isFix = r==0||r==M-1||r%r_step==0;		}
		if( isFix==1 )	{
			nodes[nFix++] = cur;
			continue;		
		}
		off = 0.0;	
		a = u[cur];
		if( r>0 )	{
			pos=GIP_RC2POS( r-1,c,ldu );
			off += fabs(a-u[pos]);
		}
		if( r<M-1 )	{
			pos=GIP_RC2POS( r+1,c,ldu );
			off += fabs(a-u[pos]);
		}		
		if( c>0 )	{
			pos=GIP_RC2POS( r,c-1,ldu );
			off += fabs(a-u[pos]);
		}
		if( c<N-1 )	{
			pos=GIP_RC2POS( r,c+1,ldu );
			off += fabs(a-u[pos]);
		}		
		PACT_ol_insert( &POL,cur,(int)off );		
	}
	}

	nz=nFix;
	reduce=nNode-n_max+1;
	while( (cur=PACT_ol_minpop(&POL) ) != -1 )	{
		if( reduce>0 )
		{	reduce--;	continue;		}
		nodes[nz++] = cur;	
	}

	PACT_ol_clear( &POL );

	return nz;
}

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		10/22/2008	
*/
void _mesh_ouput( GIP_MESH *hMesh,int flag )	{
	int i,j,k,rank,M=hMesh->M,N=hMesh->N,*tri=hMesh->tri_list,nNode=hMesh->nNode;
	int n_1,n_2,r_1,c_1,r_2,c_2,r_0,c_0,nSec;
	double *u,r,c;

	u=GIP_alloc( sizeof(double)*M*N );
	for( i = 0; i < M*N; i++ )	{
		u[i]=255;
	}
	rank=3;
	for( i = 0; i < hMesh->nEle; i++ )	{
		for( j = 0; j < rank; j++ )	{
			n_1=tri[3*i+j];		n_2=j==rank-1 ? tri[3*i] : tri[3*i+j+1] ;		
			r_1=hMesh->n_cord[n_1];			c_1=hMesh->n_cord[n_1+nNode];
			r_2=hMesh->n_cord[n_2];			c_2=hMesh->n_cord[n_2+nNode];
			u[r_1*N+c_1]=0.0;
			nSec=MAX( fabs(r_2-r_1),fabs(c_1-c_2) );
			if( nSec <= 1 )	continue;
			for( k = 1; k < nSec; k++ )	{
				r=r_1+1.0*(r_2-r_1)*k/nSec;		r_0=(int)(r+0.5);
				c=c_1+1.0*(c_2-c_1)*k/nSec;		c_0=(int)(c+0.5);
//				u[r_0*N+c_0]=0.0;				
			}
		}
	}

	GIP_save_IplImage_d( M,N,GIP_path_(_T("mesh_tri"),hMesh->nNode),u,N,0x0 );
	GIP_free( u );
}

/*
[**]
	泪水滴下，洗去尘封，现出往日的笑容...
[**]

	Copyright 2008-present, Grusoft

	v0.1	cys
		10/22/2008	
*/
void GIP_MESH_gen( int M,int N,int ldu,int n_max,GIP_FLOATING *u,GIP_MESH *hMesh,int flag )	{
	int i,r,c,*temp,*nodes,nz,ne_dt;
	double *dt_cord,*n_val;

	hMesh->M=M;		hMesh->N=N;
	hMesh->n_max=n_max;

	temp=GIP_alloc( sizeof(int)*M*N ),		nodes=temp;
	nz = _MESH_node_order( M,N,ldu,n_max,u,nodes );
	ASSERT( nz<=(n_max+(M+N)*2) && nz>=3 );

	dt_cord = GIP_alloc( sizeof(double)*2*nz );		hMesh->n_cord=dt_cord;
	n_val=GIP_alloc( sizeof(double)*nz );			hMesh->n_val=n_val;
	for( i = 0; i < nz; i++ )	{
	   r = GIP_POS2R(nodes[i],ldu);		c = GIP_POS2C(nodes[i],ldu);
	   hMesh->n_val[i]=u[nodes[i]];
	   dt_cord[i] = r;
	   dt_cord[i+nz] = c;
	}
//	Delau_Tri_1( nz,dt_cord,&ne_dt,&hMesh->tri_list );

	hMesh->nNode=nz;
	hMesh->nEle=ne_dt/3;

	GIP_free( temp );	

	_mesh_ouput( hMesh,0x0 );
//	hMesh->G_ptr=GIP_alloc(sizeof(int)*(hMesh->G_dim+1) );
}

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		10/22/2008	
*/	
void GIP_MESH_clear( GIP_MESH *hMesh )	{
	if( hMesh->n_cord!=GIP_NULL )
	{	GIP_free(hMesh->n_cord);		hMesh->n_cord=GIP_NULL;		}
	if( hMesh->tri_list!=GIP_NULL )
	{	GIP_free(hMesh->tri_list);		hMesh->tri_list=GIP_NULL;		}
}



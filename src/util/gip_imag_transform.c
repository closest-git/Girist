#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include <math.h>
#include <stdio.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "gip_util.h"
#include "gip_imag_transform.h"

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/3/2008	
*/
void GIP_Binarize( GIP_OBJECT *hSrc,GIP_OBJECT *hDest,double thresh,int flag )	{
	int r,c,pos,M=hSrc->M,N=hSrc->N;
	double a;
	ASSERT( hSrc->M==M && hSrc->N==N );

	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,N );
		a =hSrc->data[pos];
		hDest->data[pos] = a>thresh? 1.0 : 0.0;
	}
	}
}

/*
	注意:
		hDest已定义

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/3/2008	
*/
void GIP_OBJ_sub( GIP_OBJECT *hSrc,GIP_OBJECT *hDest,int r_0,int c_0,int dr,int dc,int flag )	{
	int r,c,pos,r_1=r_0+dr,c_1=c_0+dc;
	double a;

	ASSERT( hSrc!=hDest );
	ASSERT( r_1<=hSrc->M && c_1<=hSrc->N && r_0>=0 && c_0 >= 0);
	ASSERT( hDest->M==dr && hDest->N==dc );

	for( r = r_0; r < r_1; r++ )	{		
	for( c = c_0; c < c_1; c++ )	{
		pos = GIP_RC2POS( r,c,hSrc->N );
		a =hSrc->data[pos];
		pos = GIP_RC2POS( r-r_0,c-c_0,hDest->N );
		hDest->data[pos] = a;
	}
	}
}

/*
	注意:
		hDest已定义

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/6/2009	
*/
void GIP_OBJ_shift( GIP_OBJECT *hSrc,GIP_OBJECT *hDest,int dr,int dc,int flag )	{
	int r,c,pos,M=hSrc->M,N=hSrc->N;
	GIP_FLOATING a,*s_data=hSrc->data,*d_data=hDest->data;

	ASSERT( hSrc!=hDest );
	ASSERT( hDest->M==M && hDest->N==N );

	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( (r+dr+M)%M,(c+dc+N)%N,N );
		a =s_data[pos];
		pos = GIP_RC2POS( r,c,N );
		d_data[pos] = a;
	}
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/3/2008	
*/
void GIP_OBJ_init( GIP_OBJECT *hObj,int M,int N,int flag )	{
	hObj->M=M;			hObj->N=N;
	hObj->data=GIP_alloc( sizeof(GIP_FLOATING)*M*N );
}	

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/13/2008	
*/
GIP_OBJECT *GIP_OBJ_create( int M,int N,GIP_FLOATING *src,int flag )	{
	int i;

	GIP_OBJECT *hObj=GIP_alloc( sizeof(GIP_OBJECT) );		ASSERT( hObj != GIP_NULL  );
	hObj->M=M;			hObj->N=N;
	hObj->data=GIP_alloc( sizeof(GIP_FLOATING)*M*N );
	if( src != GIP_NULL )	{
		for( i = 0; i < M*N; i++ )	hObj->data[i]=src[i];
	}
	if( BIT_TEST( flag,GIP_OBJ_2ZERO ) )	{
		for( i = 0; i < M*N; i++ )	hObj->data[i]=0.0;
	}

	return hObj;
}	

/*
	扩充GIP_OBJ_create,可生成指定类型的mat

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/14/2008	
*/
GIP_OBJECT *GIP_MAT_create( int M,int N,GIP_FLOATING *src,int type,int flag )	{
	int i;
	GIP_OBJECT *hObj=GIP_alloc( sizeof(GIP_OBJECT) );
	hObj->M=M;			hObj->N=N;
	hObj->data=GIP_alloc( sizeof(GIP_FLOATING)*M*N );
	switch( type )	{
	default:
		if( src != GIP_NULL )	{
			for( i = 0; i < M*N; i++ )	hObj->data[i]=src[i];
		}
		break;
	}


	return hObj;
}	

/*
	Canny Edge Detector
	并返回g_I,g_rad(边上每个点的强度与方向)

	注意：
		1 暂不包括smooth
	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/17/2009			
*/
int GIP_canny_( GIP_OBJECT *hObj,int I_1,int I_0,double g_t0,double g_t1,int *mask,double *g_I,double *g_rad,int flag )	{
	int pos,M=hObj->M,N=hObj->N,r,c,MN=M*N,dgr,nz,adj[8],M_cut,e_nz;
	GIP_FLOATING *data=hObj->data;
	double gr,gc,w;

	ASSERT( I_0<I_1 && g_t0<g_t1 );
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos=GIP_RC2POS( r,c,N );
//		if( pos==2367 )
//			w=0;
		if( data[pos]>I_1 || data[pos]<I_0 )		
		{	mask[pos]=-1;			}
		if( BIT_TEST( flag,G_MASK_3 ) )	{
			gr = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_R,0x0 );
			gc = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_C,0x0 );
		}else	{
			gr = -GIP_Mask_5_( pos,M,N,N,data,G_MASK_FIVE_R,0x0 );
			gc = -GIP_Mask_5_( pos,M,N,N,data,G_MASK_FIVE_C,0x0 );
		}
		g_I[pos] = (sqrt( gr*gr+gc*gc ));
		g_rad[pos] =  (atan2( gr,gc ));
	}
	}

	e_nz=0;
	M_cut=M/10;
	for( r = M_cut; r < M-M_cut; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		if( mask[pos]<0 )continue;
		if( g_I[pos] < g_t0 )	
			continue;
		w = fabs(g_rad[pos])*57.295779513082;
		dgr = (22.5<w&&w<=67.5) ? 1 : (67.5<w&&w<=112.5) ? 2 : (112.5<w&&w<=157.5) ? 3 : 0;
		if( dgr==0 && (g_I[pos]<g_I[pos-1] || g_I[pos]<g_I[pos+1]) )		//nonmaximum suppression
			continue;
		if( dgr==3 && (g_I[pos]<g_I[pos-N-1] || g_I[pos]<g_I[pos+N+1]) )
			continue;
		if( dgr==2 && (g_I[pos]<g_I[pos-N] || g_I[pos]<g_I[pos+N]) )
			continue;
		if( dgr==1 && (g_I[pos]<g_I[pos-N+1] || g_I[pos]<g_I[pos+N-1]) )
			continue;
		if( g_I[pos] >= g_t1 )	{
			mask[pos] = 1;			e_nz++;				
		}else
			mask[pos] = 10;
	}
	}
/*
#ifdef _DEBUG
	I_ouput( M*N,mask,".\\trace\\2.dat" );
#else
	I_ouput( M*N,mask,".\\trace\\1.dat" );
#endif
*/
	do{
		nz = 0;			
		for( r = 1; r < M-1; r++ )	{		
		for( c = 1; c < N-1; c++ )	{
			pos=GIP_RC2POS( r,c,N );
			if( mask[pos]!=10 )		continue;
			mask[pos]=1;
			if( GIP_adjacent_2( M, N,N,mask,r,c,adj,GIP_ADJACENT_M )>0 )	{
				mask[pos]=1;		nz++;		
			}else
				mask[pos]=10;		
		}
		}
		e_nz += nz;
	}while( nz>0 );

	if( 0 )	{
		for( r = 1; r < M-1; r++ )	{		
		for( c = 1; c < N-1; c++ )	{
			pos=GIP_RC2POS( r,c,N );
			if( mask[pos]==1 )	
				data[pos]=255.0;
			else
				data[pos]=0.0;
		}
		}
//		D_ouput( M*N,data,".\\trace\\1.dat" );
//		GIP_save_IplImage_d( M,N,".\\trace\\canny_edge.jpg",data,N,0x0 );
	}

	return e_nz;
}

/*
	edge detector
	[ref]: SUPPORTING RANGE AND SEGMENT-BASED HYSTERESIS THRESHOLDING IN EDGE DETECTION
	
	效率 **-
		核心模块，must be faster,当图像较大时，效率较低
		M_cut需进一步确定
		依赖I_1,I_0,g_t0,g_t1多个参数 
	健壮 ***+			
		要优于GIP_canny_，对比见 2_22_CASIA_1,2].info
	稳定 ***+			
		g_rad之计算有误差，似乎影响不大。例如：
				|174 165 174|
				|173 168 173|
				|172 168 172|

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/20/2009			
	v0.2	cys			增加_45,_135方向的梯度计算
		3/27/2009			
*/
int GIP_lsu_( GIP_OBJECT *hObj,int I_1,int I_0,double g_t0,double g_t1,int *mask,double *g_I,double *g_rad,int flag )	{
	int M=hObj->M,N=hObj->N,MN=M*N,dgr,nz,adj[8],M_cut,e_nz,search;
	int i,k,r,c,r1,c1,pos,pos_1,*arr_10,nz_1=0,nz_2=0;
//	int neibor[4][2]={	{-1,1},{-N+1,N-1},{-N,N},{-N-1,N+1}	};
	int neibor[4][2]={	{-1,1},{-N-1,N+1},{-N,N},{-N+1,N-1}	};
	GIP_FLOATING *data=hObj->data;
	double gr,gc,g_1,g_2,rad_1,rad_2,*w,d[2],dr,dc,s,s_g,g_sum,support,I_off,rad_off;

	w = GIP_alloc( sizeof(double)*MN );
	support = 0.1*MIN(M,N);

	ASSERT( I_0<I_1 && g_t0<g_t1 );
	I_off=0.0,			rad_off=0.0;
	GIP_MEMCLEAR( g_I,sizeof(double)*MN );
	GIP_MEMCLEAR( g_rad,sizeof(double)*MN );
	arr_10=GIP_alloc( sizeof(int)*MN );		
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		if( pos==7937 )
			i=0;
		if( data[pos]>I_1 || data[pos]<I_0 )		
		{	mask[pos]=-1;			}
//		gr = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_R,0x0 );
//		gc = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_C,0x0 );
		gr = -(data[pos+N-1]+data[pos+N]*2+data[pos+N+1]-data[pos-N-1]-data[pos-N]*2-data[pos-N+1])/8;
		gc = -(data[pos-N+1]-data[pos-N-1]+data[pos+1]*2-data[pos-1]*2+data[pos+N+1]-data[pos+N-1])/8;
		g_1 = sqrt( gr*gr+gc*gc );
		rad_1 = (atan2( gr,gc ));
		if( 0 )	{
			g_I[pos] = g_1;
			g_rad[pos] = rad_1;
		}else	{
//			gr = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_135,0x0 );
//			gc = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_45,0x0 );
			gr = -(-data[pos-N]-data[pos-N+1]*2+data[pos-1]-data[pos+1]+data[pos+N-1]*2+data[pos+N])/8;
			gc = -(-data[pos-N-1]*2-data[pos-N]-data[pos-1]+data[pos+1]+data[pos+N]+data[pos+N+1]*2)/8;
			g_2 = sqrt( gr*gr+gc*gc );
			rad_2 =  atan2( gr,gc )+PI/4;
			if( rad_2>PI )
				rad_2-=2*PI;

			g_I[pos] = (g_1+g_2)/2.0;
			g_rad[pos] = (rad_1+rad_2)/2.0;
			if( 1 )	{		//!!testing!!
				I_off+=fabs( g_1-g_2);
				s=fabs( rad_1-rad_2);
				s=MIN( s,2*PI-s );
				rad_off+=s;
			}
		}
	}
	}

	I_off/=MN;		rad_off /= MN;

	search=support;
	e_nz=0;
	M_cut=M/10;
	GIP_MEMCLEAR( w,sizeof(double)*MN ); 
	g_t0 = g_t1/10;
	for( r = M_cut; r < M-M_cut; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		if( pos==7937 )
			i=0;
		if( mask[pos]<0 )			continue;
		if( g_I[pos] < g_t0 )	
			continue;
//		s = fabs(g_rad[pos])*57.295779513082;
		s = g_rad[pos]*57.295779513082;
		if( s>=0 )
			dgr = (22.5<s&&s<=67.5) ? 1 : (67.5<s&&s<=112.5) ? 2 : (112.5<s&&s<=157.5) ? 3 : 0;
		else
			dgr = (-22.5>s&&s>=-67.5) ? 3 : (-67.5>s&&s>=-112.5) ? 2 : (-112.5>s&&s>=-157.5) ? 1 : 0;
		if( g_I[pos]<=g_I[pos+neibor[dgr][0]] || g_I[pos]<=g_I[pos+neibor[dgr][1]])		//nonmaximum suppression
			continue;
/*		if( dgr==0 && (g_I[pos]<=g_I[pos-1] || g_I[pos]<=g_I[pos+1]) )		//nonmaximum suppression
			continue;
		if( dgr==3 && (g_I[pos]<=g_I[pos-N-1] || g_I[pos]<=g_I[pos+N+1]) )
			continue;
		if( dgr==2 && (g_I[pos]<=g_I[pos-N] || g_I[pos]<=g_I[pos+N]) )
			continue;
		if( dgr==1 && (g_I[pos]<=g_I[pos-N+1] || g_I[pos]<=g_I[pos+N-1]) )
			continue;*/
		w[pos]=255.0;
		if( g_I[pos] >= g_t1 )	{
			mask[pos] = 1;			e_nz++;				
		}else	{
			mask[pos] = 10;
			arr_10[nz_2++] = pos;
		}
	}
	}
//	if( 0 )
//		GIP_save_IplImage_d( M,N,".\\trace\\lsu_edge_0.jpg",w,N,0x0 );

	if( 0 )	{
		for( r = M_cut; r < M-M_cut; r++ )	{		
		for( c = 1; c < N-1; c++ )	{
			pos=GIP_RC2POS( r,c,N );
			if( mask[pos]!=10 )			continue;
			d[0] = d[1] = 0.0;		
			g_sum = g_I[pos];
			dr = sin( g_rad[pos] );				dc=cos( g_rad[pos] );	
	//		dr = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_R,0x0 );
	//		dc = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_C,0x0 );
			for( i = 0; i < 2; i++ )	{	
				for( k = 1; k < search; k++ )	{
					GIP_pixel_shift_( r,c,dr,dc,k,&r1,&c1 );
					pos_1=GIP_RC2POS( r1,c1,N );
					if( !GIP_RC_VALID(r1,0,M-1) || !GIP_RC_VALID(c1,0,N-1) || g_I[pos_1]>=g_I[pos] /*|| mask[pos_1]>0*/ )	{
						break;
					}else	{
						d[i] = k;		//k==1 ? 0.0 : sqrt( (r-r1)*(r-r1)+(c-c1)*(c-c1) );		
						g_sum += g_I[pos_1];
					}
				}
				dr*=-1;			dc*=-1;
			}
		//	s = MAX( d[0],d[1] );
		//	s_g = (g_I[pos]-g_sum/(1+d[0]+d[1]) );
		//	w[pos] = s_g*s;		//log(s_g+1.0)*s;
			if( d[0]+d[1]+1>support )
			{	mask[pos]=1;			}
		}
		}
	}

	do{
		nz = 0;			
/*		for( r = 1; r < M-1; r++ )	{		
		for( c = 1; c < N-1; c++ )	{
			pos=GIP_RC2POS( r,c,N );
			if( mask[pos]!=10 )		continue;
			mask[pos]=1;
			if( GIP_adjacent_2( M, N,N,mask,r,c,adj,GIP_ADJACENT_M )>0 )	{
				mask[pos]=1;		nz++;		
			}else
				mask[pos]=10;		
		}
		}
		e_nz += nz;*/
		for( i = nz_1; i < nz_2; i++ )	{
			pos=arr_10[i];	
			if( mask[pos]!=10 )		continue;
			r=GIP_POS2R( pos,N );			c=GIP_POS2C( pos,N );
			mask[pos]=1;
			if( GIP_adjacent_2( M, N,N,mask,r,c,adj,GIP_ADJACENT_M )>0 )	{
				mask[pos]=1;		nz++;		
			}else
				mask[pos]=10;				
		}
	}while( nz>0 );

	if( 0 )	{
		for( r = 1; r < M-1; r++ )	{		
		for( c = 1; c < N-1; c++ )	{
			pos=GIP_RC2POS( r,c,N );
			if( mask[pos]==1 )	
				w[pos]=255.0;
			else
				w[pos]=0.0;
		}
		}
//		D_ouput( M*N,data,".\\trace\\1.dat" );
//		GIP_save_IplImage_d( M,N,".\\trace\\lsu_edge.jpg",w,N,0x0 );
	}

	GIP_free( arr_10 );
	GIP_free( w );
	return e_nz;
}

/*
	注意：
		downsample可利用原有的data存储

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/12/2008	
*/
void GIP_MATRIX_downsample( GIP_MATRIX *hMat,int s_0,int s_1,int step,int flag )	{
	int isRowOri=1,r,c,M=hMat->M,N=hMat->N,nz,pos_0;
	GIP_FLOATING *data=hMat->data;

	ASSERT( s_0>=0 && step>=1 );
	nz=0;
	if( BIT_TEST( flag,GIP_MATRIX_COL) )	{
		ASSERT( s_1<N );
		for( r=0; r<M; r++ )	{
		for( c=s_0; c<=s_1; c+=step )	{
			pos_0=GIP_RC2POS( r,c,N );
			ASSERT( pos_0 >= nz );
			data[nz] = data[pos_0];
			nz++;
		}
		}
		hMat->N = (int)floor((s_1-s_0)/step+1);
	}else if( BIT_TEST( flag,GIP_MATRIX_ROW) ){
		ASSERT( s_1<M );
		for( r=s_0; r<=s_1; r+=step )	{
		for( c=0; c<N; c++ )	{
			pos_0=GIP_RC2POS( r,c,N );
			ASSERT( pos_0 >= nz );
			data[nz] = data[pos_0];
			nz++;
		}
		}
		hMat->M = (int)floor((s_1-s_0)/step+1);
	}

	ASSERT( nz==hMat->M*hMat->N );
}

/*
	注意:	
		1 要注意(N,M)<pad的情况

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/13/2008	
*/
void GIP_MATRIX_symextend( GIP_MATRIX *hMat, GIP_MATRIX *hSrc,int pad,int *k_M,int *k_N,int flag )	{
	int i,j,pos,r,c,M=hSrc->M,N=hSrc->N;
	GIP_FLOATING *u_pad,*data=hSrc->data;
//	keep = floor((fl + size(x) - 1) / 2);
	*k_M = (int)floor( (pad+M)/2 );			*k_N = (int)floor( (pad+N)/2 ); 
//	y = padarray(x, [(fl - 1) (fl - 1)], 'symmetric', 'both');
	u_pad=GIP_alloc( sizeof(GIP_FLOATING)*(M+2*pad)*(N+2*pad) );		

	for( i=0; i < M+2*pad; i++ )	{
	//	r = i<pad ? pad-i : i<M+pad ? i-pad : M-(i-(M+pad) );		
		r = i<pad ? (pad-i-1)%M : i<M+pad ? i-pad : M-1-(i-(M+pad))%M;		
		for( j=0; j < N+2*pad; j++ )	{
			pos=i*(N+2*pad)+j;
	//		c = j<pad ? pad-j : j<N+pad ? j-pad : N-(j-(N+pad) );		
			c = j<pad ? (pad-j-1)%N : j<N+pad ? j-pad : N-1-(j-(N+pad))%N;		
			ASSERT( r<M && c<N && r>=0 && c>=0 );
			u_pad[pos] = data[r*N+c];
		}
	}

	hMat->data=u_pad;
	hMat->M=M+2*pad;
	hMat->N=N+2*pad;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/13/2008	
*/
void GIP_MATRIX_fill( GIP_MATRIX *hMat,GIP_MATRIX *hFill,int m_0,int n_0,int flag )	{
	int r,c,pos,m_1,n_1;
	GIP_FLOATING *data=hMat->data,*dataF=hFill->data;

	ASSERT( m_0>=0 && n_0>=0 );
	ASSERT( m_0+hFill->M <= hMat->M && n_0+hFill->N <= hMat->N );

	m_1=m_0+hFill->M;			n_1=n_0+hFill->N;			
	for( r=m_0; r<m_1; r++ )	{
	for( c=n_0; c<n_1; c++ )	{
		pos=GIP_RC2POS( r,c,hMat->N );
		data[pos] = *dataF;
		dataF++;
	}
	}

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/5/2008	
*/
void GIP_OBJ_trans( GIP_OBJECT *hObj,int type,int flag )	{
	int i,M=hObj->M,N=hObj->N,r,c,ldu=N,pos;
	GIP_FLOATING *u=hObj->data,u1,u2,u3,u4,a,a_max,a_min,*temp;

	temp=GIP_alloc( sizeof(double)*M*N );
	switch( type )	{
	case GIP_TRANS_LAPLACIAN:
		for( r=0; r<M; r++ )	{
		for( c=0; c<N; c++ )	{
			pos=GIP_RC2POS( r,c,ldu );
			u1=c==0 ? u[pos+1] : u[pos-1];
			u2=c==N-1 ? u[pos-1] : u[pos+1];
			u3=r==0 ? u[pos+ldu] : u[pos-ldu];
			u4=r==M-1 ? u[pos-ldu] : u[pos+ldu];
			temp[pos] = u1+u2+u3+u4-u[pos]*4;
		}
		}
		break;
	case GIP_TRANS_GRADIENT:
		for( r=0; r<M; r++ )	{
		for( c=0; c<N; c++ )	{
			pos=GIP_RC2POS( r,c,ldu );
			u1=c==0 ? u[pos+1] : u[pos-1];
			u2=c==N-1 ? u[pos-1] : u[pos+1];
			u3=r==0 ? u[pos+ldu] : u[pos-ldu];
			u4=r==M-1 ? u[pos-ldu] : u[pos+ldu];
			temp[pos] = fabs(u1-u2)/2+fabs(u3-u4)/2;
		}
		}
		break;
	default:
		ASSERT(0);
		break;
	}

	if( BIT_TEST(flag,GIP_TRANS_ABS ) )	{
		for( i = 0; i < M*N; i++ )	{
			a = fabs(temp[i]);
			u[i]=a;
		}
	}else	{
		for( i = 0; i < M*N; i++ )	u[i]=temp[i];
	}
	if( 1 )	{
		a_max=DBL_MIN;			a_min=DBL_MAX;
		for( i = 0; i < M*N; i++ )	{
			if( a_min > u[i] )	{
				a_min=u[i];
			}
			if( a_max < u[i] )	{
				a_max=u[i];
			}		
		}
	}

	GIP_free( temp );
	if( 1 )	{
		char sTitle[80];
		_stprintf( sTitle,".\\trace\\trans_%d_1.jpg",M*N );	
//		GIP_save_IplImage_d( hObj->M,hObj->N,sTitle,hObj->data,hObj->N,0x0 );
	}
}

/*
	Histogram Equalization( Contrast Enhancement )
	temp[MN*2]

	v0.1	cys
		1/7/2009		
		2/19/2009		消除若干bug		
*/
int GIP_Histo_Equal(  GIP_OBJECT *hObj,int L,int type,int *temp,int flag )	{
	int i,it,M=hObj->M,N=hObj->N,MN=M*N,ret=GIP_OK;
	int *n,*s,isAlloc=temp==GIP_NULL,I_0,I_1,nz_0=-1,no_0;
	GIP_FLOATING *data=hObj->data,a;

	if( isAlloc )	temp = GIP_alloc( sizeof(int)*MN*2 );
	n=temp;		s=temp+MN;
	for( i = 0; i < MN; i++ )	n[i]=0;
//	I_0=INT_MAX,			I_1=0;
	for( i = 0; i < MN; i++ )	{
		it = GIP_INTENSITY(data[i]);		
//		I_0=MIN(I_0,it);			I_1=MAX(I_1,it);
		n[it]++;
	}
	for( i = 0; i < L; i++ )	{
		if( n[i]!=0 )	{
			nz_0=n[i];		no_0=i;
			break;
		}
	}
	if( nz_0==MN )	{
		ret = -(no_0+1);
		goto GIP_exit;
//		GIP_ERROR( _T("GIP_Histo_Equal"),-1, sInfo );
	}
	s[0]=n[0];
	for( i = 1; i < L; i++ )	{
		s[i]=s[i-1]+n[i];
	}
	for( i = no_0; i < L; i++ )	{
		a=(L-1)*1.0*(s[i]-nz_0)/(MN-nz_0);					
		s[i]=GIP_INTENSITY(a);
		ASSERT( s[i]>=0 && s[i] < L );
	}
	I_0=INT_MAX,			I_1=0;
	for( i = 0; i < MN; i++ )	{
		it = GIP_INTENSITY(data[i]);
		data[i] = s[it];
		I_0=MIN(I_0,data[i]);			I_1=MAX(I_1,data[i]);
	}	
	ASSERT( I_0==0 && I_1==L-1 );
//	GIP_save_IplImage_d( hObj->M,hObj->N,".\\trace\\histo.jpg",hObj->data,hObj->N,0x0 );
GIP_exit:	;
	if( isAlloc )		GIP_free( temp );

	return ret;
}

/*
	Adaptive Contrast Enhancement
	[REF]: A Fast and Adaptive Method for Image Contrast Enhancement 

	ISSUE-NEEDED
		1 w0的自动确定及W的范围 [0..L]

	v0.1	cys
		1/7/2008		
*/
void GIP_Contrast_Enhance(  GIP_OBJECT *hObj,int L,int type,int flag )	{
	int i,pos,r,c,M=hObj->M,N=hObj->N,MN=M*N,ldu=N;
	GIP_FLOATING *data=hObj->data,a,*D_temp,*I_a,*I_0,*I_1;
	double C=0.95,w_0=50.0,W;

	D_temp = GIP_alloc( sizeof(GIP_FLOATING)*MN*3 );
	I_a=D_temp;		I_0=D_temp+MN;		I_1=D_temp+MN*2;	
	GIP_MEMCOPY( I_a,data,sizeof(GIP_FLOATING)*MN );
	GIP_MEMCOPY( I_0,data,sizeof(GIP_FLOATING)*MN );
	GIP_MEMCOPY( I_1,data,sizeof(GIP_FLOATING)*MN );

	for( r=0; r<M; r++ )	{		//propagation in single row
		for( c=1; c<N; c++ )	{
			pos=GIP_RC2POS( r,c,ldu );
			I_a[pos]=(1-C)*I_a[pos]+C*I_a[pos-1];	
			if( I_0[pos-1] < I_0[pos] )		I_0[pos]=(1-C)*I_0[pos]+C*I_0[pos-1];	
			if( I_1[pos-1] > I_1[pos] )		I_1[pos]=(1-C)*I_1[pos]+C*I_1[pos-1];	
		}
		for( c=N-2; c>=0; c-- )	{
			pos=GIP_RC2POS( r,c,ldu );
			I_a[pos]=(1-C)*I_a[pos]+C*I_a[pos+1];	
			if( I_0[pos+1] < I_0[pos] )		I_0[pos]=(1-C)*I_0[pos]+C*I_0[pos+1];	
			if( I_1[pos+1] > I_1[pos] )		I_1[pos]=(1-C)*I_1[pos]+C*I_1[pos+1];	
		}	
	}
	for( c=0; c<N; c++ )	{		//propagation in single column
		for( r=1; r<M; r++ )	{
			pos=GIP_RC2POS( r,c,ldu );
			I_a[pos]=(1-C)*I_a[pos]+C*I_a[pos-ldu];	
			if( I_0[pos-ldu] < I_0[pos] )		I_0[pos]=(1-C)*I_0[pos]+C*I_0[pos-ldu];	
			if( I_1[pos-ldu] > I_1[pos] )		I_1[pos]=(1-C)*I_1[pos]+C*I_1[pos-ldu];	
		}
		for( r=M-2; r>=0; r-- )	{
			pos=GIP_RC2POS( r,c,ldu );
			I_a[pos]=(1-C)*I_a[pos]+C*I_a[pos+ldu];	
			if( I_0[pos+ldu] < I_0[pos] )		I_0[pos]=(1-C)*I_0[pos]+C*I_0[pos+ldu];	
			if( I_1[pos+ldu] > I_1[pos] )		I_1[pos]=(1-C)*I_1[pos]+C*I_1[pos+ldu];	
		}	
	}

	for( i = 0; i < MN; i++ )	{
		a = data[i];
		ASSERT( I_a[i]>= I_0[i] && I_a[i]<= I_1[i] );
		ASSERT( a>= I_0[i] && a<= I_1[i] );
		if( a <= w_0 )	{
			W = w_0-sqrt(w_0*w_0-a*a);
		}else	{
			a=(255-w_0)*(255-w_0)+(255-a)*(255-a);
			W = w_0+sqrt(a);
		}
		a = W*(data[i]-I_0[i])/(I_1[i]-I_0[i]);
		data[i] = MIN( GIP_INTENSITY(a),255 );
	}	

//	GIP_save_IplImage_d( hObj->M,hObj->N,".\\trace\\enhance.jpg",hObj->data,hObj->N,0x0 );

	GIP_free( D_temp );
}

/*
	Edge Detection

	ISSUE-NEEDED:
		线条过粗

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/12/2008	
*/
void GIP_EdgeMap( GIP_OBJECT *hObj,GIP_OBJECT *hEdge,double thrsh,int type,int flag )	{
	int i,M=hObj->M,N=hObj->N,MN=M*N,ldu=N;
	double gx,gy,g;
	GIP_FLOATING *data=hObj->data,*edge=hEdge->data;

	for( i = 0; i < MN; i++ )	{
		gx = GIP_Mask_3_( i,M,N,N,data,G_MASK_SOBEL_R,0x0 );
		gy = GIP_Mask_3_( i,M,N,N,data,G_MASK_SOBEL_C,0x0 );
		g = sqrt( gx*gx+gy*gy );
		if( g > thrsh )	{
			edge[i]=0.0;
		}else
			edge[i]=1.0;
	}

//	GIP_save_IplImage_d( hEdge->M,hEdge->N,".\\trace\\edge.jpg",hEdge->data,hEdge->N,0x0 );

}

/*
	从edgemap中提取每条边的信息
	注意：
		利用了edgemap中线段用0表示的特性。

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/13/2008	
*/
int GIP_Edge2List( GIP_OBJECT *hEdge,int *E_ptr,int *E_pos,int *I_buffer,int flag )	{
	int i,pos,r,c,M=hEdge->M,N=hEdge->N,MN=M*N,ldu=N;
	int nEdge=0,*data,nz=0,nBrother,nAdj,adj[8],l_max;

	data=GIP_alloc( sizeof(int)*MN );
	for( i = 0; i < MN; i++ )	{
		data[i]=(int)(hEdge->data[i]);
	}

	l_max = 0;
	E_ptr[0]=0;
	for( i = 0; i < MN; i++ )	{
		if( data[i]!=0 )	continue;

		nEdge++;
		pos=i;
		nBrother=1;
		while( 1 )	{
			r=GIP_POS2R( pos,ldu );		c=GIP_POS2C( pos,ldu );
			ASSERT( GIP_RC_VALID( r,0,M-1 ) && GIP_RC_VALID( c,0,N-1 ));
			if( nBrother==1 )
				nAdj=GIP_adjacent_2( M,N,N,data,r,c,adj,GIP_ADJACENT_M );	
			else
				nAdj=0;
			E_pos[nz]=pos;			nz++;
			data[pos]=nEdge;			
			if( nAdj==0 )	break;
			pos=adj[0];			
			nBrother=nz==1?1:nAdj;
		}
		E_ptr[nEdge]=nz;
		if( E_ptr[nEdge]>E_ptr[nEdge-1] )
			l_max = E_ptr[nEdge]-E_ptr[nEdge-1];
	}

	GIP_free( data );

	return nEdge;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/31/2008	
*/
int GIP_circle_pos( GIP_OBJECT *hObj,double c_r,double c_c,double c_R,int *arrPos,int flag )	{
	int i,pos,r,c,r1,c1,M=hObj->M,N=hObj->N,nz=0;
	double off_1,off_2;	

	c=(int)(c_R+0.5);
	r=0;		
	while(1)	{		//testing
		int off_r[8]={ r,r,-r,-r,c,c,-c,-c };
		int off_c[8]={ c,-c,c,-c,r,-r,r,-r };
		if( r>c )
			break;
		//mirror to 8 point
		for( i = 0; i < 8; i++ )	{
			r1=c_r+off_r[i];	c1=c_c+off_c[i];
			if( r1<0||r1>=M || c1<0||c1>=N )	
				continue;
			pos = GIP_RC2POS( r1,c1,N );
			arrPos[nz++] = pos;
		}
		r++;
		off_1 = fabs(r*r+c*c-c_R*c_R);
		off_2 = fabs(r*r+(c-1)*(c-1)-c_R*c_R);
		if( off_1>=off_2 )
			c--;
	}

	return nz;
}

/*
	从edgemap中提取每条边的信息
	注意：
		1 需替换为midpoint circle algorithm。
		2 允许圆出界的情况

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/13/2008	
*/
void GIP_Circle( GIP_OBJECT *hObj,int c_r,int c_c,double c_R,GIP_IMAG_TYPE type )	{
	int i,pos,r,c,r1,c1,M=hObj->M,N=hObj->N,MN=M*N,ldu=N;
	GIP_FLOATING *data=hObj->data,back,fore;
	double off_1,off_2;

//	ASSERT( c_r+c_R<M && c_r-c_R>=0 && c_c+c_R<N && c_c-c_R>=0 );
	switch( type )	{
	case GIP_IMAG_EDGE:
		back=0.0,		fore=1.0;
		break;
	case GIP_IMAG_ORIGIN:
		back=-1.0,		fore=0.0;
		break;
	default:
		ASSERT( 0 );
		break;
	}
	if( back != -1.0 )
		for( i = 0; i < MN; i++ )		data[i]=back;
	if( 1 )	{
		int *arrPos=GIP_alloc( sizeof(int)*MN ),nz;
		nz = GIP_circle_pos( hObj,c_r,c_c,c_R,arrPos,0x0 );
		for( i = 0; i < nz; i++ )
			data[arrPos[i]]=fore;
		GIP_free( arrPos );
	}else	{
		c=(int)(c_R+0.5);
		r=0;		
		while(1)	{		//testing
			int off_r[8]={ r,r,-r,-r,c,c,-c,-c };
			int off_c[8]={ c,-c,c,-c,r,-r,r,-r };
			if( r>c )
				goto FINISH;
			//mirror to 8 point
			for( i = 0; i < 8; i++ )	{
				r1=c_r+off_r[i];	c1=c_c+off_c[i];
				if( r1<0||r1>=M || c1<0||c1>=N )	continue;
				pos = GIP_RC2POS( r1,c1,N );
				data[pos]=fore;
			}
			r++;
			off_1 = fabs(r*r+c*c-c_R*c_R);
			off_2 = fabs(r*r+(c-1)*(c-1)-c_R*c_R);
			if( off_1>=off_2 )
				c--;
		}
		}
FINISH:
	;
//	GIP_save_IplImage_d( hObj->M,hObj->N,".\\trace\\circle.jpg",hObj->data,hObj->N,0x0 );
}

/*
	双曲线：	(c-h)*(c-h)=a(r-k)

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/17/2009	
*/
void GIP_Parabolic( GIP_OBJECT *hObj,double a,double h,double k,GIP_IMAG_TYPE type )	{
	int i,pos,r,c,M=hObj->M,N=hObj->N,MN=M*N,ldu=N;
	GIP_FLOATING *data=hObj->data,back,fore;
	double w;

	if( a==0.0 )
		return;
//	ASSERT( c_r+c_R<M && c_r-c_R>=0 && c_c+c_R<N && c_c-c_R>=0 );
	switch( type )	{
	case GIP_IMAG_EDGE:
		back=0.0,		fore=1.0;
		break;
	case GIP_IMAG_ORIGIN:
		back=-1.0,		fore=0.0;
		break;
	default:
		ASSERT( 0 );
		break;
	}
	if( back != -1.0 )
		for( i = 0; i < MN; i++ )		data[i]=back;

	for( c=0; c<N; c++ )	{
		w=(c-h)*(c-h)*1.0/a+k;
		r = G_DOUBLE2INT( w );
		if( r<0 || r>=M )		continue;
//		e_nz++;
		pos=GIP_RC2POS( r,c,N );
		data[pos]=fore;
/*		while( r < M )	{
			pos=GIP_RC2POS( r,c,N );
			data[pos]=fore;
			r++;
		}*/
	}
	;
}

/*
	基于梯度方向的校验

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/20/2008	
	v0.2	cys		返回measure
		3/13/2009	
	v0.3	cys		有时g_I,g_rad,mask已知
		3/29/2009	
*/
double GIP_Circle_GMeasure( GIP_OBJECT *hObj,int c_r,int c_c,double c_R,double xita,
						   GIP_FLOATING *g_I,GIP_FLOATING*g_rad,int *mask,int flag )	{
	int i,pos,r,c,M=hObj->M,N=hObj->N,MN=M*N,ldu=N,ret=0x0;
	int r1,c1,*E_ptr,*E_pos,nEdge,nz,nSec,isMatch=1,nz_1;
	GIP_FLOATING *data=hObj->data,*e_data,back=1.0,fore=0.0;
	GIP_OBJECT *hEdge=GIP_NULL;
	double off_1,off_2,*g,gr,gc,a_1,a_2,a,s_match=-1.0;

	ASSERT( g_I!=GIP_NULL && g_rad!=GIP_NULL );
	xita*=1.5;	//考虑到取整，拟合等因素，thrsh放大1.5倍

	if( BIT_TEST( flag,GIa_INNER_OBJ) &&
		!(c_r+c_R<M && c_r-c_R>=0 && c_c+c_R<N && c_c-c_R>=0) )	{
		s_match=-1.0;		goto END;
	}
	hEdge = GIP_OBJ_create( M,N,GIP_NULL,0x0 );
	e_data=hEdge->data;
	for( i = 0; i < MN; i++ )		e_data[i]=back;

	c=(int)(c_R+0.5);				r=0;
	nz = 0;
	while(1)	{		//testing
		int off_r[8]={ r,r,-r,-r,c,c,-c,-c };
		int off_c[8]={ c,-c,c,-c,r,-r,r,-r };
		if( r>c )
			break;
		//mirror to 8 point
		for( i = 0; i < 8; i++ )	{
			r1=c_r+off_r[i];	c1=c_c+off_c[i];
			if( r1<0||r1>=M || c1<0||c1>=N )	continue;
			pos = GIP_RC2POS( r1,c1,N );
			if( e_data[pos]!=fore )
				nz++;
			e_data[pos]=fore;
		}
		r++;
		off_1 = fabs(r*r+c*c-c_R*c_R);
		off_2 = fabs(r*r+(c-1)*(c-1)-c_R*c_R);
		if( off_1>=off_2 )
			c--;
	}
	E_ptr=GIP_alloc(sizeof(int)*MN*2);		E_pos=E_ptr+MN;	
	g=GIP_alloc( sizeof(double)*MN );
	nEdge = GIP_Edge2List( hEdge,E_ptr,E_pos,GIP_NULL,0x0 );
	ASSERT( E_ptr[nEdge]==nz );
	nz = 0;			nSec=0;				nz_1=0;
	for( i = 0; i < E_ptr[nEdge]; i++ )	{
		pos = E_pos[i];
		r=GIP_POS2R( pos,N ),		c=GIP_POS2C( pos,N );
/*		a_1 = atan2( r-c_r,c-c_c );
		gr = GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_R,0x0 );
		gc = GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_C,0x0 );
		g[i] = sqrt( gr*gr+gc*gc );
		a_2 = atan2( gr,gc );*/
		if( mask!=GIP_NULL && mask[pos]<0 )	{
			continue;
		}
		g[i]=g_I[pos];
		a_1 = atan2( c_r-r,c_c-c );
		a_2=g_rad[pos];
		a = fabs(a_1-a_2);				
		if( a>2*PI )		
			a -= 2*PI;
		a = MIN(a,2*PI-a);
		if( a<xita )	{	
			if( mask!=GIP_NULL && mask[pos]==1 )
				nz_1++;
			isMatch=1;		nz++;
		}	else	{
			if( isMatch==1 )
				nSec++;
			isMatch = 0;
		}
	}
//	s_match = nz_1<=nz/10 ? nz_1*1.0/E_ptr[nEdge] : (nz+nz_1)*0.5/E_ptr[nEdge];
	s_match = nz*1.0/E_ptr[nEdge];
/*	if( nz*1.0/E_ptr[nEdge]<thrsh )	{
		ret=-1;		
	}
	for( i = 0; i < nEdge; i++ )	{
		g_mean = 0.0;
		nz = E_ptr[i+1]-E_ptr[i];
		for( j = E_ptr[i]; j < E_ptr[i+1]; j++ )	{
			g_mean += g[j];
		}
		g_mean /= nz;
		g_devia = 0.0;
		for( j = E_ptr[i]; j < E_ptr[i+1]; j++ )	{
			g_devia += (g[j]-g_mean)*(g[j]-g_mean);
		}	
		g_devia = sqrt( g_devia/nz );
	}*/

	GIP_free( E_ptr );			GIP_free( g );
	GIP_OBJ_clear( hEdge );		GIP_free( hEdge );
END:	
	return s_match;
//	GIP_save_IplImage_d( hObj->M,hObj->N,".\\trace\\circle.jpg",hObj->data,hObj->N,0x0 );
}

/*
	only for GIP_Smooth_RecurAniso
	equation (4) from [REF]:Recursive Implementation of Anisotropic Filtering

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/9/2009	
*/
void _formular_recur_3( int N,GIP_FLOATING* v_1,GIP_FLOATING *v_2,int step,GIP_FLOATING s1,GIP_FLOATING s2,GIP_FLOATING s3,int type )	{
	int i,pos;
	GIP_FLOATING a,b_1,b_2,b_3;

	switch( type )	{
	case 0:			//local average
		for( i=3; i<N; i++ )	{
			a = v_1[i*step];
			pos = i*step;
			v_2[pos]=a+s1*(v_2[pos-step]-a)+s2*(v_2[pos-2*step]-a)+s3*(v_2[pos-3*step]-a);
		}
		break;
	case 1:			//To compute the local minimum
		for( i=3; i<N; i++ )	{
			a = v_1[i*step];
			pos = i*step;
			b_1 = a>v_2[pos-step] ? v_2[pos-step]-a :0.0;
			b_2 = a>v_2[pos-2*step] ? v_2[pos-2*step]-a :0.0;
			b_3 = a>v_2[pos-3*step] ? v_2[pos-3*step]-a :0.0;
			v_2[pos]=a+s1*b_1+s2*b_2+s3*b_3;
		}
		break;
	case 2:			//To compute the local maximum
		for( i=3; i<N; i++ )	{
			a = v_1[i*step];
			pos = i*step;
			b_1 = a>v_2[pos-step] ? a-v_2[pos-step] : 0.0;
			b_2 = a>v_2[pos-2*step] ? a-v_2[pos-2*step] : 0.0;
			b_3 = a>v_2[pos-3*step] ? a-v_2[pos-3*step] : 0.0;
			v_2[pos]=a+s1*b_1+s2*b_2+s3*b_3;
			ASSERT( v_2[pos] > 0.0 );
		}
		break;
	default:
		ASSERT( 0 );
		break;
	}
}

/*
	Gaussian filter需针对不同的情况选择不同的q
	对于边界比较清楚,纹理较弱,q取小值(如MMU iris database)
	对于边界比较模糊,纹理较强,q取大值(如CASIA)

	出现了一个非常奇怪的bug		2/10/2009

	[REF]
	1	Recursive Implementation of Anisotropic Filtering (Final Report) 
	2	Recursive implementation of the Gaussian filter 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/9/2009	
*/
void GIP_Smooth_RecurAniso( GIP_OBJECT *hObj,GIP_FLOATING *out,double q,int flag )	{
	int i,r,c,M=hObj->M,N=hObj->N,MN=M*N,ldu=N;
	GIP_FLOATING *data=hObj->data,a,*D_temp,*cur,*l_avg,*l_max,*l_min,*w;
	double s1,s2,s3;

	if( 1 )	{
		double b0,b1,b2,b3;
//		double b0,b1,b2,b3,q=1.0;
		b0=1.57825+2.44413*q+1.4281*q*q+0.422205*q*q*q;
		b1=		  2.44413*q+2.85619*q*q+1.26661*q*q*q;
		b2=					-1.4281*q*q-1.26661*q*q*q;
		b3=								0.422205*q*q*q;

		s1=b1/b0;		s2=b2/b0;		s3=b3/b0;	
	}else	{
		s1=2.36565;		s2=-1.89709;		s3=0.51601;		//for q=5.0
	}

	D_temp = GIP_alloc( sizeof(GIP_FLOATING)*MN*4 );
	l_avg=D_temp;		l_max=D_temp+MN;		l_min=D_temp+MN*2;	
	w=D_temp+MN*3;	
	GIP_MEMCOPY( l_avg,data,sizeof(GIP_FLOATING)*MN );
	GIP_MEMCOPY( l_min,data,sizeof(GIP_FLOATING)*MN );
	GIP_MEMCOPY( l_max,data,sizeof(GIP_FLOATING)*MN );

	for( c=0; c<N; c++ )	{		//propagation in single column
		cur = l_avg+c;
		w[0]=cur[0];	w[N]=cur[N];	w[2*N]=cur[2*N];	
		_formular_recur_3( M,cur,w,N,s1,s2,s3,0 );		
		_formular_recur_3( M,w+(M-1)*N,cur+(M-1)*N,-N,s1,s2,s3,0 );	

/*		cur = l_min+c;
		w[0]=cur[0];	w[N]=cur[N];	w[2*N]=cur[2*N];	
		_formular_recur_3( M,cur,w,N,s1,s2,s3,1 );		
		_formular_recur_3( M,w+(M-1)*N,cur+(M-1)*N,-N,s1,s2,s3,1 );		

		cur = l_max+c;
		w[0]=cur[0];	w[N]=cur[N];	w[2*N]=cur[2*N];	
		_formular_recur_3( M,cur,w,N,s1,s2,s3,2 );		
		_formular_recur_3( M,w+(M-1)*N,cur+(M-1)*N,-N,s1,s2,s3,2 );	*/		//非常奇怪的bug,将2改为1即可避免内存泄露	
	}	
	for( r=0; r<M; r++ )	{		//propagation in single row
		cur = l_avg+r*N;
		w[0]=cur[0];	w[1]=cur[1];	w[2]=cur[2];	
		_formular_recur_3( N,cur,w,1,s1,s2,s3,0 );		
		_formular_recur_3( N,w+N-1,cur+N-1,-1,s1,s2,s3,0 );		

	/*	cur = l_min+r*N;
		w[0]=cur[0];	w[1]=cur[1];	w[2]=cur[2];	
		_formular_recur_3( N,cur,w,1,s1,s2,s3,1 );		
		_formular_recur_3( N,w+N-1,cur+N-1,-1,s1,s2,s3,1 );	

		cur = l_max+r*N;
		w[0]=cur[0];	w[1]=cur[1];	w[2]=cur[2];	
		_formular_recur_3( N,cur,w,1,s1,s2,s3,2 );		
		_formular_recur_3( N,w+N-1,cur+N-1,-1,s1,s2,s3,2 );	*/	
	}

//	GIP_save_IplImage_d( M,N,".\\trace\\RAG_filter.jpg",l_avg,N,0x0 );
	for( i = 0; i < MN; i++ )	{
		a = data[i];
//		ASSERT( l_max[i]>= l_avg[i] && l_avg[i]>=l_min[i] );
//		ASSERT( a>= I_0[i] && a<= I_1[i] );
		if( a >= l_avg[i] )	{
			a = l_max[i];
		}else	{
			a = l_min[i];
		}
		a = l_avg[i];
		out[i] = MIN( GIP_INTENSITY(a),255 );	
	}	
//	GIP_save_IplImage_d( hObj->M,hObj->N,".\\trace\\RAG_filter.jpg",out,hObj->N,0x0 );

	GIP_free( D_temp );
}

/*
	q为控制参数

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/13/2008	
*/
void GIP_Smooth( GIP_OBJECT *hObj,int type,GIP_FLOATING *d_temp,double q,int flag )	{
	int i,j,pos,r,c,M=hObj->M,N=hObj->N,MN=M*N,ldu=N;
	int off[9]={-ldu-1,-ldu,-ldu+1,-1,0,1,ldu-1,ldu,ldu+1};
	double w[9]={1,2,1,2,4,2,1,2,1},w_sum=16;
//	double LAP_[9]={-1,-1,-1,-1,11,-1,-1,-1,-1};
	double LAP_[9]={0,-1,0,-1,5,-1,0,-1,0};
	GIP_FLOATING *data=hObj->data,a,*vec;

	vec=d_temp+MN;
	for( i = 0; i < MN; i++ )	d_temp[i]=0.0;
	switch( type )	{
	case GIP_SMOOTH_MEDIAN:
		for( r=1; r<M-1; r++ )	{		
		for( c=1; c<N-1; c++ )	{
			pos=GIP_RC2POS( r,c,ldu );
			for( i = 0; i < 9; i++ )	vec[i] = data[pos+off[i]];
			for( i = 0; i < 9; i++ )	{
				for( j = i+1; j < 9; j++ )	{
					if( vec[i]>vec[j] )
					{	a=vec[i];	vec[i]=vec[j];		vec[j]=a;	}
				}			
			}
			d_temp[pos]=vec[5];
		}
		}
		break;
	case GIP_SMOOTH_GAUSSIAN:
		break;
	case GIP_SMOOTH_LAPLACIAN:
		for( r=1; r<M-1; r++ )	{		
		for( c=1; c<N-1; c++ )	{
			pos=GIP_RC2POS( r,c,ldu );
			a = 0.0;
			for( i = 0; i < 9; i++ )	
				a += data[pos+off[i]]*LAP_[i];
			d_temp[pos]=a;
		}
		}
		break;
	case GIP_SMOOTH_ANISO:
		GIP_Smooth_RecurAniso( hObj,d_temp,q,flag );
		break;
	default:
		for( r=1; r<M-1; r++ )	{		
		for( c=1; c<N-1; c++ )	{
			pos=GIP_RC2POS( r,c,ldu );
			a = 0.0;
			for( i = 0; i < 9; i++ )	{
				a += data[pos+off[i]]*w[i];
			}
			d_temp[pos]=a/w_sum;
		}
		}
		break;
	};

	for( i = 0; i < MN; i++ )	data[i]=d_temp[i];
//	GIP_save_IplImage_d( hObj->M,hObj->N,".\\trace\\smooth.jpg",hObj->data,hObj->N,0x0 );
}

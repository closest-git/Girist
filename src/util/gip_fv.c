#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include "../PCH/GIP_def.h"
#include "gip_util.h"
#include "gip_fv.h"

/*
	v0.1	cys
		3/24/2009
*/
void _dump_16_( )	{
	unsigned int i,j;
	int t;

	G_PRINTF( _T("\r\n") );
	for( i = 0; i < USHRT_MAX; i++ )	{
		t = 0;
		for( j=0; j < 32; j++ )	{
			t += BIT_TEST( i,BIT_FLAG_1<<j ) ? 1.0:0.0;
		}
		G_PRINTF( _T("%d,"),t );
		if( (i+1)%64==0 )
			G_PRINTF( _T("\r\n") );
		//	no[t*64+nz[t]]=i;
	//	nz[t]++;
	}

	G_PRINTF( _T("\r\n") );
}

/*
	ISSUE-NEEDED:
	1 why subfeatures is 781 in [REF]:Robust Iris Feature Extraction and Matching

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/3/2009	
*/
void G_FVector_init( G_FVector_ *hFV,int M,int N,int nShift,int type )	{
	int nz_M,nz_N,off;

	GIP_MEMCLEAR( hFV,sizeof(G_FVector_) );
	switch( type )	{
	default:			//double overlap
		hFV->sf_M=8;			hFV->sf_N=12;
		hFV->M_step = hFV->sf_M/2;
		hFV->N_step = hFV->sf_N/2;
		break;
	}
	
	hFV->nShift = nShift;
	off=hFV->sf_M-hFV->M_step;			nz_M = (M-off)/hFV->M_step;
	off=hFV->sf_N-hFV->N_step;			nz_N = (N-off)/hFV->N_step;
	hFV->bpSF = 8;
	hFV->nzSF = nz_M*nz_N;
	hFV->M_step = 1;					hFV->N_step = 1;

//	hFV->type = G_FV_COEFFICIENT;		
//	hFV->type = G_FV_ROI;
//	hFV->type = G_FV_DCT_CODE;
//	hFV->type = G_FV_MASEK;
//	hFV->type = G_FV_CYS;
//	hFV->type = G_FV_CYS_1;
	hFV->type = G_FV_CYS_HALF;
	switch( hFV->type )	{
	case G_FV_ROI:
		hFV->sf_M=M;			hFV->sf_N=N;
		hFV->unit_size=sizeof(double);
		hFV->V_len = M*N;
		break;
	case G_FV_MALI:
		break;
	case G_FV_MASEK:
		hFV->unit_size=sizeof(char);
		hFV->sf_M=M;			hFV->sf_N=N;
		hFV->V_len = M*N*2/8;
		break;
	case G_FV_CYS:
	case G_FV_CYS_HALF:		//half of cys
		hFV->M_step = 4;					hFV->N_step = 1;
		hFV->unit_size=sizeof(char);
		hFV->sf_M=M/hFV->M_step;			hFV->sf_N=N;
		hFV->V_len = hFV->sf_M*hFV->sf_N/8;
		break;
	case G_FV_CYS_1:
		hFV->M_step = 4;					hFV->N_step = 1;
		hFV->unit_size=sizeof(char);
		hFV->sf_M=M/hFV->M_step;			hFV->sf_N=N;
		hFV->V_len = hFV->sf_M*hFV->sf_N*2/8;
		break;
	case G_FV_COEFFICIENT:
		hFV->bpSF=2;
		hFV->unit_size=sizeof(double);
		hFV->V_len = nz_M*nz_N*hFV->bpSF;
		break;
	default:
		hFV->unit_size=sizeof(char);
		hFV->V_len = hFV->bpSF*hFV->nzSF/8;
	}
	ASSERT( hFV->V==GIP_NULL );
	hFV->V= GIP_alloc( hFV->unit_size*hFV->V_len );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/10/2009	
*/
void G_FVector_set_step( G_FVector_ *hFV,int M_step,int N_step )	{
	int M=hFV->sf_M*hFV->M_step,N=hFV->sf_N*hFV->N_step;
	if( M_step==hFV->M_step && N_step==hFV->N_step )
		return;

	hFV->sf_M=M/M_step;
	hFV->sf_N=N/N_step;
	hFV->M_step=M_step;				hFV->N_step=N_step;
	if( hFV->type==G_FV_CYS || hFV->type==G_FV_CYS_HALF )
		hFV->V_len = hFV->sf_M*hFV->sf_N/8;
	else
		hFV->V_len = hFV->sf_M*hFV->sf_N*2/8;

	ASSERT( hFV->V!=GIP_NULL );
	GIP_free( hFV->V );
	hFV->V= GIP_alloc( hFV->unit_size*hFV->V_len );

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/3/2009	
*/
void G_FVector_clear( G_FVector_ *hFV )	{
	if( hFV->V != GIP_NULL )		{	
		GIP_free( hFV->V );				hFV->V=GIP_NULL; 
	}
	if( hFV->mask != GIP_NULL )		{	
		GIP_free( hFV->mask );			hFV->mask=GIP_NULL; 
	}
	GIP_MEMCLEAR( hFV,sizeof(G_FVector_) ); 
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/7/2009		
*/
double _euclid_distance_( int v_len,double *V1,double *V2,int flag )	{
	int i;
	double dist=0.0,a;

	for( i = 0; i < v_len; i++ )	{
		a = V1[i]-V2[i];
		dist += a*a;
	}

	dist = sqrt( dist/v_len );

	return dist;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/3/2009	
*/
double _Hamming_distance_1_( int v_len,char *V1,char *V2,int bpSF,int nzSF,char *temp )	{
	double a,sum=1.0,hd=0.0;
	int byte,bit,i,j;

	for( i = 0; i < v_len; i++ )
		temp[i] = V1[i]^V2[i];
	for( i = 0;	i < bpSF; i++ )	{	
		a = 0.0;
		for( j = 0; j < nzSF; j++ )	{
			byte = j;					bit=i;
//			byte = (j*bpSF+i)/8;		bit=(j*bpSF+i)%8;
			a += BIT_TEST( temp[byte],BIT_FLAG_1<<bit ) ? 1.0:0.0;
		}
		sum *= (a/nzSF);
	}
	hd = pow( sum,1.0/bpSF );		

	return hd;
}

/*	
	注意:
		1 sub_len,v_len,pos_0以byte(8 bit)为单位,off以bit为单位

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/26/2009	
*/
void G_fv_shift_( G_FV_TYPE type,int nShift,int v_len,G_U8 *V1,G_U8 *V_shift,int flag )	{
	int bit,r,c,j,cur,base,pos_1,M,pos_0,off,grid=2;
	int sub_len,sft;
	G_U8 temp,*v_cur;

	switch( type )	{
	case G_FV_CYS_HALF:		
	case G_FV_CYS:		
		sub_len=512/8;		sft=1*grid;
		break;
	case G_FV_CYS_1:		
		sub_len=512/4;		sft=2*grid;
		break;
	case G_FV_MASEK:		
		sub_len=512/4;		sft=2*grid;
		break;
	default:		
		ASSERT( 0 );
	break;
	}
	M = v_len/sub_len;			
	ASSERT( v_len%M==0 );
	base=nShift;
	for( cur=-nShift; cur<=nShift; cur++ )	{
		v_cur = V_shift+(cur+base)*v_len;
		off = -sft*cur+sub_len*8;
		for( r = 0; r < M; r++ )	{
			pos_0 = r*sub_len;
			for( c = 0; c < sub_len; c++ )	{
				temp = 0;
				for( j = 0; j < 8; j++ )	{
					pos_1 = (c*8+j-sft*cur+sub_len*8)%(sub_len*8);
					bit = pos_1%8;
					if( BIT_TEST( V1[pos_0+pos_1/8],BIT_FLAG_H>>bit ) )
						BIT_SET( temp,BIT_FLAG_H>>j );
				}
				v_cur[pos_0+c] =temp;
			}
		}
//		GIP_save_IplImage_I8( M,sub_len,GIP_path_(_T("fv_shift"),cur+base),v_cur,sub_len,0x0 );
	}

}

void _dump_8_( )	{
	int i,j,t;
	int no[9*64],nz[8];

	for( i = 0; i < 9; i++ )	nz[i]=0;
	for( i = 0; i < 0x100; i++ )	{
		t = 0;
		for( j=0; j < 8; j++ )	{
			t += BIT_TEST( i,BIT_FLAG_1<<j ) ? 1.0:0.0;
		}
		G_PRINTF( _T("%d,"),t );
		//	no[t*64+nz[t]]=i;
	//	nz[t]++;
	}

/*	for( i = 0; i < 9; i++ )	{
		for( j = 0; j < nz[i]; j++ )	{
			G_PRINTF( "case %d: ",no[i*64+j] );
			if( (j+1)%8==0 )	G_PRINTF( "\r\n" );
		}
		G_PRINTF( "hd+=%d;\r\nbreak;\r\n\r\n",i );
	}*/
	G_PRINTF( _T("\r\n") );
}

/*
	mapping of hamming distance

	v0.1	cys
		3/16/2009
*/
static int _HDD_MAP_[256]={
		0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,
		3,4,3,4,4,5,3,4,4,5,4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,
		3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,
		4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,
		3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,
		6,7,6,7,7,8
	};

/*	
	注意：
		1 V1包含多个shift之后的code，而V2只有一个code

	ISSUE-NEEDED:
		1 需进一步提高速度

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/25/2009	
*/
double _Hamming_distance_( int nShift,int v_len,G_U8 *V1,G_U8 *V2,G_U8* mask_1,G_U8* mask_2,int flag )	{
	double *a,hd=0.0,hd_min;
	int byte,bit,i,j,sft,c_1,c_2,base,pos_1,nNoiseBit;
	G_U8 temp,*v_cur,*mask_cur,noise;
	
//	_dump_8_( );

	base=nShift;
	hd_min=DBL_MAX;
	for( sft=-nShift; sft<=nShift; sft++ )	{
		v_cur = V1+(sft+base)*v_len;
		mask_cur = mask_1+(sft+base)*v_len;
		hd = 0.0;
		nNoiseBit = 0;
		for( i = 0; i < v_len; i++ )	{
//			noise = mask_cur[i] | mask_2[i];
//			nNoiseBit += bit_1[noise];
			temp = (V2[i]^v_cur[i]);
//			temp = temp&(~noise);
			hd += _HDD_MAP_[temp];
//			for( bit=0; bit < 8; bit++ )	{
//				hd += BIT_TEST( temp,BIT_FLAG_1<<bit ) ? 1.0:0.0;
//			}
		}
		hd = hd/(v_len*8-nNoiseBit);
		hd_min = MIN( hd,hd_min );
	}
/*
	a = GIP_calloc( sizeof(double),nShift*2+1 );
	for( i = 0;	i < v_len; i++ )	{
		byte = i;		
		for( j = 0; j < 8; j++ )	{
			c_1 = BIT_TEST( V1[byte],BIT_FLAG_H>>j );
			for( sft=-nShift; sft<=nShift; sft++ )	{
				pos_1 = (i*8+j+sft+v_len*8)%(v_len*8);
				bit = pos_1%8;
				c_2 = BIT_TEST( V2[pos_1/8],BIT_FLAG_H>>bit );
				a[sft+base] +=  c_1==c_2 ? 0.0:1.0;
			}
		}
	}

	hd_min=DBL_MAX;
	for( sft=-nShift; sft<=nShift; sft++ )	{
		hd = a[sft+base]/(v_len*8);
		hd_min = MIN( hd,hd_min );
	}

	GIP_free( a );*/

	return hd_min;
}

/*	
	基于hamming distance确定夹角

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/25/2009	
*/
int _Hamming_col_off_( int M,int N,G_U8 *V1,G_U8 *V2,G_U8 *temp,int flag )	{
	double *a,hd=0.0,hd_min,pos_min=-1;
	int byte,bit,i,j,sft,c_1,c_2,pos;
	G_U8 *v_col=temp,map;

	for( i = 0; i < M; i++ )	{
		pos=GIP_RC2POS( i,0,N );
		v_col[i] = V1[pos];
	}

	hd_min=DBL_MAX;
	for( j = 0; j < N; j++ )	{
		hd = 0.0;
		for( i = 0; i < M; i++ )	{
			pos=GIP_RC2POS( i,j,N );
			map = (V2[pos]^v_col[i]);
			hd += _HDD_MAP_[map];
		}
		hd = hd/(M*8);
		if( hd_min>hd )	{
			hd_min = hd;		pos_min=j;
		}
	}

	return pos_min;
}

/*	
	设定只计算left half part的_Hamming_distance_因此无需mask
	为了速度单独列出


	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/17/2009	
*/
double _Hamming_distance_cys_( int M,int N,int nShift,int v_len,G_U8 *V1,G_U8 *V2,int flag )	{
	double hd=0.0,hd_min;
	int i,j,sft,base,pos;
	G_U8 temp,*v_cur;
	
	base=nShift;
	hd_min=DBL_MAX;
	for( sft=-nShift; sft<=nShift; sft++ )	{
		v_cur = V1+(sft+base)*v_len;
		hd = 0.0;
		for( i = 0; i < M; i++ )	{
		for( j = 0; j < N/2; j++ )	{
			pos=GIP_RC2POS( i,j,N );
			temp = (V2[pos]^v_cur[pos]);
			hd += _HDD_MAP_[temp];
		}
		}
		hd = hd/(M*N*4);
/*		{			
//
//	教训：采用the product-of-sum (POS) of individual subfeature Hamming distances (HD).
//		速度慢，结果也莫名其妙
//		3/20/2009
		int byte,bit,b1,b2;
		double w;
		hd = 1.0;
		for( i = 0; i < N/2*8; i++ )	{
			w = 0.0;
			byte = i/8;		bit=i%8;
			for( j = 0; j < M; j++ )	{
				pos=GIP_RC2POS( j,byte,N );
				b1 = BIT_TEST( V2[pos],BIT_FLAG_H>>bit );
				b2 = BIT_TEST( v_cur[pos],BIT_FLAG_H>>bit );
				w+= (b1==b2 ? 0.0 : 1.0);
			}	
			w /= M;
			hd = hd*w;			
		}
		hd = pow( hd,1.0/4/N );
		}
*/
		hd_min = MIN (hd,hd_min )	;
	}

	return hd_min;
}

/*
	同时计算fv与fv[nEnry]的hamming距离，
	需要集中优化

	I_buffer[nEntry]

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/24/2009		
*/
void G_hamming_dist_( int nEntry,double *res,G_U8 *fv,G_U8 *mask,G_FVector_ *hfv,int *I_buffer,int flag )	{
	double hd=0.0;
	int v_len=hfv->V_len,bpSF=hfv->bpSF,k,nShift=hfv->nShift,ldN;
	int M,N,i,j,sft,base,pos,*m_pos=GIP_NULL,m_len;
	G_U8 temp,*V1,*V2;
	G_FV_TYPE type=hfv->type;

	m_len=0;
	switch( type )	{
	case G_FV_MASEK:		
		M=hfv->sf_M,		N=hfv->sf_N/4;
		ldN = N;							
/*		m_pos=GIP_alloc( sizeof(int)*M*N );
		M=hfv->sf_M,		N=hfv->sf_N/4;
		ldN = N;
		for( i = 0; i < M*N; i++ )	
			m_pos[i]=i;
		m_len=M*N;*/
		break;
	case G_FV_CYS:		
		M=hfv->sf_M,		N=hfv->sf_N/8;
		ldN = N;							
		break;
	case G_FV_CYS_HALF:		
		M=hfv->sf_M,		N=hfv->sf_N/8/2;
		ldN = N*2;							
		break;
	case G_FV_CYS_1:		
		M=hfv->sf_M,		N=hfv->sf_N/4;
		ldN = N;							
		break;
	case G_FV_DCT_CODE:		
		break;
	default:	
		ASSERT( 0 );
	}
	for( k = 0; k < nEntry; k++ )
		res[k] = DBL_MAX;

	base=nShift;
	for( sft=-nShift; sft<=nShift; sft++ )	{
		V1 = hfv->V+(sft+base)*v_len;
		for( k = 0; k < nEntry; k++ )	{
			V2 = fv+k*v_len;
			hd = 0.0;
			for( i = 0; i < M; i++ )	{
				pos=GIP_RC2POS( i,0,ldN );
				for( j = 0; j < N; j++ )	{
					temp = (V2[pos]^V1[pos]);
					hd += _HDD_MAP_[temp];
					pos++;
				}
			}			
/*			for( i = 0; i < m_len; i++ )	{
				pos=m_pos[i];
				temp = (V2[pos]^V1[pos]);
				hd += _HDD_MAP_[temp];
			}*/
			hd = hd/(M*N*8);
			ASSERT( hd>=0.0 && hd<=1.0 );
			res[k] = MIN (hd,res[k] )	;
//			res[k] += hd/(nShift*2+1);
		}
	}
	
//	GIP_free( m_pos );
	return;
}

/*
	注意：
		1 hfv记录了fv的各种变化
		2 v_sft与fv的区别
	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/7/2009		
*/
double G_fv_distance( G_U8 *fv,G_U8 *mask,G_FVector_ *hfv,int pos_1,int pos_2,char *temp )	{
	double dist,*dv,*base,*v_1,d,ang;
	int V_len=hfv->V_len,bpSF=hfv->bpSF,k,nShift=hfv->nShift;
	G_FV_TYPE type=hfv->type;
	G_U8 *v_sft,*m_sft;

	bpSF = 3;
	switch( type )	{
	case G_FV_MASEK:		//GIP_DISTANCE_HAMMING:
		if( 0 )	{			//!!testing!!		
			nShift=2;
			V_len=1;
			fv[pos_1*V_len]=0xa0;		fv[pos_2*V_len]=0x0b;			
		}
		v_sft = hfv->V;			m_sft = hfv->mask;
//		_code_shift_( nShift,2,V_len,fv+pos_1*V_len,v_sft,0x0 );
		dist = _Hamming_distance_( nShift,V_len,v_sft,fv+pos_2*V_len,m_sft,mask+pos_2*V_len,0x0 );
		break;
	case G_FV_CYS:		//GIP_DISTANCE_HAMMING:
		v_sft = hfv->V;			m_sft = hfv->mask;
//		ang = _Hamming_col_off_( hfv->sf_M,hfv->sf_N/4,v_sft+nShift*V_len,fv+pos_2*V_len,temp,0x0 );
		dist = _Hamming_distance_cys_( hfv->sf_M,hfv->sf_N/4,nShift,V_len,v_sft,fv+pos_2*V_len,0x0 );
//		dist = _Hamming_distance_cys_( hfv->sf_M,hfv->sf_N/8,nShift,V_len,v_sft,fv+pos_2*V_len,0x0 );
		break;
	case G_FV_DCT_CODE:		//GIP_DISTANCE_HAMMING:
		dist = _Hamming_distance_1_( V_len,fv+pos_1*V_len,fv+pos_2*V_len,hfv->bpSF,hfv->nzSF,temp );
		break;
	case G_FV_COEFFICIENT:		//GIP_DISTANCE_HAMMING:
		dv = (double*)(fv);
		dist = _euclid_distance_( V_len,dv+pos_1*V_len,dv+pos_2*V_len,0x0 );
		break;
	default:	{
		int p_1,pos,r,c,c_1,M,N,i;
		M=hfv->sf_M;		N=hfv->sf_N;
		dv = (double*)(fv);
		base = GIP_alloc( sizeof(double)*V_len*2 );		v_1=base+V_len;
		for( i = 0 ; i < V_len; i++ )		base[i]=dv[pos_1*V_len+i];
		dist = DBL_MAX;
		for( k=0; k<=0; k++ )	{
//		for( k=-5; k<=5; k++ )	{
			for( r=0; r<M; r++ )	{
			for( c=0; c<N; c++ )	{
				c_1 = (c+k+N)%N;
				p_1=GIP_RC2POS( r,c_1,N );			pos=GIP_RC2POS( r,c,N );
				v_1[p_1] = base[pos];
			}
			}
			d = _euclid_distance_( V_len,v_1,dv+pos_2*V_len,0x0 );
			dist = MIN( dist,d );
		}
		for( i = 0 ; i < V_len; i++ )		dv[pos_1*V_len+i]=base[i];
		GIP_free( base );
				}
		break;
	}

//	GIP_free( v_sft );

	return dist;
}

/*
	注意:
		1 hBase指向GIP_LIB*

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/7/2009		
*/
void G_DCRIMI_set( void *hBase,G_DCRIMI_ *hDcrimi,int D_span,float *D_inter,float *D_intra,double crr,int flag )	{
	int i,grid=0;
	double mean_a,mean_r,devian_a,devian_r,nz_a,nz_r,s,D_s;
	double f_ar,f_rr,w_a,w_r,f_ar_g=1.0e-7;
	double *f_ar_8=hDcrimi->f_ar_8,*f_rr_8=hDcrimi->f_rr_8,*hd_8=hDcrimi->hd_8;

	GIP_MEMCLEAR( hDcrimi,sizeof(G_DCRIMI_) );
	if( D_inter==GIP_NULL || D_intra==GIP_NULL )
		return;

	for( i=0; i<8; i++ )	{
		f_ar_8[i]=-1.0;			f_rr_8[i]=-1.0;				hd_8[i]=-1.0;
	}
	hDcrimi->CRR = crr;
	hDcrimi->hBase = hBase;
	D_s =1.0/D_span;		//only for hamming distance

	mean_a = 0.0;		mean_r=0;
	nz_a=0.0;			nz_r=0.0;
	for( i = 0; i <= D_span; i++ )	{
		mean_a += i*D_intra[i];
		nz_a += D_intra[i];
		mean_r += i*D_inter[i];
		nz_r += D_inter[i];

	}
	hDcrimi->nz_a = nz_a;				hDcrimi->nz_r = nz_r;
	mean_a = nz_a==0 ? 0.0 : mean_a/nz_a*D_s;
	mean_r = nz_r==0 ? 0.0 : mean_r/nz_r*D_s;
	hDcrimi->mean_a = mean_a;				hDcrimi->mean_r = mean_r;

	devian_a = 0.0;		devian_r=0;
	w_a=0.0;			w_r=0.0;
	for( i = 0; i <= D_span; i++ )	{
		w_a+=D_intra[i];		w_r+=D_inter[i];
		f_ar = w_r*1.0/nz_r;
		f_rr = (nz_a-w_a)*1.0/nz_a;
		while( f_ar>=f_ar_g && f_ar_8[grid]==-1.0 )	{
			f_ar_8[grid]=f_ar;
			f_rr_8[grid]=f_rr;				hd_8[grid]=i*1.0/D_span;
			f_ar_g*=10;		grid++;
			if( f_ar < f_ar_g )
				break;
		}
		if( f_ar>=f_rr && hDcrimi->rEER==0.0 )	{
			hDcrimi->rEER = (f_ar+f_rr)/2.0;
			hDcrimi->eer_sep = hDcrimi->sep;
		}
		if( f_ar>1.0e-4 && hDcrimi->sep==0.0 )	{
			hDcrimi->sep = i*1.0/D_span;
			hDcrimi->rFAR = f_ar;
			hDcrimi->rFRR = f_rr;
		}
		if( D_intra[i]!=0 )	{
			s = (i*D_s-mean_a);
			devian_a += s*s*D_intra[i];
		}
		s = (i*D_s-mean_r);
		devian_r += s*s*D_inter[i];	
	}
	devian_a = nz_a==0 ? 0.0 : sqrt(devian_a/nz_a);
	devian_r = nz_r==0 ? 0.0 : sqrt(devian_r/nz_r);
	hDcrimi->devia_a = devian_a;				hDcrimi->devia_r = devian_r;

	s=sqrt( (devian_a*devian_a)+(devian_r*devian_r) )/2.0;
	hDcrimi->D_ = s==0 ? 0.0 : fabs(mean_a-mean_r)/s;
}

/*
	sep改变，重新计算far,frr

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/12/2009		
*/
void G_DCRIMI_update_sep( G_DCRIMI_ *hDcrimi,int D_span,float *D_inter,float *D_intra,int flag )	{
	int i,pos;
	double w_a,w_r,D_s =1.0/D_span;

	pos = G_DOUBLE2INT(hDcrimi->sep/D_s);
	w_a=0.0;			w_r=0.0;
	for( i = 0; i <= pos; i++ )	{
		w_a+=D_intra[i];		w_r+=D_inter[i];
	}
	hDcrimi->rFAR = w_r*1.0/hDcrimi->nz_r;
	hDcrimi->rFRR = (hDcrimi->nz_a-w_a)*1.0/hDcrimi->nz_a;
}


/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/7/2009		
*/
void G_FV_info( G_FVector_ *hfv,int flag )	{
	G_PRINTF( _T("------ Feature Vector Infomation ------\r\n") );
	if( hfv==GIP_NULL )	{
		G_PRINTF( _T("!!!hfv=GIP_NULL!!!") );
		return;
	}

	switch( hfv->type )	{
	case G_FV_ROI:
		G_PRINTF( _T(" type=ROI\r\n") );
		break;
	case G_FV_MALI:
		G_PRINTF( _T(" type=MALI\r\n") );
		break;
	case G_FV_MASEK:
		G_PRINTF( _T(" type=MASEK\r\n") );
		break;
	case G_FV_CYS:
		G_PRINTF( _T(" type=CYS\r\n") );
		break;
	case G_FV_CYS_HALF:
		G_PRINTF( _T(" type=CYS_half\r\n") );
		break;
	case G_FV_CYS_1:
		G_PRINTF( _T(" type=CYS_1\r\n") );
		break;
	case G_FV_COEFFICIENT:
		G_PRINTF( _T(" type=COEFFICIENT\r\n") );
		break;
	default:
		G_PRINTF( _T(" type=UNKNOWN\r\n") );
	}	

	G_PRINTF( _T(" f_M=%d,f_N=%d,V_len=%d\r\n"),hfv->sf_M,hfv->sf_N,hfv->V_len );
	
}

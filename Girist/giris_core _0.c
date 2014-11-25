/*
	包含部分giris_core.c的代码
*/
#include <malloc.h>
#include <memory.h>
#include <time.h>
#include <FLOAT.h>
#include <math.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/gip_util.h"
#include "../util/gip_imag_transform.h"
#include "../util/gip_imag_property.h"
#include "../util/gip_fv.h"
#include "../util/gip_thread.h"
#include "../util/gip_library.h"
#include "../util/Compact_DataStruc.h"
#include "../d_ft_wt/gip_dwt.h"
#include "../d_ft_wt/gip_gabor.h"
#include "giris_core.h"


/*
	based on Hough Transform
	注意:
		1 hBinary为二进制图

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/4/2008	
*/
double _circle_fit( GIP_OBJECT *hBinary,double *R,double *c_r,double *c_c,int flag )	{
	int i,r,c,pos,M=hBinary->M,N=hBinary->N,*hough,hough_limit,h_max;
	int r_0,r_1,c_0,c_1,hc_0,hc_1,hr_0,hr_1,hR_1,hR_0,x_R,x_c,x_r,h_r,h_c,h_R;
	double w,*data=hBinary->data,a,da,err=0.0;
	
	ASSERT( M>2 && N>2 );

	r_0=M*2;	r_1=-1;
	c_0=N*2;	c_1=-1;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,N );
		if( data[pos]==1.0 )		continue;
		r_0=MIN( r_0,r);		r_1=MAX( r_1,r);		
		c_0=MIN( c_0,c);		c_1=MAX( c_1,c);		
	}
	}
//	ASSERT( r_0>0 && r_1<M-1 && c_0>0 && c_1<N-1 );
	if( BIT_TEST( flag,IRIS_CENTER_FIX_) )	{
		hR_1=MIN(r_1-r_0,c_1-c_0)/2+1;			hR_0=hR_1/2-1;
		hr_0 = *c_r;							hr_1 = *c_r+1;
		hc_0 = *c_c;							hc_1 = *c_c+1;
	}else	{
		hR_1=MIN(r_1-r_0,c_1-c_0)/2+1;			hR_0=hR_1/2-1;
		hr_0 = (r_1+r_0)/2-hR_0;				hr_1 = (r_1+r_0)/2+hR_0;
		hc_0 = (c_1+c_0)/2-hR_0;				hc_1 = (c_1+c_0)/2+hR_0;
	}

	x_R=hR_1-hR_0+1;		x_r=(hr_1-hr_0+1);			x_c=(hc_1-hc_0+1);
	hough_limit=x_R*x_r*x_c;
	hough=GIP_alloc( sizeof(int)*hough_limit );
	memset( hough,0x0,sizeof(int)*hough_limit );

	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos = GIP_RC2POS( r,c,N );
		w = data[pos]+data[pos-1]+data[pos+1]+data[pos-N]+data[pos+N];		 
		if( w==0.0 || w==5.0 )	continue;
//		_hough_3( r,c,M,N,hough );
		for( h_r=hr_0; h_r<hr_1; h_r++ )	{
		for( h_c=hc_0; h_c<hc_1; h_c++ )	{
			h_R=(int)( sqrt( (h_r-r)*(h_r-r)+(h_c-c)*(h_c-c)+0.5 ) );
//			ASSERT( hR_0<=h_R && h_R<=hR_1 );
			if( hR_0>h_R || h_R>hR_1 )	continue;
			pos=(h_R-hR_0)*x_r*x_c+(h_r-hr_0)*x_c+h_c-hc_0;
			hough[pos]++;			
		}
		}
	}
	}
	h_max=0;
	for( i = 0; i < hough_limit; i++ )	{
		if( hough[i]>h_max )	{
			h_max=hough[i];		pos=i;
		}
	}
	h_R=pos/(x_r*x_c);		h_r=(pos-h_R*x_r*x_c)/x_c;		h_c=pos-h_R*x_r*x_c-h_r*x_c;
	*R=h_R+hR_0;
	*c_r=h_r+hr_0;
	*c_c=h_c+hc_0;

	err=0.0;
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		a=sqrt( (h_r-r)*(h_r-r)+(h_c-c)*(h_c-c) );
		pos = GIP_RC2POS( r,c,N );
		if( a <= h_R )	{
			if( data[pos]==0.0 )	err+=1.0;
			data[pos]=1.0;
		}else	{
			if( data[pos]==1.0 )	err+=1.0;
			data[pos]=0.0;
		}
	}
	}
	GIP_free( hough );

	return err;
}

/*
	基于pupil的(r,c,radius)

	Copyright 2008-present, Grusoft
	v0.1	cys
		12/5/2008	
*/
void _IRIS_find_iris_( G_IRIS_ *hIris,int method )	{
	GIP_OBJECT *hObj=hIris->hObj,iris,*hOI=&iris,binary,*hBinary=&binary;
	int i,r,c,pos,M=hObj->M,N=hObj->N,M_p,N_p,n_1,n_2,nPupil;
	double r_0,c_0,r_min,c_min,a,*data,u_1,u_2,th_1,thresh,err,R;

	r_0=hIris->p_r;						c_0=hIris->p_c;
	hIris->s_r = r_0;					hIris->s_c = c_0;
	M_p=MIN(hIris->p_r,M-hIris->p_r)*1.9;
	N_p=MIN(hIris->p_c,N-hIris->p_c)*1.8;
//	M_p=MAX(2*hIris->p_radius,R);		N_p=MAX(2*hIris->p_radius,R);	
	GIP_OBJ_init( hOI,M_p,N_p,0x0 );
	r=(int)(hIris->p_r-M_p/2),			c=(int)(hIris->p_c-N_p/2);
	GIP_OBJ_sub( hObj,hOI,r,c,M_p,N_p,0x0 );
	
	GIP_OBJ_init( hBinary,M_p,N_p,0x0 );
//	GIP_canny( hOI,hBinary,0x0 );
	data = hBinary->data;
	R=hIris->p_radius;
	r_0=M_p/2;				c_0=N_p/2;
	for( r = 0; r < hOI->M; r++ )	{		
	for( c = 0; c < hOI->N; c++ )	{
		a = sqrt( (r_0-r)*(r_0-r)+(c_0-c)*(c_0-c) );
		pos=GIP_RC2POS( r,c,hOI->N ); 
		if( a<R*2.5 )	
		{	data[pos]=0.0;		continue;	}
		data[pos] = data[pos]>0 ? 1.0:0.0;
	}
	}
	err = _circle_fit( hBinary,&R,&r_0,&c_0,IRIS_CENTER_FIX_ );
	ASSERT( R/hIris->p_radius >= 2.5 );
	r_0 = r_0+hIris->p_r-M_p/2;
	c_0 = c_0+hIris->p_c-N_p/2;
	hIris->s_r = r_0;		hIris->s_c = c_0; 
	hIris->s_radius = R;

	data = hObj->data;
	for( r = 0; r < hObj->M; r++ )	{		
	for( c = 0; c < hObj->N; c++ )	{
		a=sqrt( (r_0-r)*(r_0-r)+(c_0-c)*(c_0-c) );
		pos=GIP_RC2POS( r,c,hObj->N ); 
		if( fabs(a-R)<1.0 )	
		{	data[pos]=0.0;		}
		a=sqrt( (hIris->p_r-r)*(hIris->p_r-r)+(hIris->p_c-c)*(hIris->p_c-c) );
		if( fabs(a-hIris->p_radius)<1.0 )	
		{	data[pos]=0.0;		}
	}
	}
	GIP_save_IplImage_d( hObj->M,hObj->N,GIP_path_(_T("iris"),M*N),hObj->data,hObj->N,0x0 );
}

/*
	get the center coordinates of the pupil

	效率 ***
		need fast GIP_CHT
	健壮 -			
		仅针对iris寻找pupil的情况，需要完善各种情况
	稳定 **			
		基于阈值迭代总能收敛到(直方图)波谷的假设
		基于pupil占据一定的面积(>0.05)

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/9/2006	
*/	
int _IRIS_find_pupil_1_( G_IRIS_ *hIris,int method )	{
	GIP_OBJECT *hObj=hIris->hObj;
	int i,r,c,pos,M=hObj->M,N=hObj->N,MN=M*N,n_1,n_2,nz,ret=0,nLevel=0,*arrN1,*arrPos;
	GIP_FLOATING *data = hObj->data;
	double u_1,u_2,a,thresh,th_1,epsi,I_0,I_1,*arrThresh,R,r_0,c_0;
	GIa_HANDLE *hIa=hIris->hIa;

	I_0=0.0;		I_1=hIa->mean;
	arrThresh = hIa->D_buffer;
	arrN1 = hIa->I_buffer;			arrPos=arrN1+MN;
	arrN1[0] = MN;
	while( 1 )	{
		thresh = 0.0;	nz=0;	//背景与对象的面积相近时，可取thresh为平均值 
		for( i = 0; i < MN; i++ )	{
			a = data[i];		
			if( a < I_0 || a > I_1 )		continue;
			thresh+= data[i];		nz++;
		}
		arrN1[nLevel] = arrN1[nLevel]-nz;
		arrN1[nLevel+1] = nz;			//寄存器
		if( nz*1.0/MN<0.05 )
			break;
		thresh /= nz;
		epsi = fabs(thresh*0.01);
		while( 1 )	{
			u_1=0.0;		u_2=0.0;
			n_1=0;			n_2=0;
			for( i = 0; i < MN; i++ )	{
				a = data[i];
				if( a < I_0 || a > I_1 )		continue;
	//			a = GIP_LAPLAS_( i,M_p,N_p,N_p,hObj->data,0x0 );
				if( a>thresh )
				{	u_1 += a;		n_1++;	}
				else
				{	u_2 += a;		n_2++;	}
			}
			ASSERT( n_1+n_2==nz );
			th_1 = (u_1/n_1+u_2/n_2)/2.0;
			if( fabs(th_1-thresh)<=epsi )
				break;
			thresh = th_1;
		}
		I_1 = thresh;	
		arrThresh[nLevel++] = thresh;
	}
	for( i = nLevel; i >= 0; i-- )	{
		if( arrN1[i]>MN*0.05 )
		{	n_1=i;		break;		}
	}
	n_1=2;
	thresh = arrThresh[n_1-1];		
	nz=0;
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos = GIP_RC2POS( r,c,N );
		a = data[pos];		
		if( a > thresh )		continue;
		if( data[pos-1]>thresh && data[pos+1]>thresh && data[pos-N]>thresh && data[pos+N]>thresh )
			continue;
		if( data[pos-1]>thresh || data[pos+1]>thresh || data[pos-N]>thresh || data[pos+N]>thresh 
		|| data[pos-N-1]>thresh || data[pos-N+1]>thresh || data[pos+N-1]>thresh || data[pos+N+1]>thresh )	{
			arrPos[nz++]=pos;	
		}
	}
	}
	
	GIP_save_IplImage_d( hObj->M,hObj->N,_T(".\\trace\\iris_thrsh_1.jpg"),hObj->data,hObj->N,0x0 );
	GIa_CHT( hIris->hObj,nz,arrPos,&R,&r_0,&c_0,0x0 );
	hIris->p_r = r_0;		hIris->p_c = c_0; 
	hIris->p_radius = R;

	return ret;
}

/*	
	using linear hough transform to get the eyelid line

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/17/2009	
*/	
void _eyelid_line_( int M,int N,int *mask,GIP_FLOATING *e_data,int flag )	{
	int r,c,pos,dgr_0,dgr_1,R_0,R_1,d,R,s_max,no_max=-1;
	int *lht;
	double x_1,xita,rou;

	dgr_0=-5;				dgr_1=5;		
	x_1 = dgr_1*PI/180;
//	R_1=MAX( M*sin(x_1)+N*cos(x_1),0.0 );
	R_1=M+N;
	R_0=-R_1;					
	lht=GIP_calloc( sizeof(int),(dgr_1-dgr_0)*(R_1-R_0) );
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		if( mask[pos]!=1 )	continue;
		for( d=dgr_0; d<dgr_1; d++ )	{
			xita = d*PI/180;
			rou = r*cos(xita)+c*sin(xita);
			ASSERT( rou>=R_0 && rou<=R_1 );
			R=G_DOUBLE2INT( rou );
			pos = (d-dgr_0)*(R_1-R_0)+R-R_0;
			lht[pos]++;
		}
	}
	}

	s_max=0;		no_max=-1;
	for( d=(dgr_1-dgr_0)*(R_1-R_0)-1; d>=0; d--)	{
		if( lht[d]>s_max )
		{	s_max=lht[d];		no_max=d;		}
	}
	d=no_max/(R_1-R_0);			R=no_max%(R_1-R_0);
	d+=dgr_0;					R+=R_0;

	if( 1 )	{		//print line to the edgemap
		xita = d*PI/180;
		x_1 =1.0/cos(xita);
		for( c = N*0.05; c < N*0.95; c++ )	{
			if( c%3==0 )	continue;
			r = (R-c*sin(xita))*x_1;
			if( r<0 || r>= M )		continue;
			pos=GIP_RC2POS( r,c,N );
			e_data[pos]=0.0;
		}
	}

	GIP_free( lht );
}

/*
	注意：
	1 在确定pupil之后调用，以帮助定位
	2 暂不考虑rotation

	[REF]: Comparison of eyelid and eyelash detection algorithms for performance improvement of iris recognition

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/16/2009		
*/
void _IRIS_find_eyelid_( G_IRIS_ *hIris,int M,int N,GIP_FLOATING *data,double *g_I,double *g_rad,int *mask,G_IRIS_REGION region,int flag )	{
	GIP_OBJECT *hEdge;
	int MN=M*N,I_0,I_1,r,c,pos,e_nz,isNewMask=0;
	int h,k,h_0,h_1,k_0,k_1,*vote,a_no,v_max,v_no,a_span;
	GIP_FLOATING *e_data;
	double g_t1=10.0,g_t0,w,a,a_0,a_1,p_R=hIris->p_radius,a_step,xita,q;
//	float *g_I,*g_rad;

	isNewMask = mask==GIP_NULL;
	hEdge = GIP_OBJ_create( M,N,data,0x0 );
	e_data = hEdge->data;
	q=1.0;

	if( isNewMask==1 )	{
		ASSERT( region==IRIS_NOISE_EYELID );
		ASSERT( g_I==GIP_NULL );						ASSERT( g_rad==GIP_NULL );		
		g_I = GIP_alloc( sizeof(double)*MN*2 );			g_rad=g_I+MN;
		mask=GIP_alloc( sizeof(int)*MN );				GIP_MEMCLEAR( mask,sizeof(int)*MN );
		I_0=100;			I_1=255;
		g_t1=15.0,			g_t0=g_t1/4;
//		GIP_Smooth( hEdge,GIP_SMOOTH_MEDIAN,hIris->hIa->D_buffer,q,0x0 );
		GIP_Smooth( hEdge,GIP_SMOOTH_ANISO,hIris->hIa->D_buffer,q,0x0 );
		GIP_canny_( hEdge,I_1,I_0,g_t0,g_t1,mask,g_I,g_rad,G_MASK_3 );
	}
	h_0=N*0.3;		h_1=N*0.7;		
	k_0=M*0.2;		k_1=M*0.9;
	a_0=250;		a_span=40;		a_step=5;	a_1=a_0+a_step*a_span;		
	if( region==IRIS_UPPER_EYELID )	{		//(c-h)*(c-h)=a*(r-k);
		h_0 = hIris->p_c-p_R/2;		h_1 = hIris->p_c+p_R/2;
		k_0=M*0.1;								k_1=MIN( M*0.9,hIris->p_r );
	}else if( region==IRIS_LOWER_EYELID )	{
		h_0 = hIris->p_c-p_R/2;		h_1 = hIris->p_c+p_R/2;
		k_0=hIris->p_r;						k_1=M*0.9;
		a_0=-600;		a_span=40;			a_1=a_0+a_step*a_span;		
	}else if( region==IRIS_NOISE_EYELID )	{
		h_0 = N*0.7;					h_1 = N*0.8;
		k_0 = -M*0.5;					k_1=M;
		a_0=100;		a_span=40;			a_1=a_0+a_step*a_span;	
		hIris->nzNoise = 0;
	}
	vote = GIP_calloc( sizeof(int),(h_1-h_0)*(k_1-k_0)*a_span );
	e_nz = 0;
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		if( mask[pos]!=1 )	continue;
	//	if( fabs((r-hIris->p_r)*(r-hIris->p_r)+(c-hIris->p_c)*(c-hIris->p_c)-p_R*p_R)<2.0 )
	//		continue;
	//	if( c<hIris->p_c-p_R*2 || c>hIris->p_c+p_R*2 )
	//		continue;
		e_data[pos]=255.0;
		e_nz++;
		for( h=h_0; h<h_1; h++ )	{
		for( k=k_0; k<k_1; k++ )	{
			if( r<=k )		continue;
			a = (c-h)*(c-h)/(r-k);
			if( a<a_0 || a>=a_1 )		continue;
			xita = atan2( -a,2*(c-h) );
			if( fabs(xita-g_rad[pos]) > 0.50 ) 
				continue;
			a_no = (a-a_0)/a_step;
//			pos_1 = ((h-h_0)*(k_1-k_0)+k-k_0)*a_span+a_no;
			vote[((h-h_0)*(k_1-k_0)+k-k_0)*a_span+a_no]++;
		}
		}
		
OUTER_LOOP:	;	
	}
	}
//		G_PRINTF( "e_nz=%d\r\n",e_nz );
	v_max = 0;			v_no=-1;
	for( k = (h_1-h_0)*(k_1-k_0)*a_span-1; k >= 0; k-- )	{
		if( vote[k]>v_max )	{
			v_max = vote[k];		v_no=k;
		}
	}
	if( v_max<6 )	{
	}else	{
		a_no = v_no%a_span;				pos=v_no/a_span;
		h=pos/(k_1-k_0)+h_0;			k=pos%(k_1-k_0)+k_0;
		a = a_0+a_no*a_step;
		hIris->ul_a = a;			hIris->ul_h=h;			hIris->ul_k=k;		
	}
	e_nz = 0;
	if( region==IRIS_NOISE_EYELID && hIris->ul_a!=0.0 )	{
		for( c=0; c<N; c++ )	{
			w=(c-h)*(c-h)*1.0/a+k;
			r = G_DOUBLE2INT( w );
			if( r < 0 )	{
				hIris->nzNoise += M;
			}else if( r>=M )	{		
			} else	{
				e_nz++;
				hIris->nzNoise += (M-r);
			}
		}
		if( v_max<e_nz/4 )	{		//verify
			hIris->ul_a = 0.0;			hIris->ul_h=0.0;			hIris->ul_k=0.0;
			hIris->nzNoise = 0.0;
		}
	}
	
	GIP_Parabolic( hEdge,hIris->ul_a,hIris->ul_h,hIris->ul_k,GIP_IMAG_ORIGIN );
	GIP_save_IplImage_d( hEdge->M,hEdge->N,_T(".\\trace\\roi_edge_1.jpg"),hEdge->data,hEdge->N,0x0 );
	if( isNewMask==1 )		{
		GIP_free( mask );			
		GIP_free( g_I );			
	}
	GIP_free( vote );
	GIP_OBJ_clear(hEdge);		GIP_free( hEdge );

}
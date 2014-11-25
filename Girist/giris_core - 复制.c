/*
ISSUE-NEEDED:
	1 采用contour integral确定半径。
	2 "Unsupervised image segmentation using a simpleMRFmodel with a newimplementation scheme"
	3 "A Study of Segmentation and Normalization for Iris Recognition Systems"
	4 "LOWER ORDER CIRCLE AND ELLIPSE HOUGH TRANSFORM"
	  "CIRCULAR OBJECT DETECTION USING A MODIFIED HOUGH TRANSFORM"
	5 Suppressing eyelash interference by dilation operator
		"An effective iris location method with high robustness"
	6 "Rotation Compensated Human Iris Matching"
	7 "ROBUST IRIS FEATURE EXTRACTION AND MATCHING"
	8 "AN IRIS RECOGNITION SYSTEM USING PHASE-BASED IMAGE MATCHING"
	9 " A Review of Iris Recognition Algorithms "

*/
#include <malloc.h>
#include <memory.h>
#include <time.h>
#include <FLOAT.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/gip_util.h"
#include "../util/gip_imag_transform.h"
#include "../util/gip_imag_property.h"
#include "../util/gip_fv.h"
#include "../util/gip_library.h"
#include "../d_ft_wt/gip_dwt.h"
#include "../d_ft_wt/gip_gabor.h"
#include "giris_core.h"
#include "highgui.h"

/*
	针对不同的库暂采用不同的参数
	(I)	MMU iris database
		gs_Iris_G_0=8.0	
		g_II_s0=0.1,	g_II_s1=0.6
	(II) CNKI iris database
		gs_Iris_G_0=5.0	
		g_II_s0=0.1,	g_II_s1=0.6
*/
static double _ToG_= 5.0;					//thrshold of gradient for finding iris 
static double _ToI_0=0.1,_ToI_1=0.7;			//thrshold of intensity for finding iris 

//static GIa_HANDLE g_s_hIa;
static GIP_GABOR gabor_1,gabor_2;
static int M_normal=48,N_normal=512,g_nzSample;

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/1/2008	
*/
void G_IRIS_clear( G_IRIS_ *hIris )	{
	if( hIris == GIP_NULL )	return;

//	if( hIris->temp != GIP_NULL )	
//	{	GIP_free( hIris->temp );	hIris->temp=GIP_NULL; }
	if( hIris->hObj != GIP_NULL )		{	
		GIP_OBJ_clear( hIris->hObj );	
		GIP_free( hIris->hObj );	hIris->hObj=GIP_NULL;

		GIP_OBJ_clear( &(hIris->normal) );	
	}
	if( hIris->hMask != GIP_NULL )		{	
		GIP_OBJ_clear( hIris->hMask );	
		GIP_free( hIris->hMask );	hIris->hMask=GIP_NULL;
	}

	if( hIris->hIa != GIP_NULL )		{	
		GIa_HANDLE_clear( hIris->hIa );	
		GIP_free( hIris->hIa );		hIris->hIa=GIP_NULL; 
	}
}

/*
	注意：
		1 需要simple filtering

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/1/2008	
	v0.2	cys		
		2/3/2009	
*/
void G_IRIS_init_( GIP_IMAGE *hIMG,G_IRIS_ *hIris )	{
	int M,N,MN,ret=GIP_OK,t_len;
	GIP_FLOATING *val=GIP_NULL;
	double beta;
	GIP_OBJECT *hObj=GIP_NULL;
	GIa_HANDLE *hIa;

	memset( hIris,0x0,sizeof(G_IRIS_) );
	hIris->hObj=GIP_alloc( sizeof(GIP_OBJECT) );			hObj=hIris->hObj;
	GIP_import_IplImage( hIMG,&hObj->M,&hObj->N,&hObj->data,0x0 );
	M=hObj->M;			N=hObj->N;
	MN=M*N;
	t_len=MAX((M+2)*(N+2)*6,N);

	hIris->hMask=GIP_OBJ_create( hObj->M,hObj->N,GIP_NULL,0x0 );

	hIris->hIa=GIP_alloc( sizeof(GIa_HANDLE) );		//&(g_s_hIa);		
	hIa=hIris->hIa;
	GIa_HANDLE_init( hIa,hObj,0x0 );

//	GIP_Histo_Equal( hObj,256,0x0,GIP_NULL,0x0 );
//	GIP_Contrast_Enhance( hObj,256,0x0,GIP_NULL,0x0 );

//	hIris->alg = CGM_SS_DIRECT | CGM_A_AT | CGM_METHOD_MEDIAN;				//CGM_SS_ITER;	CGM_METHOD_DUAL
//	GIP_OBJ_init( &(hIris->normal),64,512,0x0 );
	GIP_OBJ_init( &(hIris->normal),48,512,0x0 );							//80->48

//	printf( "  {_init_}: M=%d,N=%d,mean=%g,devia=%g,contrast=%g\r\n",hObj->M,hObj->N,
//		hIa->mean,hIa->deviation,hIa->contrast );

GIP_exit:
	;
}

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
	char sTitle[80];

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
	sprintf( sTitle,".\\trace\\iris_%d_1.jpg",M*N );	
	GIP_save_IplImage_d( hObj->M,hObj->N,sTitle,hObj->data,hObj->N,0x0 );
	
/*	if( err < err_min )	{	
	}*/

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
	
	GIP_save_IplImage_d( hObj->M,hObj->N,".\\trace\\iris_thrsh_1.jpg",hObj->data,hObj->N,0x0 );
	GIa_CHT( hIris->hObj,nz,arrPos,&R,&r_0,&c_0,0x0 );
	hIris->p_r = r_0;		hIris->p_c = c_0; 
	hIris->p_radius = R;

	return ret;
}

/*
	for _IRIS_find_iris_2_

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/8/2009	
*/	
void _Arc_template_( int M,int N,double R_0,double R_1,int **A_ptr_,int **A_off_,int **A_R_,int *temp,int flag )	{
	GIP_OBJECT *hArc;
	int A_nz,*A_ptr,*A_off,*A_R,i,pos,dgr,MN=M*N,DGREE=361,r,c;
	GIP_FLOATING *A_data;
	double g,gr,gc;

	hArc=GIP_OBJ_create( M,N,GIP_NULL,0x0 );
	A_nz = 0;
	A_data=hArc->data;
	for( i = 0; i < DGREE; i++ )			temp[i]=0;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,N );
		gr = r-M/2.0;			gc = c-N/2.0;
		g = G_DOUBLE2INT( sqrt( gr*gr+gc*gc ) );
		if( g<R_0 || g>=R_1 )	
			A_data[pos]=1.0;
		else	{
			A_data[pos]= 0.0;
			A_nz++;
			dgr = G_RADIEN2DEGREE( atan2(gr,gc) );
			temp[dgr]++;
		}
	}
	}
//	GIP_save_IplImage_d( hArc->M,hArc->N,".\\trace\\ring.jpg",hArc->data,hArc->N,0x0 );

	A_ptr=GIP_alloc( sizeof(int)*(DGREE+1+A_nz*2) );		A_off=A_ptr+DGREE+1;		A_R=A_off+A_nz;
	*A_ptr_=A_ptr;			*A_off_=A_off;				*A_R_=A_R;
	A_ptr[DGREE]=A_nz;
	for( i = DGREE-1; i >= 0; i-- )	{
		A_ptr[i]=A_ptr[i+1]-temp[i];		temp[i]=A_ptr[i];
	}
	ASSERT( A_ptr[0]==0 );
	pos=GIP_RC2POS( M/2,N/2,N );
	for( i = 0; i < MN; i++ )	{
		if( A_data[i]==1.0 )		continue;
		gr=GIP_POS2R(i,N)-M/2.0;		gc=GIP_POS2C(i,N)-N/2.0;	
		dgr = G_RADIEN2DEGREE( atan2(gr,gc) );		
		A_R[temp[dgr]] = G_DOUBLE2INT( sqrt(gr*gr+gc*gc)-R_0 );
		A_off[temp[dgr]] = i-pos;
		temp[dgr]++;
	}
	for( i = 0; i < DGREE; i++ )			ASSERT( temp[i]=A_ptr[i+1] );

	GIP_OBJ_clear( hArc );				GIP_free( hArc );
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
	int MN=M*N,I_0,I_1,r,c,pos,e_nz,dgr,adj[8],isNewMask=0;
	int h,k,h_0,h_1,k_0,k_1,*vote,a_no,v_max,v_no,a_span;
	GIP_FLOATING *e_data;
	double gr,gc,g_t1=10.0,g_t0,w,a,a_0,a_1,p_R=hIris->p_radius,a_step,xita;
//	float *g_I,*g_rad;

	isNewMask = mask==GIP_NULL;
	g_t1=40.0,		g_t0=g_t1/4;
	I_0=0;			I_1=50;
	hEdge = GIP_OBJ_create( M,N,data,0x0 );
	e_data = hEdge->data;
//	GIP_Smooth( hEdge,GIP_SMOOTH_MEDIAN,hIris->hIa->D_buffer,0x0 );
//	GIP_Smooth( hEdge,GIP_SMOOTH_ANISO,hIris->hIa->D_buffer,0x0 );

	if( isNewMask==1 )	{
		ASSERT( g_I=GIP_NULL );							ASSERT( g_rad=GIP_NULL );		
		g_I = GIP_alloc( sizeof(float)*MN*2 );			g_rad=g_I+MN;
		mask=GIP_alloc( sizeof(int)*MN );
		GIP_canny_( hEdge,I_1,I_0,g_t0,g_t1,mask,g_I,g_rad,G_MASK_3 );
	}
//	GIP_save_IplImage_d( hEdge->M,hEdge->N,".\\trace\\roi_edge.jpg",hEdge->data,hEdge->N,0x0 );
	h_0=N*0.3;		h_1=N*0.7;		
	k_0=M*0.2;		k_1=M*0.9;
	a_0=200;		a_span=40;		a_step=5;	a_1=a_0+a_step*a_span;		
	if( 1 )	{
		h_0 = hIris->p_c-p_R/2;		h_1 = hIris->p_c+p_R/2;
	}
	if( region==IRIS_UPPER_EYELID )	{
		k_0=M*0.1;								k_1=MIN( M*0.9,hIris->p_r );
	}else if( region==IRIS_LOWER_EYELID )	{
		k_0=hIris->p_r;						k_1=M*0.9;
		a_0=-600;		a_span=40;			a_1=a_0+a_step*a_span;		
	}
	vote = GIP_calloc( sizeof(int),(h_1-h_0)*(k_1-k_0)*a_span );
	e_nz = 0;
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		if( mask[pos]!=1 )	continue;
		if( fabs((r-hIris->p_r)*(r-hIris->p_r)+(c-hIris->p_c)*(c-hIris->p_c)-p_R*p_R)<2.0 )
			continue;
	//	if( c<hIris->p_c-p_R*2 || c>hIris->p_c+p_R*2 )
	//		continue;
		e_data[pos]=255.0;
		e_nz++;
		for( k=k_0; k<k_1; k++ )	{
		for( a_no=0; a_no<a_span; a_no++ )	{
			a = a_0+a_no*a_step;
			if( a*(r-k)<0 )
				continue;
			w = sqrt( a*(r-k) );
			h = G_DOUBLE2INT( c+w );
			if( h<h_0 || h>=h_1 )		
				continue;
	//		xita = atan2( a,-2*(c-h) );
	//		if( fabs(xita-g_rad[GIP_RC2POS( r,c,N )]) > 0.30 ) 
	//			continue;
			pos = ((h-h_0)*(k_1-k_0)+k-k_0)*a_span+a_no;
			vote[pos]++;
			h = G_DOUBLE2INT( c-w );
			if( h<h_0 || h>=h_1 )		
				continue;
	//		xita = atan2( a,-2*(c-h) );
	//		if( fabs(xita-g_rad[GIP_RC2POS( r,c,N )]) > 0.30 ) 
	//			continue;
			pos = ((h-h_0)*(k_1-k_0)+k-k_0)*a_span+a_no;
			vote[pos]++;
		}
		}
		
OUTER_LOOP:	;	
	}
	}
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
	/*	e_nz=0;
		for( c=0; c<N; c++ )	{
			w=(c-h)*(c-h)*1.0/a+k;
			r = G_DOUBLE2INT( w );
			if( r<0 || r>=M )		continue;
			e_nz++;
			pos=GIP_RC2POS( r,c,N );
			e_data[pos]=0.0;
			while( r < M )	{
				pos=GIP_RC2POS( r,c,N );
				e_data[pos]=0.0;			data[pos]=0.0;
				r++;
			}
		}*/
//		hIris->lid_k=k;		
//		w = sqrt( a*(M-k) );
//		hIris->lid_h0=h-w;		hIris->lid_h1=h+w;		
//		hIris->nzNoise = e_nz;
	}
	
//	GIP_save_IplImage_d( hEdge->M,hEdge->N,".\\trace\\roi_edge_1.jpg",hEdge->data,hEdge->N,0x0 );
	if( isNewMask==1 )		{
		GIP_free( mask );			
		GIP_free( g_I );			
	}
	GIP_free( vote );
	GIP_OBJ_clear(hEdge);		GIP_free( hEdge );

}

/*
	检测iris的内外半径
	gradient operator + circular hash transform

	[REF]:	Optimizations in Iris Recognition
算法
	1 确定pupil的圆心,半径(基于pupil的圆度最高)
	2 确定sclera-boundary的圆心,半径
		利用与pupil圆心比较近的特点
	3 

参数：
	1 iris的内、外径取值范围 (0.05-0.5)*min(M,N)
	2 edge的梯度方向尚可，检测角可取为20度以内。
	3 edge的梯度值变化很大，取>=10
	4 图像强度变化较大，取I_0=0,I_1=200

ISSUE-NEEDE:
	1 现有的edge过粗，需要类canny edge-detector以减少edge数目。
	2 可尝试edge-preserve smooth filter，以增强边缘，减弱iris自身纹理
		"Recursive Implementation of Anisotropic Filtering"
	3 基于histogram确定更合理的I_0,I_1
	4 GIP_Circle_Verify的参数需要根据iris被遮挡的情况确定
	5 光照分析，以避免光照对边界梯度的影响

注意:
`	1 hObj不同于hIris->hObj
	2 edge direction与gradient vector垂直

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/17/2009	
*/	
int _IRIS_find_iris_2_( G_IRIS_ *hIris,int I_0,int I_1,double s_p,int method )	{
	GIa_HANDLE *hIa=hIris->hIa;
	GIP_OBJECT *hObj_0=hIris->hObj,*hEdge,*hObj,*hMask=hIris->hMask;
	int i,j,k,r,c,r1,c1,pos,v_pos,M=hObj_0->M,N=hObj_0->N,MN=M*N,ret=0,e_nz=0,dgr,dgr_0,dgr_1;
	int A_nz,*A_ptr,*A_off,*A_R,*temp=hIa->I_buffer,R_0,R_1,R,xita_deg,*mask,adj[8];
	GIP_FLOATING *data,*e_data,*g_I,*g_rad,gr,gc;
	float *vote,*v_w,v_max,v,g,g_t0,g_t1,rad,xita=PI/18,t,w;
	clock_t start=clock( );

	g_t1=_ToG_;				g_t0=g_t1/3;			
	xita_deg=G_RADIEN2DEGREE( xita )*2;
	R_0=(int)(MIN(M,N)*s_p);			R_1=R_0*10;
	_Arc_template_( M,N,R_0,R_1,&A_ptr,&A_off,&A_R,temp,0x0 );
//	R_0=MIN(M,N)*s_p;			R_1=R_0+10;

	start=clock( );
	hObj=GIP_OBJ_create( M,N,hObj_0->data,0x0 );
//	GIP_Smooth( hObj,GIP_SMOOTH_MEDIAN,hIris->hIa->D_buffer,0x0 );
	GIP_Smooth( hObj,GIP_SMOOTH_ANISO,hIris->hIa->D_buffer,0x0 );
	data = hObj->data;

	hEdge = GIP_OBJ_create( M,N,data,0x0 );
//	hEdge = GIP_OBJ_create( M,N,GIP_NULL,GIP_OBJ_2ZERO );
	e_data = hEdge->data;
	g_I = GIP_alloc( sizeof(GIP_FLOATING)*MN*2 );			g_rad=g_I+MN;
	mask = GIP_calloc( sizeof(int),MN );
	for( r = 2; r < M-2; r++ )	{		
	for( c = 2; c < N-2; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		if( data[pos]>I_1 || data[pos]<I_0 )		
		{	mask[pos]=-1;		continue;	}
		if( 1 )	{
			gr = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_R,0x0 );
			gc = -GIP_Mask_3_( pos,M,N,N,data,G_MASK_SOBEL_C,0x0 );
		}else	{
			gr = -GIP_Mask_5_( pos,M,N,N,data,G_MASK_FIVE_R,0x0 );
			gc = -GIP_Mask_5_( pos,M,N,N,data,G_MASK_FIVE_C,0x0 );
		}
		g_I[pos] = sqrt( gr*gr+gc*gc );
		g_rad[pos] = atan2( gr,gc );
	}
	}
	for( r = M*0.1; r < M*0.9; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		if( r==156 && c==107 )
			i=0;
		if( mask[pos]<0 )continue;
		if( g_I[pos] < g_t0 )	
			continue;
	//	w = G_RADIEN2DEGREE( fabs(g_rad[pos]) );
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
			mask[pos] = 1;					e_data[pos]=0.0;
		}else
			mask[pos] = 10;
	}
	}
	do{
		e_nz = 0;
		for( r = 1; r < M-1; r++ )	{		
		for( c = 1; c < N-1; c++ )	{
			pos=GIP_RC2POS( r,c,N );
			if( mask[pos]!=10 )		continue;
			mask[pos]=1;
			if( GIP_adjacent_2( M, N,N,mask,r,c,adj,GIP_ADJACENT_M )>0 )	{
//			if( mask[pos-1]==1 || mask[pos+1]==1 || mask[pos-N]==1 || mask[pos+N]==1 ){
				mask[pos]=1;		e_nz++;		e_data[pos]=0.0;
			}else
				mask[pos]=10;		
		}
		}
	}while( e_nz>0 );	

//	_eyelid_line_( M,N,mask,e_data,0x0 );
	GIP_save_IplImage_d( hEdge->M,hEdge->N,".\\trace\\CHT_edge.jpg",hEdge->data,hEdge->N,0x0 );
//voting			
	vote = GIP_calloc( sizeof(float),(MN+1)*(R_1-R_0) );		v_w=vote+MN*(R_1-R_0);
	for( i = 0; i < R_1-R_0; i++ )	{
		v_w[i] = atan2(1.0,(i+R_0) );
	}
	e_nz = 0;
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		if( mask[pos]!=1 )	continue;
		w = (1.0-(data[pos]-I_0)/(I_1-I_0))/R_0;			//偏向暗边界
		e_nz++;
		rad = g_rad[pos];
		dgr_0 = G_RADIEN2DEGREE( rad-xita );			
		for( j = dgr_0; j < dgr_0+xita_deg; j++ )	{
			dgr = j%360;
			for( k=A_ptr[dgr];	k<A_ptr[dgr+1];	k++ )	{
				v_pos = pos+A_off[k];
				if( v_pos<0 || v_pos>=MN )	continue;
//				e_data[v_pos]=0.0;
				v_pos = v_pos*(R_1-R_0)+A_R[k];
//				vote[v_pos]+=v_w[A_R[k]];		
				vote[v_pos]+=v_w[A_R[k]]+w;		
//				vote[v_pos]+=1.0+w;		
			}
		}
OUTER_LOOP:	;	//	break;
	}
	}

	v_max=0.0;
	for( i = MN*(R_1-R_0)-2; i > 0; i-- )	{
		v = vote[i];//+vote[i-1]+vote[i+1];
		if( v>v_max )	{
			v_max = v;		pos=i;
		}
	}
	R=pos%(R_1-R_0);			R+=R_0;		
	pos=pos/(R_1-R_0);
	r1=GIP_POS2R(pos,N);			c1=GIP_POS2C(pos,N);	
	if( GIP_Circle_Verify( hObj,r1,c1,R,xita,0.55,&(hIris->p_match),GIa_INNER_OBJ ) != 0 )	{	
		printf( "\t!!!Failed:\tpupile(%d,%d,%d),p_match=%4.2g\r\n",r1,c1,R,hIris->p_match );
		ret = GIRIST_ERROR_CIRCLE_VERIFY;		goto QUIT;	
	}
	hIris->p_r=r1;				hIris->p_c=c1;
	hIris->p_radius=R;
	_IRIS_find_eyelid_( hIris,M,N,data,g_I,g_rad,mask,IRIS_UPPER_EYELID,0x0 );
//	_IRIS_find_eyelid_( hIris,M,N,data,mask,IRIS_LOWER_EYELID,0x0 );

	v_max=0.0;
	k = MIN( R*1.5,R_1)-R_0;
	for( r = r1-R*0.1; r < r1+R*0.1; r++ )	{		
	for( c = c1-R*0.1; c < c1+R*0.1; c++ )	{
		i=GIP_RC2POS( r,c,N );
		for( j=k;	j < R_1-R_0; j++ )	{
			v_pos = i*(R_1-R_0)+j;
			v = vote[v_pos];//+vote[i-1]+vote[i+1];
			if( v>v_max )	{
				v_max = v;		pos=v_pos;
			}
		}
	}
	}
	R=pos%(R_1-R_0);			
	R+=R_0;		
	pos=pos/(R_1-R_0);
	r=GIP_POS2R(pos,N);			c=GIP_POS2C(pos,N);	
	hIris->s_r=r;				hIris->s_c=c;
	hIris->s_radius=R;
	t = (clock()-start)*1.0/CLOCKS_PER_SEC;
	GIP_MEMCOPY( e_data,hIris->hObj->data,sizeof(GIP_FLOATING)*MN );
	GIP_Parabolic( hEdge,hIris->ul_a,hIris->ul_h,hIris->ul_k,GIP_IMAG_ORIGIN );
	GIP_Circle( hEdge,hIris->p_r,hIris->p_c,hIris->p_radius,GIP_IMAG_ORIGIN );
	GIP_Circle( hEdge,hIris->s_r,hIris->s_c,hIris->s_radius,GIP_IMAG_ORIGIN );
	GIP_save_IplImage_d( hEdge->M,hEdge->N,".\\trace\\CHT.jpg",hEdge->data,hEdge->N,0x0 );
	if( GIP_Circle_Verify( hObj,r,c,R,xita,0.2,&(hIris->s_match),0x0 ) != 0 )	{	
		printf( "\t!!!Failed:\tsclera(%d,%d,%d),s_match=%4.2g\r\n",r,c,R,hIris->s_match );
		ret = GIRIST_ERROR_CIRCLE_VERIFY;		goto QUIT;	
	}
//	printf( "  {find_iris_2}: e_nz=%.3g,R_0=%g,R_1=%g,p_match=%.3g,s_match=%.3g,time=%g\r\n",e_nz*1.0/MN,hIris->p_radius,hIris->s_radius,
//		hIris->p_match,hIris->s_match,t );
QUIT:
	GIP_free( mask );
	GIP_free( g_I );
	GIP_OBJ_clear( hEdge );				GIP_free( hEdge );
	GIP_OBJ_clear( hObj );				GIP_free( hObj );
	GIP_free( vote );			GIP_free( A_ptr );

	return ret;
}

/*
	Circular Object Detecor on Global Symmetry

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/14/2006	
*/	
int _IRIS_find_pupil_3_( G_IRIS_ *hIris,int method )	{
	GIP_OBJECT *hObj=hIris->hObj;
	int i,r,c,pos,M=hObj->M,N=hObj->N,MN=M*N,ret=0;
	GIP_FLOATING *data = hObj->data;
	GIa_HANDLE *hIa=hIris->hIa;

	GIa_Circle_Detect_GPV( hObj,200,0x0,hIa->D_buffer,hIa->I_buffer,0x0 );

	return ret;
}

/*
	get the center coordinates of the pupil.

	ISSUE-NEEDED:
		1 需要更好的确定（r_0,c_0 )，可尝试梯度法

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/1/2008	
*/
void _IRIS_find_pupil_( G_IRIS_ *hIris,int method )	{
	GIP_OBJECT *hObj=hIris->hObj,pupil,*hPupil=&pupil,binary,*hBinary=&binary;
	int i,r,c,pos,M=hObj->M,N=hObj->N,M_p,N_p,n_1,n_2;
	double r_0,c_0,r_min,c_min,a,*data,u_1,u_2,th_1,thresh,err,R;
	char sTitle[80];
	GIa_REGION *hRegion;

	if( 1 )	{
		GIP_centre( hObj->M,hObj->N,hObj->N,hObj->data,hIris->hIa->I_buffer,&r_0,&c_0,0x0 );
		hIris->p_r = r_0;		hIris->p_c = c_0;

		M_p=MIN(150,M/2);		N_p=MIN(120,N/2);	
		GIP_OBJ_init( hPupil,M_p,N_p,0x0 );
		r=(int)(hIris->p_r-M_p/2),			c=(int)(hIris->p_c-N_p/2);
		GIP_OBJ_sub( hObj,hPupil,r,c,M_p,N_p,0x0 );		
		GIP_OBJ_trans( hPupil,GIP_TRANS_GRADIENT,0x0/*GIP_TRANS_ABS*/ );		
		GIP_GlobalThresh_2( hPupil,&thresh,0x0 );

		GIP_OBJ_init( hBinary,M_p,N_p,0x0 );
		GIP_Binarize( hPupil,hBinary,thresh,0x0 );
	//	sprintf( sTitle,".\\trace\\iris_bin_%d.bmp",M_p*N_p );	
	//	GIP_save_IplImage_d( hBinary->M,hBinary->N,sTitle,hBinary->data,hBinary->N,0x0 );
		err = _circle_fit( hBinary,&R,&r_0,&c_0,0x0 );
	}else	{
//		query_iris_pupil( hIris,&hRegion,0,0x0 );
//		err = GIa_circle_fit( hRegion,&R,&r_0,&c_0,0x0 );
	}
	r_0 = r_0+hIris->p_r-M_p/2;
	c_0 = c_0+hIris->p_c-N_p/2;
	hIris->p_r = r_0;		hIris->p_c = c_0; 
	hIris->p_radius = R;

	data = hObj->data;
	for( r = 0; r < hObj->M; r++ )	{		
	for( c = 0; c < hObj->N; c++ )	{
		a=sqrt( (r_0-r)*(r_0-r)+(c_0-c)*(c_0-c) );
		pos=GIP_RC2POS( r,c,hObj->N ); 
		if( fabs(a-R)<1.0 )	
		{	data[pos]=0.0;		}
	}
	}
	sprintf( sTitle,".\\trace\\iris_%d_1.jpg",M*N );	
	GIP_save_IplImage_d( hObj->M,hObj->N,sTitle,hObj->data,hObj->N,0x0 );
/*	if( err < err_min )	{	
	}*/

}

/*
	1 The remaping of raw co-ordinates (x,y) to the non concentric polar co-ordinate system
	2 (80,512)->(48,512)
	3 Image Enhancement

	注意：
		1 iris的内圆(pupil)总是在图形内部，但外圆有可能越界

	ISSUE-NEEDED
		1 the gray levels are adjusted by removing the peak illumination caused by light sources reflecting from the eye
		2 estimating and subtracting the slowly varying background illumination
		[REF] Image Normalization for Illumination Compensation in Facial Images
		3 标记被污染区域

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/7/2008	
	v0.2	cys
		1/19/2009
		增加enhancement
*/
void _IRIS_normalize_( G_IRIS_ *hIris,int method )	{
	GIP_OBJECT *hNormal=&(hIris->normal),*hObj=hIris->hObj;
	int i,r,c,M=hNormal->M,N=hNormal->N,r_o,c_o,pos_o,pos,M_0,I_0;
	GIP_FLOATING *data,*N_data=hNormal->data;
	double sR,sX,r_p,c_p,r_s,c_s,a,xita;
	char sTitle[80];
	GIa_HANDLE *hIa=hIris->hIa;

	M_0=80;
	I_0 = hIris->hIa->histo.I_0;
//	xita = atan2( hIris->p_r-hIris->s_r,hIris->p_c-hIris->s_c );
	data = GIP_alloc( sizeof(GIP_FLOATING)*M_0*N );
	for( r = 0; r < M_0; r++ )	{		
	for( c = 0; c < N; c++ )	{
		sX = 2*PI*c/N;
		sR = r*1.0/M_0;
//		r_p = hIris->p_r+hIris->p_radius*cos( sX );
//		c_p = hIris->p_c+hIris->p_radius*sin( sX );
		r_p = hIris->p_r+hIris->p_radius*sin( sX );
		c_p = hIris->p_c+hIris->p_radius*cos( sX );
		ASSERT( r_p>=0 && r_p<hObj->M && c_p>=0 && c_p<hObj->N );
		r_s = hIris->s_r+hIris->s_radius*sin( sX );
		c_s = hIris->s_c+hIris->s_radius*cos( sX );
//		r_s = hIris->s_r+hIris->s_radius*cos( sX );
//		c_s = hIris->s_c+hIris->s_radius*sin( sX );
		pos= GIP_RC2POS( r,c,N );	
		r_o = r_p+(r_s-r_p)*sR;
		c_o = c_p+(c_s-c_p)*sR;
		if( r_o<=0 || r_o>=hObj->M-1 || c_o<=0 || c_o>=hObj->N-1 )		//需作标记
			data[pos] = I_0;
		else	{
			if( (c_o-hIris->ul_h)*(c_o-hIris->ul_h)>hIris->ul_a*(r_o-hIris->ul_k) && r<M)
				hIris->nzNoise++;
	/*		a = 0.0;
			for( i = 0; i < 4; i++ )	{
				r_o = (int)(r_p+(r_s-r_p)*sR+r_off_[i]);
				c_o = (int)(c_p+(c_s-c_p)*sR+c_off_[i]);
				pos_o = GIP_RC2POS( r_o,c_o,hObj->N );	
				a +=  hObj->data[pos_o];
			}
			data[pos] = G_DOUBLE2INT( a/4.0 );*/
			a = GIP_sub_pixel( hObj->M,hObj->N,hObj->N,hObj->data,r_p+(r_s-r_p)*sR,c_p+(c_s-c_p)*sR,0x0 );
			data[pos] = G_DOUBLE2INT( a );
		}
		if( 1 )				
			if( r==M && c%3==0)		data[pos]=0.0;			
	}
	}
	if( 0 )	{			//	!!testing!!	
		sprintf( sTitle,".\\trace\\iris_normal_%d.bmp",g_nzSample );
		GIP_save_IplImage_d( M_0,N,sTitle,data,N,0x0 );
	}
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos= GIP_RC2POS( r,c,N );		
		N_data[pos] = data[pos];
	}
	}
	GIP_free( data );

//	GIP_save_IplImage_d( hNormal->M,hNormal->N,".\\trace\\iris_normal_0.jpg",hNormal->data,hNormal->N,0x0 );

	GIP_Histo_Equal( hNormal,256,0x0,hIa->I_buffer,0x0 );
//	GIP_save_IplImage_d( hNormal->M,hNormal->N,".\\trace\\iris_normal_1.jpg",hNormal->data,hNormal->N,0x0 );
//	GIP_Contrast_Enhance( hNormal,256,0,0 );

	if( 1 )	{			//	!!testing!!	
		sprintf( sTitle,".\\trace\\iris_roi_%d.bmp",g_nzSample );
		GIP_save_IplImage_d( hNormal->M,hNormal->N,sTitle,hNormal->data,hNormal->N,0x0 );
	}

	GIa_HANDLE_clear( hIa );
	GIa_HANDLE_init( hIa,hNormal,0x0 );
//	ASSERT( hIa->contrast==255 );

//	GIP_mean_deviation( hNormal->M,hNormal->N,hNormal->N,hNormal->data,&sR,&sX );

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/16/2008	
*/
void _IRIS_Feature_Vector( G_IRIS_ *hIris,G_FVector_ *hFV,int method )	{
	int V_len,g_M,g_N,nz=0,r,c,grid=8,pos,M,N;
	GIP_OBJECT *hROI,*hGabor_1,hGabor_2,mY_1,mY_2;
	double mean,devia;
	char *V;

/*	hROI=GIP_OBJ_create( M_normal,N_normal,GIP_NULL,0x0 );
	for( r = 0; r < M_normal; r++ )	{
	for( c = 0; c < N_normal; c++ )	{
			}
	}*/
	hROI=&(hIris->normal);			
	M=hROI->M*0.75;				N=hROI->N;
	ASSERT( M==M_normal && N==N_normal );
//	GIP_filter_ft( M,N,&mY_1,hROI,gabor_1.hMat,0x0 );	
	GIP_Convol_( &mY_1,hROI,gabor_1.hMat,0x0 );	
	GIP_Convol_( &mY_2,hROI,gabor_2.hMat,0x0 );	

	ASSERT( hFV->V_len==M_normal/grid*N_normal/grid*4 );
//	hIris->V = GIP_alloc( sizeof(double)*hIris->V_len );		
	V =	hFV->V;	
	for( r = 0; r < M; r+=grid )	{		
	for( c = 0; c < N; c+=grid )	{
		pos= GIP_RC2POS( r,c,hROI->N );
		GIP_mean_deviation( grid,grid,hROI->N,mY_1.data+pos,&mean,&devia );
		V[nz++]=mean;			V[nz++]=devia;
		GIP_mean_deviation( grid,grid,hROI->N,mY_2.data+pos,&mean,&devia );
		V[nz++]=mean;			V[nz++]=devia;
	}
	}
//	ASSERT( nz==hIris->V_len );

	GIP_OBJ_clear( &mY_1 );			
	GIP_OBJ_clear( &mY_2 );			
}

/*
	binary code from zero crossings of adjacent patch vectors

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/2/2009	
*/
void _binary_zero_crossings_( char *pvc,int v_len,double *pv,int pos_1,int pos_2,double *temp )	{
	int i,flag=0x01;
	ASSERT( sizeof(char)==1 );

	v_len=MIN( 9,v_len );
	for( i = 0; i < v_len; i++ )	{
		temp[i] = pv[pos_1+i]-pv[pos_2+i];
	}
	for( i = 0; i < v_len-1; i++ )	{
		if( /*temp[i]==0.0 ||*/ temp[i]*temp[i+1]<0.0 )
			BIT_SET( *pvc,flag );
		flag = flag<<1;
	}
}

static double _w_DCT_12[144];
/*
	only for _IRIS_Feature_DCT_ with 1_D vector len=2
	d_temp[N*2]

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/2/2009	
*/
void _DCT_12_( double *patch,double *d_temp )	{
	int k,n,N=12,isVerify=1;
	double s,*w=_w_DCT_12,err,*old;

	if( isVerify )		{		//!!testing!!
		old = d_temp+N;
		for( n = 0; n < N; n++ )	old[n]=patch[n];
	}
	for( k = 0; k < N; k++ )	{
		s = k==0.0 ? 2.0/N/sqrt(2.0) : 2.0/N;
		for( n = 0; n < N; n++ )	{
			w[k*N+n] = s*cos( PI*k*(2*n+1)/2/N );
//			printf( "%lf ",w[k*N+n] );
		}	
//		printf( "\r\n" );
	}
	for( k = 0; k < N; k++ )	{
		d_temp[k] = 0.0;
		for( n = 0; n < N; n++ )	d_temp[k] += patch[n]*w[k*N+n];
	}
	for( n = 0; n < N; n++ )
		patch[n] = d_temp[n];
	if( isVerify )		{		//!!testing!!
		err = 0.0;
		for( n = 0; n < N; n++ )	{
			d_temp[n] = 0.0;
			for( k = 0; k < N; k++ )	d_temp[n] += patch[k]*w[k*N+n];
			s = d_temp[n]*N/2-old[n];
			err += s*s;
		}
		err = sqrt( err );
		ASSERT( err < FLT_EPSILON );
	}
}

/*
	only for !!testing!!	

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/4/2009	
*/
void 	__testing_zero_crossing( int p_M,int p_N,int step,double *pv,int flag )	{
	char sTitle[80];
	int N=step,r,i,pos_1,pos_2;
	double *zc,*temp;

	sprintf( sTitle,".\\trace\\DCT_%d.dat",g_nzSample );
	D_ouput( p_M*p_N,pv,sTitle );

	zc = GIP_calloc( (p_M+1)*N,sizeof(double) );
	temp = zc+p_M*N;
	for( r = 0; r <p_M+1; r++ )	{
		pos_1=r*N;				
		pos_2=r*p_N;
		for( i = 0; i < N; i++ )	{
			zc[pos_1+i] = pv[pos_2+i]-pv[pos_2+p_N+i];
		}
	}
	D_ouput( p_M*N,zc,sTitle );

	GIP_MEMCLEAR( zc,p_M*N*sizeof(double) );
	for( r = 0; r <p_M-1; r++ )	{
		pos_1=r*p_N;				pos_2=(r+1)*p_N;
		for( i = 0; i < N; i++ )	{
			temp[i] = pv[pos_1+i]-pv[pos_2+i];
		}
		pos_1=r*N;
		for( i = 0; i < N-1; i++ )	{
			if( /*temp[i]==0.0 ||*/ temp[i]*temp[i+1]<0.0 )
				zc[pos_1+i]=1.0;
		}
	}
	sprintf( sTitle,".\\trace\\zero_crossing_%d.bmp",g_nzSample );
	GIP_save_IplImage_d( p_M,N,sTitle,zc,N,0x0 );

	GIP_free( zc );
}

/*
	With 781 such subfeatures, the final feature length was 2343 bits.
	[REF] Robust Iris Feature Extraction and Matching
	[REF] DCT-based Iris Recognition

	注意:
		hROI不一定等于hIris->normal,可能是变换后的结果

	ISSUE-NEEDED:
		1 需要和G_FVector_相匹配

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/2/2009	
*/
void _IRIS_Feature_DCT_( G_IRIS_ *hIris,GIP_OBJECT *hROI,G_FVector_ *hFV,int method )	{
	int i,j,p_len,no,r,c,pos_1,pos_2,M,N,p_M,p_N,step_M,step_N,ldu;
	char *codelets,*pvc,*V; 
	GIP_OBJECT *hBand;
	GIP_FLOATING *data,*B_data;
	double *patch,*pv,*d_temp,a,*dv;		//mean,devia,

	M=hROI->M;				N=hROI->N;
	ldu=N;					//ldu始终不变
	ASSERT( M==M_normal && N==N_normal );
	hBand=GIP_OBJ_create( (M/4-1)*8,N,GIP_NULL,0x0 );
	data=hROI->data,		B_data=hBand->data;
	for( i = 0; i < hBand->M; i++ )	{
		pos_1 = GIP_RC2POS( i,0,ldu );
		no = (i/8)*4+i%8;
		pos_2 = GIP_RC2POS( no,0,ldu );
		for( j = 0; j < ldu; j++ )	{
			B_data[pos_1+j] = data[pos_2+j];
		}
	}
//	GIP_save_IplImage_d( hBand->M,hBand->N,".\\trace\\iris_band_.bmp",hBand->data,hBand->N,0x0 );

	step_M=8;						step_N=6;
	p_len = step_N*2;			//double-overlap horizontally 
	p_M=hBand->M/step_M;			p_N=hBand->N/step_N-1;
	pv = GIP_alloc( sizeof(double)*((p_M*p_N+2)*p_len) );
	d_temp = pv+p_M*p_N*p_len;
//DCT transform
	for( r = 0; r < p_M; r++ )	{
	for( c = 0; c < p_N; c++ )	{
		patch = pv+(r*p_N+c)*p_len;
		pos_1 = GIP_RC2POS( r*step_M,c*step_N,ldu );
		for( i = 0; i < p_len; i++ )	{
			a = 0.0;
			for( j = 0; j < step_M; j++ )		a+=B_data[pos_1+i+ldu*j];
			patch[i] = a/p_M;
		}
		_DCT_12_( patch,d_temp );				//1D-DCT
	}
	}
	GIP_OBJ_clear( hBand );					GIP_free( hBand );
//	__testing_zero_crossing( p_M,p_N*p_len,p_len,pv,0x0 );
	switch( hFV->type )	{
	case G_FV_ROI:
		GIP_fourier_spectrum( hROI->M,hROI->N,hROI,hFV->V,0x0 );
//		GIP_MEMCOPY( hFV->V,hROI->data,sizeof(double)*hROI->M*hROI->N );
		break;
	case G_FV_COEFFICIENT:
		dv = (double*)(hFV->V);
		for( r = 0; r < p_M; r++ )	{
		for( c = 0; c < p_N; c++ )	{
			for( i = 0; i < hFV->bpSF; i++ )	{
				dv[(r*p_N+c)*hFV->bpSF+i] = pv[(r*p_N+c)*p_len+i];
			}
		}
		}
//		GIP_MEMCOPY( hFV->V,pv,sizeof(double)*p_M*p_N*p_len );
		break;
	default:
		codelets = GIP_alloc( sizeof(char)*p_M*p_N );
	//binary code from zero crossings of adjacent patch vectors
		for( i = 0; i < p_M; i++ )	{
		for( j = 0; j < p_N; j++ )	{
			pvc = codelets + i*p_N+j;
			*pvc=0.0;
			pos_1 = (i*p_N+j)*p_len;
			if( j+1<p_N )	{
				pos_2 = pos_1+p_len;
				_binary_zero_crossings_( pvc,p_len,pv,pos_1,pos_2,d_temp );
			}
			if( i+1<p_M )	{
				pos_2 = pos_1+p_N*p_len;
				_binary_zero_crossings_( pvc,p_len,pv,pos_1,pos_2,d_temp );
			}
		}
		}
		ASSERT( hFV->V_len==p_M*p_N );				//p_M*p_N*3/8
	//	hIris->V = GIP_alloc( sizeof(char)*hIris->V_len );		
		V =	hFV->V;	
		for( i = 0; i < p_M*p_N; i++ )	{
			V[i] = codelets[i];
		}
	//	GIP_save_IplImage_I8( p_M,p_N,".\\trace\\DCT_V_.bmp",V,p_N,0x0 );	
		GIP_free( codelets );
		break;
	}
	GIP_free( pv );
}

/*
	Gaussian filter to smooth the iris image
	
{
}
	*/

/*
[**]
	小马和豪哥，兰博和教练之间的关系：师徒 父子般的感情 及同志
[**]	

	A testing code

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/1/2008	
*/
G_IRIS_ *G_IRIS_process( GIP_IMAGE *hIMG,GIP_IMAGE *hIMG_0,int model,int flag,int *ret )	{
	double *temp,res,tol,dT=0.5,s,rmse,snr,v_1,v_2,*v;
	G_IRIS_ iris,*hIris;
	GIP_OBJECT *hObj,obj_0;
	GIP_IMAGE *hIMG_core;
	int M,N,isTrace=1,nKey,difus_step=1,i,n_limit=90000,I_0,I_1;
	char sTitle[80],sModel[80];

	*ret=GIP_OK;
//	GIP_init_heapinfo( );
	if( 0/*hIMG->height*hIMG->width>=n_limit*/ )	{
		s=sqrt(n_limit*1.0/hIMG->height/hIMG->width);
		M=hIMG->height*s;		N=hIMG->width*s;
		M &= -(1<<2);		N &= -(1<<2);
		hIMG_core = GIP_IMG_resize( hIMG_0,M,N,0x0 );
		GIP_import_IplImage( hIMG_core,&(obj_0.M),&(obj_0.N),&(obj_0.data),0x0 );
	
		cvReleaseImage( &hIMG_core );
		hIMG_core = GIP_IMG_resize( hIMG,M,N,0x0 );
		printf( "RESIZE:\tM:%d->%d,N:%d->%d\r\n",hIMG->height,M,hIMG->width,N );		
	}	else	{
		hIMG_core = hIMG;
		GIP_import_IplImage( hIMG_0,&(obj_0.M),&(obj_0.N),&(obj_0.data),0x0 );
	}
	hIris=GIP_alloc( sizeof(G_IRIS_) );
	G_IRIS_init_( hIMG_core,hIris );
	hObj=hIris->hObj;
	I_0=hIris->hIa->histo.I_0;			I_1=hIris->hIa->histo.I_1;
/*	if( _IRIS_find_pupil_(hIris,0x0)==0 && _IRIS_find_iris_(hIris,0)	{
		printf( "Find iris based on center\r\n" );
	}else */ 
//	if( _IRIS_find_iris_2_( hIris,0,200,0.05,0x0 )!=0 )	{			//based on circular hash transform
	if( _IRIS_find_iris_2_( hIris,I_0+(I_1-I_0)*_ToI_0,I_0+(I_1-I_0)*_ToI_1,0.05,0x0 )!=0 )	{			//based on circular hash transform
		*ret=GIRIST_ERROR_LOCATION;
		goto GIP_exit;		
	}
	_IRIS_normalize_( hIris,0x0 );

//	printf( "\r\nFINISH:\tRMSE=%g,SNR=%g,v_1=%g,v_2=%g\r\n",rmse,snr,v_1,v_2	);

GIP_exit:
	GIP_OBJ_clear( &obj_0 );

	return hIris;
}

/*
[**]
	快过年了,嘻嘻
[**]	

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/21/2009	
*/
void G_IRIS_save_training( GIP_LIB *hLib,G_FVector_ *fv,char *sLibPath,char *sample,int flag )	{
	int nEntry=hLib->nEntry,nSample=hLib->nSample,N;
	FILE *fp;

	fp=fopen( sLibPath,"wb" );
	ASSERT( fp!=0 );

	fwrite( fv,sizeof(G_FVector_),1,fp );
	fwrite( hLib,sizeof(GIP_LIB),1,fp );
/*	fwrite( hLib->W_ptr,sizeof(int),hLib->nClass+1,fp ); 
	fwrite( hLib->W_no,sizeof(int),hLib->nSample,fp ); 
*/	fwrite( hLib->arrEntry,sizeof(GIP_LIB_entry),nEntry,fp ); 

//	N=hLib->fv_n_0;
	N=fv->V_len;
	fwrite( sample,fv->unit_size,N*nSample,fp ); 

	fclose( fp );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/21/2009	
*/
void G_IRIS_load_training( GIP_LIB *hLib,G_FVector_ *fv,char *sLibPath,char **sample,int flag )	{
	int nEntry,nSample,N;
	FILE *fp;

	fp=fopen( sLibPath,"rb" );
	ASSERT( fp!=0 );

	fread( fv,sizeof(G_FVector_),1,fp );	
	fv->V= GIP_alloc( fv->unit_size*fv->V_len );
	fread( hLib,sizeof(GIP_LIB),1,fp );				
	hLib->W_fv=GIP_NULL;	
	hLib->W_ptr=GIP_NULL;				hLib->W_no=GIP_NULL;	
/*	hLib->W_ptr=GIP_alloc( sizeof(int)*(hLib->nClass+1) );		
	hLib->W_no=GIP_alloc( sizeof(int)*hLib->nSample );		
	fread( hLib->W_ptr,sizeof(int),hLib->nClass+1,fp ); 
	fread( hLib->W_no,sizeof(int),hLib->nSample,fp );*/ 
	nEntry=hLib->nEntry,	nSample=hLib->nSample,	N=hLib->fv_n_0;
	hLib->arrEntry=GIP_alloc( sizeof(GIP_LIB_entry)*nEntry );
	fread( hLib->arrEntry,sizeof(GIP_LIB_entry),nEntry,fp ); 

	N=fv->V_len;
	*sample=GIP_alloc( fv->unit_size*N*nSample );
	fread( *sample,fv->unit_size,N*nSample,fp ); 

	fclose( fp );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/24/2008	
*/
void _sort_lib_( GIP_LIB *hLib,double *s,int *temp,int flag )	{
	int i,j,nSample = hLib->nSample,*no,t;
	GIP_LIB_entry *hEntry;
	double a;
	FILE *fp;

	no=GIP_alloc( sizeof(int)*nSample );
	for( i=0; i < nSample; i++ )
		no[i] = i;
	for( i=0; i < nSample; i++ )	{
		for( j=i+1; j < nSample; j++ )	{
			if( s[i]<s[j] ){
				a=s[j];		s[j]=s[i];		s[i]=a;	
				t=no[j];	no[j]=no[i];	no[i]=t;	
			}
		}	
	}

	fp = fopen( ".\lib_sort.info","w" );
	for( i=0; i < nSample; i++ )	{
		t = no[i];
		hEntry = hLib->arrEntry+t;
		fprintf( fp,"%s\r\n",hEntry->sFilePath );
	}
	fclose( fp );

	GIP_free( no );
}

/*
[**]
	圣诞快乐,嘻嘻
[**]	

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/24/2008	
*/
void G_IRIS_training( char *sTrainPath,char *sLibPath,int flag )	{
	int nSample=0,N,m,i,nRet,iscolor=0x0,nOK=0,ret,cls_old=-1;		
	G_IRIS_ *hIris;
	GIP_LIB lib;
	G_FVector_ *hfv=&(lib.fv);
	GIP_IMAGE *src_image,*img;
	GIP_LIB_entry *hEntry;
	GIa_HANDLE* hIa;
	clock_t start=clock( ),start_0;
	double s=1.0/CLOCKS_PER_SEC,*w_sort,s_N=1.0/48/512;
	char *sample;

	GIP_init_heapinfo( );

//	GIP_GABOR_init( &gabor_1,M_normal,N_normal,3,1.5,0x0 );
//	GIP_GABOR_init( &gabor_2,M_normal,N_normal,4.5,1.5,0x0 );
	start_0=clock( );
	if( 0 )	{
	//	G_IRIS_load_training( &lib,&fv,"..\\iris_training_1.data",&sample,0x0 );
		GIP_LIB_load( &lib,sLibPath, 0x0 );
	}else	{
//		N=M_normal*N_normal/16;
		nRet = GIP_LIB_init( &lib,sTrainPath,0x0 );			ASSERT(nRet==0x0);
		G_FVector_init( hfv,M_normal,N_normal,0x0 );
//		lib.fv_n_0 = N;
		nSample=lib.nEntry;
		N = hfv->V_len;
	//	nSample=1;
		sample=(char*)GIP_alloc( hfv->unit_size*nSample*N );		
		w_sort=GIP_alloc( sizeof(double)*nSample );
		lib.W_fv=sample;
		printf( " no  R_0 R_1 p_mt s_mt \t mean devia noise\r\n" );
		for( i = 0; i < nSample; i++ )	{
			g_nzSample = nOK;
			hEntry = lib.arrEntry+i;
			if( hEntry->cls!=cls_old )	{
				printf( "%s\r\n",hEntry->sFilePath );
				cls_old = hEntry->cls;
			}
			start=clock( );
			src_image=cvLoadImage(hEntry->sFilePath,iscolor);		ASSERT( src_image!=0x0 );
			img = cvCloneImage( src_image );
//			hObj = GIP_IMG2OBJ( src_image,0x0 );
//			GIP_MASK_testing( hObj,0x0 );
//			GIP_OBJ_clear( hObj );		GIP_free( hObj );
			hIris = G_IRIS_process( img,src_image,0x0,flag,&ret );
			hIa = hIris->hIa;
			if( ret==GIP_OK )	{
				hEntry->status=1;
				//	_IRIS_Feature_Vector( hIris,&fv,0x0 );
				_IRIS_Feature_DCT_( hIris,&(hIris->normal),hfv,0x0 );
				memcpy( sample+nOK*N*hfv->unit_size,hfv->V,hfv->unit_size*N );
				w_sort[nOK] = hIris->s_match;
				nOK++;
			}
			printf( "%-4d %3g %3g %4.2g %4.2g\t %-5.4g %-5.3g %-5.2g\r\n",i+1,hIris->p_radius,hIris->s_radius,
				hIris->p_match,hIris->s_match,hIa->mean,hIa->deviation,hIris->nzNoise*s_N );	
			G_IRIS_clear( hIris );
			GIP_free( hIris );		hIris=GIP_NULL;
			cvReleaseImage( &src_image );
			cvReleaseImage( &img );
//			printf( "I_%d:\ttime=%g\r\n",i+1,(clock( )-start)*s );	
//			break;
		}
		lib.nSample = nOK;
		printf( "\r\nFINISH:fail=%g(%d),time=%g\r\n",(nSample-nOK)*1.0/nSample,nSample-nOK,(clock( )-start_0)*s );
	//	G_IRIS_save_training( &lib,&fv,"..\\iris_training_1.data",sample,0x0 );
	//	GIP_LIB_reducing( &lib,sample, 0x0 );
		_sort_lib_( &lib,w_sort,GIP_NULL,0x0 );
		GIP_free( w_sort );
		GIP_LIB_scan( &lib,0x0 );
		GIP_LIB_save( &lib,sLibPath, 0x0 );
	}
	if( 1 )
		GIP_LIB_discriminate( &lib,0x0 );

//	G_FVector_clear( &fv );
	GIP_LIB_clear( &lib );
//	GIP_free( sample );	
//	GIP_GABOR_clear( &gabor_1 );			GIP_GABOR_clear( &gabor_2 );

	GIP_exit_heapinfo( );

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/3/2009	
*/
double _euclid_distance_shift_( int v_len,double *V1,double *V2,int sft_0,int stf_1 )	{
	int i,sft;
	double dist=0.0,a;

	for( sft=sft_0; sft<=stf_1; sft++ )	{
		for( i = 0; i < v_len; i++ )	{
			a = V1[i]-V2[i];
			dist += a*a;
		}
	}
	
	dist = sqrt( dist/v_len );

	return dist;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/3/2009	
*/
void _IRIS_matching_( G_IRIS_ *hIris,G_FVector_ *hFv,GIP_LIB *hLib,char *sample,int flag )	{
	int i,nSample,no_min,V_len=hFv->V_len,bpSF=hFv->bpSF;
	int d_type=GIP_DISTANCE_EUCLID;
	char *V=hFv->V,*temp;
	double dist,dist_min=DBL_MAX,*dv_1,*dv_2;
	GIP_OBJECT *hROI,*hNormal;

	bpSF = 1;
	temp=GIP_alloc( sizeof(char)*V_len );
	nSample = hLib->nSample;

	hNormal = &(hIris->normal);
	hROI = GIP_OBJ_create( hNormal->M,hNormal->N,hNormal->data,0x0 );
	GIP_OBJ_shift( hNormal,hROI,0,-4,0x0 );
	_IRIS_Feature_DCT_( hIris,hROI,hFv,0x0 );

	d_type=hFv->type==G_FV_DCT_CODE ? GIP_DISTANCE_HAMMING : GIP_DISTANCE_EUCLID;
	for( i = 0; i < nSample; i++ )	{
		switch( d_type )	{
		case GIP_DISTANCE_HAMMING:
//			dist = _Hamming_distance_( V_len,V,sample+i*V_len,bpSF,hFv->nzSF,temp );
			break;
		default:
			dv_1=(double*)(V);			dv_2=(double*)(sample)+i*V_len;
			dist = _euclid_distance_shift_( V_len,dv_1,dv_2,-10,10 );
		}
		if( dist<dist_min )	{
			no_min=i;		dist_min = dist;
		}
	}

	GIP_OBJ_clear( hROI );			GIP_free( hROI );
	GIP_free( temp );
}

/*

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/3/2009	
*/
void G_IRIS_matching( char *sFilePath,char *sLibPath,int flag )	{
	int nSample=0,N,m,iscolor=0x0,ret;	//CV_LOAD_IMAGE_GRAYSCALE=0;;	
	G_IRIS_ *hIris;
	G_FVector_ fv;
	GIP_LIB lib;
	GIP_IMAGE *src_image,*img;
	clock_t start=clock( );
	double s=1.0/CLOCKS_PER_SEC;
	char *sample;

	GIP_init_heapinfo( );

	G_IRIS_load_training( &lib,&fv,"..\\iris_training_1.data",&sample,0x0 );

	src_image=cvLoadImage(sFilePath,iscolor);		ASSERT( src_image!=0x0 );
	img = cvCloneImage( src_image );
	hIris = G_IRIS_process( img,src_image,0x0,flag,&ret );
	if( ret==GIP_OK )	{
		//	_IRIS_Feature_Vector( hIris,&fv,0x0 );
		_IRIS_matching_( hIris,&fv,&lib,sample,0x0 );
	}

	G_IRIS_clear( hIris );
	GIP_free( hIris );		hIris=GIP_NULL;
	cvReleaseImage( &src_image );
	cvReleaseImage( &img );
	printf( "M_:\ttime=%g\r\n",(clock( )-start)*s );	

	GIP_LIB_clear( &lib );
	GIP_free( sample );	
	G_FVector_clear( &fv );			//G_FVector_clear( &fv_lib );

	GIP_exit_heapinfo( );

}
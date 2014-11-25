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
	针对不同的库暂采用不同的参数
	(I)	MMU iris database
		gs_Iris_G_0=8.0	
		g_II_s0=0.1,	g_II_s1=0.6
		G_MASK_3
	(II) CNKI iris database
		gs_Iris_G_0=5.0	
		g_II_s0=0.1,	g_II_s1=0.6
		G_MASK_5

static double _ToM_s=0.2,_ToM_p=0.6;		//thrshold of match ratio of sclera and pupil 
static double _ToG_= 5.0;					//thrshold of gradient for finding iris 
static double _ToI_0=0.1,_ToI_1=0.7;		//thrshold of intensity for finding iris 
*/

static GIP_GABOR gabor_1,gabor_2;
static int g_nzSample;
static int g_nF_1=0,g_nF_2=0,g_nF_3=0;

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/1/2008	
*/
void G_IRIS_clear( G_IRIS_ *hIris )	{
	if( hIris == GIP_NULL )	return;

	if( hIris->hObj != GIP_NULL )		{	
		GIP_OBJ_clear( hIris->hObj );	
		GIP_free( hIris->hObj );	hIris->hObj=GIP_NULL;
		GIP_OBJ_clear( &(hIris->normal) );	
	}
	if( hIris->g_I != GIP_NULL )		{	
		GIP_free( hIris->g_I );			hIris->g_I=GIP_NULL;
	}
	if( hIris->g_rad != GIP_NULL )		{	
		GIP_free( hIris->g_rad );		hIris->g_rad=GIP_NULL;
	}
	if( hIris->mask != GIP_NULL )		{	
		GIP_free( hIris->mask );		hIris->mask=GIP_NULL;
	}
/*	if( hIris->hMask != GIP_NULL )		{	
		GIP_OBJ_clear( hIris->hMask );	
		GIP_free( hIris->hMask );	hIris->hMask=GIP_NULL;
	}
*/
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
void G_IRIS_init_( double *control,GIP_OBJECT *hObj_0,G_IRIS_ *hIris )	{
	int M,N,MN,ret=GIP_OK,t_len;
	int M_normal=(int)(control[GIRIS_M_NORMAL]),N_normal=(int)(control[GIRIS_N_NORMAL]);
	GIP_FLOATING *val=GIP_NULL;
	double beta;
	GIP_OBJECT *hObj=GIP_NULL;
	GIa_HANDLE *hIa;

	memset( hIris,0x0,sizeof(G_IRIS_) );
//	hIris->hObj=GIP_alloc( sizeof(GIP_OBJECT) );			hObj=hIris->hObj;
//	GIP_import_IplImage( hIMG,&hObj->M,&hObj->N,&hObj->data,0x0 );
	hIris->hObj=GIP_OBJ_create( hObj_0->M,hObj_0->N,hObj_0->data,0x0 );
	hObj=hIris->hObj;
	M=hObj->M;			N=hObj->N;
	MN=M*N;
	t_len=MAX((M+2)*(N+2)*6,N);

//	hIris->hMask=GIP_OBJ_create( hObj->M,hObj->N,GIP_NULL,GIP_OBJ_2ZERO );
	hIris->mask=GIP_alloc( sizeof(int)*MN );
	hIris->g_I=GIP_alloc( sizeof(GIP_FLOATING)*MN );
	hIris->g_rad=GIP_alloc( sizeof(GIP_FLOATING)*MN );

	hIris->hIa=GIP_alloc( sizeof(GIa_HANDLE) );		//&(g_s_hIa);		
	hIa=hIris->hIa;
	GIa_HANDLE_init( hIa,hObj,M_normal*N_normal*2,0x0 );

//	GIP_Histo_Equal( hObj,256,0x0,GIP_NULL,0x0 );
//	GIP_Contrast_Enhance( hObj,256,0x0,GIP_NULL,0x0 );

//	hIris->alg = CGM_SS_DIRECT | CGM_A_AT | CGM_METHOD_MEDIAN;				//CGM_SS_ITER;	CGM_METHOD_DUAL
//	GIP_OBJ_init( &(hIris->normal),64,512,0x0 );
	GIP_OBJ_init( &(hIris->normal),M_normal,N_normal,0x0 );							//80->48

//	G_PRINTF( "  {_init_}: M=%d,N=%d,mean=%g,devia=%g,contrast=%g\r\n",hObj->M,hObj->N,
//		hIa->mean,hIa->deviation,hIa->contrast );

GIP_exit:
	;
}

/*
	for ._find_iris_2_

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
	for( i = 0; i < DGREE; i++ )			ASSERT( temp[i]==A_ptr[i+1] );

	GIP_OBJ_clear( hArc );				GIP_free( hArc );
}

/*
	传入圆心的范围及半径的范围，基于CHT及GIP_Circle_GMeasure确定圆的位置

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/31/2009	
*/
void  _IRIS_circle_refine_( double *control,G_IRIS_ *hIris,double R_0,double R_1,double r_0,double c_0,int grid,int type )	{
	GIP_OBJECT *hObj=hIris->hObj,*hEdge;
	int i,r,c,i_r,i_c,R,pos,M=hObj->M,N=hObj->N,flag,*mask=hIris->mask,*R_sum;
	int hR_1,hR_0,h_R,e_nz,nz,isDark=0,isNoX=0;
	double _ToM_,match_0,match_h,*w,match,a,a_2,h_r,h_c,dist;
	double g_t1=control[GIRIS_T_GRAD],xita=control[GIRIS_GRAD_BIAS];			
	GIP_FLOATING *g_I=hIris->g_I,*g_rad=hIris->g_rad,*data=hObj->data,*e_data;
	int I_0=hIris->hIa->histo.I_0,I_1=hIris->hIa->histo.I_1,*arrPos;

	if( type==IRIS_PUPIL_BOUNDARY )	{
		match_0=hIris->p_match;
		_ToM_ = control[GIRIS_T_MTCH_P];
		isDark=1;
		flag=GIa_INNER_OBJ;
	}else	{
		match_0=hIris->s_match;
		_ToM_ = control[GIRIS_T_MTCH_S];
		flag=0x0;
		isNoX = 1;
	/*	for( r = 1; r < M-1; r++ )	{	//锦上添花之技巧	
		for( c = 1; c < N-1; c++ )	{
			w = sqrt( (r-hIris->p_r)*(r-hIris->p_r)+(c-hIris->p_c)*(c-hIris->p_c) );
			if( w<hIris->p_radius*1.4 )	{
				pos=GIP_RC2POS( r,c,N );
				mask[pos]=-3;
			}
		}
		}*/
	}
	hR_0=INT_MAX;				hR_1=INT_MIN;			
	hEdge=GIP_OBJ_create( M,N,data,0x0 );
	e_data = hEdge->data;
	e_nz = 0;
	w = hIris->hIa->D_buffer;
	arrPos=hIris->hIa->I_buffer;		R_sum=arrPos+M*N;
	for( i = 0; i < MAX(M,N); i++ )		R_sum[i]=0;
	for( i = 0; i < M*N; i++ )	w[i]=0.0;
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos= GIP_RC2POS( r,c,N );
		if( pos==5126 )			//!!tesing!!
			i=0;
		if( mask[pos]!=1 )			continue;
		dist = sqrt( (r-r_0)*(r-r_0)+(c-c_0)*(c-c_0) );
		if( dist<R_0 || dist>R_1 )	
			continue;
		if( isNoX==1 && fabs(c-c_0)<dist*xita*2 )
			continue;
		a_2 = atan2( r_0-r,c_0-c );
		a = fabs(g_rad[pos]-a_2);				
		if( a>2*PI )a -= 2*PI;
		a = MIN( a,2*PI-a ); 
		if( a>xita*4 )	
			continue;
		hR_0 = MIN( hR_0,dist-grid );		hR_1 = MAX( hR_1,dist+grid );		
		for( i = dist-grid; i < dist+grid; i++ )
			R_sum[i]++;
		e_data[pos]=255.0;				e_nz++;
		if( isDark==1 )
			w[pos] = (1.0+(I_1-data[pos])*1.0/(I_1-I_0));
		else
			w[pos] = 1.0;
	}
	}
	ASSERT( e_nz!= 0 );
	GIP_save_IplImage_d( hEdge->M,hEdge->N,_T(".\\trace\\refine_edge.jpg"),hEdge->data,hEdge->N,0x0 );

	hR_0 = MAX( hR_0,R_0 );		hR_1 = MIN( hR_1,R_1 );
	match_h = 0.0;
	for( R=hR_0; R<=hR_1; R++ )	{
	//	if( R_sum[R]<e_nz/10 )
	//		continue;
		nz = GIP_circle_pos( hObj,r_0,c_0,R,arrPos,0 );
		for( i_r=-grid; i_r<=grid; i_r++ )	{
		for( i_c=-grid; i_c<=grid; i_c++ )	{
			match = 0.0;
			for( i = 0; i < nz; i++ )	{
				r = GIP_POS2R(arrPos[i],N);			c=GIP_POS2C(arrPos[i],N);
				r += i_r;							c+=i_c;
				if( !GIP_RC_VALID( r,0,M-1) || !GIP_RC_VALID( c,0,N-1) )
					continue;
				pos = GIP_RC2POS( r,c,N );			ASSERT( pos>=0 && pos<M*N );
				if( type==IRIS_PUPIL_BOUNDARY )
					match += w[pos]/R;
				else
					match += w[pos];
			}
			if( match_h<match )	{
				h_R=R;		h_r=r_0+i_r;		h_c=c_0+i_c;
				match_h = match;
			}			
		}
		}		
	}
	if( match_h > match_0 )	{
		if( type==IRIS_PUPIL_BOUNDARY )	{
			hIris->p_match=GIP_Circle_GMeasure(hObj,h_r,h_c,h_R,xita,g_I,g_rad,GIP_NULL,flag );		
			hIris->p_radius=h_R;		
			hIris->p_r=h_r;				hIris->p_c=h_c;	
		}else	{
			hIris->s_match=GIP_Circle_GMeasure(hObj,h_r,h_c,h_R,xita,g_I,g_rad,GIP_NULL,flag );		
			hIris->s_radius=h_R;			
			hIris->s_r=h_r;				hIris->s_c=h_c;	
		}
	}
	GIP_OBJ_clear( hEdge );			GIP_free( hEdge );
}

/*
	传入圆心的范围及半径的范围，基于CHT及GIP_Circle_GMeasure确定圆的位置

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2009	

void  _IRIS_circle_refine_( double *control,G_IRIS_ *hIris,double R_0,double R_1,double r_0,double c_0,double grid,int type )	{
	GIP_OBJECT *hObj=hIris->hObj,*hEdge;
	int i,r,c,R,pos,M=hObj->M,N=hObj->N,flag,hough_limit,*mask=hIris->mask;
	int hc_0,hc_1,hr_0,hr_1,hR_1,hR_0,x_R,x_c,x_r,h_r,h_c,h_R,e_nz,*map,isDark=0,isNoX=0;
	double _ToM_,match_0,match_h,radius,w,*hough,match,g_max,a,a_2;
	double g_t1=control[GIRIS_T_GRAD],xita=control[GIRIS_GRAD_BIAS];			
	GIP_FLOATING *g_I=hIris->g_I,*g_rad=hIris->g_rad,*data=hObj->data,*e_data;
	int I_0=hIris->hIa->histo.I_0,I_1=hIris->hIa->histo.I_1,*arrPos,s_loop;

	if( type==IRIS_PUPIL_BOUNDARY )	{
		match_0=hIris->p_match;
		_ToM_ = control[GIRIS_T_MTCH_P];
		isDark=1;
		flag=GIa_INNER_OBJ;
	}else	{
		match_0=hIris->s_match;
		_ToM_ = control[GIRIS_T_MTCH_S];
		flag=0x0;
		isNoX = 1;
	}
	hR_0=(int)(R_0);				hR_1=(int)(R_1+0.5);			
	hr_0=(int)(r_0-grid);			hr_1 =(int)(r_0+grid+0.5);
	hc_0=(int)(c_0-grid);			hc_1 =(int)(c_0+grid+0.5);
	x_R=hR_1-hR_0+1;				x_r=(hr_1-hr_0+1);			x_c=(hc_1-hc_0+1);
	hough_limit=x_R*x_r*x_c;
	hough=GIP_calloc( sizeof(double),hough_limit );

	g_max = 0.0;
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos= GIP_RC2POS( r,c,N );
		w = sqrt( (r-r_0)*(r-r_0)+(c-c_0)*(c-c_0) );
		if( w<hR_0 || w>hR_1 )	
			continue;
		g_max = MAX( g_max,g_I[pos] );
	}
	}
	hEdge=GIP_OBJ_create( M,N,data,0x0 );
	e_data = hEdge->data;
	e_nz = 0;
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos= GIP_RC2POS( r,c,N );
		if( pos==5126 )			//!!tesing!!
			i=0;
		w = sqrt( (r-r_0)*(r-r_0)+(c-c_0)*(c-c_0) );
		if( w<hR_0 || w>hR_1 )	
			continue;
		if( isNoX==1 && fabs(c-c_0)<w*xita*2 )
			continue;
//		if( g_I[pos]<g_max*0.5 )	
//			continue;
		if( mask[pos]!=1 )			continue;
		a_2 = atan2( r_0-r,c_0-c );
		a = fabs(g_rad[pos]-a_2);				
		if( a>2*PI )a -= 2*PI;
		a = MIN( a,2*PI-a ); 
		if( a>xita*4 )	
			continue;
		e_data[pos]=255.0;
		w = 1.0;			
		if( isDark==1 )				//偏向暗边界
			w += (I_1-data[pos])*1.0/(I_1-I_0);
//		+w = g_rad[pos];								//偏向强边界
//		w = 1.0;
		w *= fabs(cos(a_2 ));						//水平衰减
		e_nz++;
		for( h_r=hr_0; h_r<hr_1; h_r++ )	{
		for( h_c=hc_0; h_c<hc_1; h_c++ )	{
			h_R=G_DOUBLE2INT( sqrt( (h_r-r)*(h_r-r)+(h_c-c)*(h_c-c) ) );
			if( hR_0>h_R || h_R>hR_1 )	continue;
			pos=(h_R-hR_0)*x_r*x_c+(h_r-hr_0)*x_c+h_c-hc_0;
			hough[pos]+=w/h_R;			
		}
		}
	}
	}
	ASSERT( e_nz!= 0 );
	GIP_save_IplImage_d( hEdge->M,hEdge->N,_T(".\\trace\\refine_edge.jpg"),hEdge->data,hEdge->N,0x0 );
//	GIP_Circle( hEdge,h_r,h_c,hR_0,GIP_IMAG_ORIGIN );

	s_loop = hough_limit/10;
	arrPos = GIP_alloc( sizeof(int)*(s_loop+hR_1+1) );
	map=arrPos+s_loop;
	for( i = 0; i < hR_1+1; i++ )		map[i]=0;
	GIP_arr_top( hough,hough_limit,s_loop,arrPos,0x0 );
	match_h = 0.0;
	i = -1;
	while( ++i < s_loop )	{
		pos = arrPos[i];
		R=pos/(x_r*x_c);		r=(pos-R*x_r*x_c)/x_c;		c=pos-R*x_r*x_c-r*x_c;
		R+=hR_0;				r+=hr_0;					c+=hc_0;
	//	if( map[R]==1 )
	//		continue;
		map[R]=1;
		if( R==20 && type==IRIS_SCLERA_BOUNDARY )
			pos=0;
		match = GIP_Circle_GMeasure(hObj,r,c,R,xita,g_I,g_rad,mask,flag );
		if( match_h<match )	{
			match_h = match;
			h_R=R;		h_r=r;		h_c=c;
		}		
	}
	GIP_free( arrPos );

	if( match_h > match_0 )	{
		if( type==IRIS_PUPIL_BOUNDARY )	{
			hIris->p_radius=h_R;		hIris->p_match=match_h;		
			hIris->p_r=h_r;				hIris->p_c=h_c;	
		}else	{
			hIris->s_radius=h_R;		hIris->s_match=match_h;		
			hIris->s_r=h_r;				hIris->s_c=h_c;	
		}
	}
	GIP_OBJ_clear( hEdge );			GIP_free( hEdge );
}
*/

/*
	注意:
		1 z是相对值，即radius=z+R_0

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/19/2009	
*/
void  _IRIS_circle_search_( double *control,G_IRIS_ *hIris,GIP_OBJECT *hROI,PACT_order_list *l_vote,int R_1,int R_0,int z[4],int type )	{
	float match_0,match_8,match,radius,radius_8,_ToM_;
	int s_loop=0,*map,R,r,c,r_8,c_8,v_pos,pos,M=hROI->M,N=hROI->N,flag,z_2,z_3;
	int nz,nzCand,*ptr,*ind,*mask=hIris->mask;
	double xita=control[GIRIS_GRAD_BIAS],w;			//PI/18;

	nzCand = l_vote->nz;
	ptr=GIP_alloc( sizeof(int)*(nzCand+l_vote->compact+1) );		ind=ptr+l_vote->compact+1;
	PACT_ol_2ccs( l_vote,ptr,ind,0x0 );

	if( type==IRIS_PUPIL_BOUNDARY )	{
		match_0=match_8=hIris->p_match;
		radius=hIris->p_radius;		r=hIris->p_r;			c=hIris->p_c;	
		_ToM_ = control[GIRIS_T_MTCH_P];
		flag=GIa_INNER_OBJ;
//		z_2=R_0;		z_3=R_1;			
		z_2=z[0];		z_3=z[1];				
	}else	{
		for( r = 1; r < M-1; r++ )	{		
		for( c = 1; c < N-1; c++ )	{
			w = sqrt( (r-hIris->p_r)*(r-hIris->p_r)+(c-hIris->p_c)*(c-hIris->p_c) );
			if( w<hIris->p_radius*1.4 )	{
				pos=GIP_RC2POS( r,c,N );
				mask[pos]=-3;
			}
		}
		}
		match_0=match_8=hIris->s_match;
		radius=hIris->s_radius;		r=hIris->s_r;			c=hIris->s_c;
		z_2=z[2];		z_3=z[3];			
		_ToM_ = control[GIRIS_T_MTCH_S];
		flag=0x0;
	}
	map = GIP_calloc( sizeof(int),(R_1+1) );
	map[(int)(radius)] = 1;
	nz=nzCand-1;
	while( s_loop++<20 )	{
		while( nz>=0 )	{
			v_pos = ind[nz];
			if( (--nz)<0 )
				goto _LOOP_EXIT_;
			v_pos += l_vote->base;
			R=v_pos%(R_1-R_0)+R_0;
			if( (R-R_0>=z_2&&R-R_0<=z_3) && map[R]==0 )		//radius!=R)
				break;
		}
		radius=R;					map[R] = 1;
		pos=v_pos/(R_1-R_0);
		r=GIP_POS2R(pos,N);			c=GIP_POS2C(pos,N);	
//		ASSERT( r>=r1-grid && r<=r1+grid && c>=c1-grid && c<=c1+grid );
		match=GIP_Circle_GMeasure(hROI,r,c,radius,xita,hIris->g_I,hIris->g_rad,mask,flag );
		if( match > match_8 )	{
			match_8 = match;			radius_8=radius;				
			r_8=r;						c_8=c;	
		}
		if( match > _ToM_*2 )
		{	break;		}
	}
_LOOP_EXIT_:
	GIP_free( map );
	if( match_8 > match_0 )	{
		if( type==IRIS_PUPIL_BOUNDARY )	{
			hIris->p_radius=radius_8;		hIris->p_match=match_8;		
			hIris->p_r=r_8;					hIris->p_c=c_8;	
		}else	{
			hIris->s_radius=radius_8;		hIris->s_match=match_8;		
			hIris->s_r=r_8;					hIris->s_c=c_8;	
		}
	}
	GIP_free( ptr );
}

/*
	顺序计算histogram,smooth等变换，
	并计算对应的gradient
	注意：
		hObj的内容被改变

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/25/2009	
*/

void _IRIS_transform( double *control,G_IRIS_ *hIris,int flag )	{
	GIP_OBJECT *hObj=hIris->hObj;
	int i,M=hObj->M,N=hObj->N,MN=M*N,*mask;
	double q=control[GIRIS_SMOOTH_Q],_ToI_0=control[GIRIS_T_INTN_0],_ToI_1=control[GIRIS_T_INTN_1],g_t0,g_t1;
	int I_0,I_1,low,high;
	GIa_HANDLE *hIa = hIris->hIa;

	if( BIT_TEST( flag,GIP_M_HISTO ) )	{
		GIP_Histo_Equal( hObj,256,0x0,hIris->hIa->I_buffer,0x0 );
		GIa_HANDLE_clear( hIa );
		GIa_HANDLE_init( hIa,hObj,0,0x0 );
	}
	I_0=hIa->histo.I_0,						I_1=hIa->histo.I_1;
	low=(int)(I_0+(I_1-I_0)*_ToI_0);		high=(int)(I_0+(I_1-I_0)*_ToI_1);


	if( BIT_TEST( flag,GIP_M_SMOOTH ) )	{
	//	GIP_Smooth( hObj,GIP_SMOOTH_MEDIAN,hIris->hIa->D_buffer,q,0x0 );
		GIP_Smooth( hObj,GIP_SMOOTH_ANISO,hIris->hIa->D_buffer,q,0x0 );
	}

//	g_I = GIP_alloc( sizeof(GIP_FLOATING)*MN*2 );			g_rad=g_I+MN;
//	mask = GIP_calloc( sizeof(int),MN );
	mask = hIris->mask;
	for( i = 0; i < MN; i++ )
		mask[i] = 0;

	g_t1=(float)(control[GIRIS_T_GRAD]);		//_ToG_;				
//	g_t0=g_t1/3;			
	g_t0=g_t1/2;			
//	GIP_canny_( hObj,I_1,I_0,g_t0,g_t1,mask,g_I,g_rad,G_MASK_5 );
	GIP_lsu_( hObj,high,low,g_t0,g_t1,mask,hIris->g_I,hIris->g_rad,G_MASK_5 );
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
	6 存在先找到外边界的可能,需采用双峰检测以同时确定内 外径
	7 _IRIS_find_eyelid_效果不理想

注意:
`	1 hObj不同于hIris->hObj
	2 edge direction与gradient vector垂直
	3 返回_IRIS_on_vote_

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/17/2009	
*/	
int _IRIS_find_iris_2_( double *control,G_IRIS_ *hIris,int I_0,int I_1,int flag )	{
	GIa_HANDLE *hIa=hIris->hIa;
	GIP_OBJECT *hObj_0=hIris->hObj,*hEdge,*hObj;
	int i,j,k,r,c,r1,c1,pos,v_pos,M=hObj_0->M,N=hObj_0->N,MN=M*N,ret,e_nz=0,dgr,dgr_0;
	int *A_ptr,*A_off,*A_R,*temp=hIa->I_buffer,R_0,R_1,xita_deg,*mask=hIris->mask;
	int ldR,pos_1,pos_2,grid,z[6],v_pos_0,v_pos_1;
	GIP_FLOATING *data,*e_data,*g_I,*g_rad;
	float *vote,*v_w,v_max,v,rad,t,w,radius_1,radius_2,s_1=3.5,s_2=1.5,v_s=1000;
	clock_t start=clock( );
	double _ToM_s=control[GIRIS_T_MTCH_S],_ToM_p=control[GIRIS_T_MTCH_P],match,xita=control[GIRIS_GRAD_BIAS];//PI/18
	double s_p=control[GIRIS_T_R0],q=control[GIRIS_SMOOTH_Q];
	PACT_order_list l_vote;

	start=clock( );

	ret=GIP_OK;
/*	g_t1=(float)(control[GIRIS_T_GRAD]);		//_ToG_;				
//	g_t0=g_t1/3;			
	g_t0=g_t1/2;			
*/	xita_deg=G_RADIEN2DEGREE( xita )*2;
	R_0=(int)(MIN(M,N)*s_p);			R_1=R_0*10;
	_Arc_template_( M,N,R_0,R_1,&A_ptr,&A_off,&A_R,temp,0x0 );
	g_I=hIris->g_I,		g_rad=hIris->g_rad;
	GIP_MEMCLEAR( &l_vote,sizeof(PACT_order_list) );

	hObj=GIP_OBJ_create( M,N,hObj_0->data,0x0 );
	data = hObj->data;
	hEdge = GIP_OBJ_create( M,N,data,0x0 );
	e_data = hEdge->data;
//	_eyelid_line_( M,N,mask,e_data,0x0 );
//voting			
	vote = GIP_calloc( sizeof(float),(MN+1)*(R_1-R_0) );		v_w=vote+MN*(R_1-R_0);
	for( i = 0; i < R_1-R_0; i++ )	{
		v_w[i] = (float)( atan2(1.0,(i+R_0) ));
	}
	e_nz = 0;
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos=GIP_RC2POS( r,c,N );
//		if( hMask->data[pos]!=0 )	
//			continue;
		if( mask[pos]!=1 )			continue;
//		if( mask[pos]<0 )			continue;
		e_data[pos]=0.0;
//		w = (float)(1.0-(data[pos]-I_0)/(I_1-I_0))/R_0;			//偏向暗边界
		w = (float)(I_1-data[pos])/(I_1-I_0);
		ASSERT( w>= 0.0 );
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
				vote[v_pos]+=v_w[A_R[k]]+w/(A_R[k]+R_0);		
//				vote[v_pos]+=1.0+w;		
			}
		}
	}
	}
	GIP_save_IplImage_d( hEdge->M,hEdge->N,_T(".\\trace\\CHT_edge.jpg"),hEdge->data,hEdge->N,0x0 );
	if( e_nz==0 )	{
		G_PRINTF( _T("\t!!!e_nz=%d!!!\r\n"),e_nz );
		ret = GIRIST_ERROR_BOUNDARY;
		goto QUIT;
	}

//	ret = _IRIS_on_vote_( control,hIris,hObj,vote,R_1,R_0,0x0 );
	ldR=R_1-R_0;
	v_max=0.0;				pos_1=-1;
	for( r = 1; r < M-1; r++ )	{	
	for( c = 1; c < N-1; c++ )	{
		i=GIP_RC2POS( r,c,N );
		for( j=0;	j < R_1-R_0; j++ )	{
			v_pos = i*(R_1-R_0)+j;
	//		v = vote[v_pos];
			v = (vote[v_pos]+vote[v_pos-ldR]+vote[v_pos+ldR])/3;//+vote[v_pos-ldR*N]+vote[v_pos+ldR*N];
			if( v>v_max )	{
				v_max = v;		pos_1=v_pos;
			}
		}
	}
	}

	radius_1=(float)(pos_1%(R_1-R_0)+R_0);	
	if( radius_1<MIN(M,N)*0.1 )	{		//针对类似UBIRIS中的特小pupil
		s_1=4.5;		s_2=2.5;
	}
	pos=pos_1/(R_1-R_0);
	r1=GIP_POS2R(pos,N);						c1=GIP_POS2C(pos,N);	
	grid = MAX( (int)(radius_1*0.05),2 );
	z[0]=MAX(radius_1*0.8,radius_1-5)-R_0;		
	if( z[0]<0 )		z[0]=0;
	z[1]=MIN(radius_1*1.2,radius_1+5)-R_0;
	z[2]=(int)(radius_1*s_2-R_0);				z[3]=(int)(radius_1*s_1-R_0+1.0);
	z[4]=(int)(radius_1/4.5-R_0);				z[5]=(int)(radius_1/2.5-R_0+1.0);

	v_pos_0=GIP_RC2POS( r1-grid,c1-grid,N )*(R_1-R_0);
	v_pos_1=GIP_RC2POS( r1+grid,c1+grid+1,N )*(R_1-R_0);
	PACT_ol_init( &l_vote,(int)((v_max*v_s)*1.2),v_pos_1-v_pos_0,0x0 );
	l_vote.base = v_pos_0;
	v_max=0.0;				pos_2=-1;
	for( r = r1-grid; r <= r1+grid; r++ )	{	//尽量缩小范围，两者偏心并不厉害，也有助于normalize	
	for( c = c1-grid; c <= c1+grid; c++ )	{
		i=GIP_RC2POS( r,c,N );
		for( j=0;	j < R_1-R_0; j++ )	{
			if( !(j>=z[0]&&j<=z[1]) && !(j>=z[2]&&j<=z[3]) && !(j>=z[4]&&j<=z[5]))
				continue;
			v_pos = i*(R_1-R_0)+j;
//			v = vote[v_pos];
			v = (vote[v_pos]+vote[v_pos-ldR]+vote[v_pos+ldR])/3;//+vote[i-1]+vote[i+1];
			if( v==0.0 )
				continue;
			PACT_ol_insert( &l_vote,v_pos-v_pos_0,(int)(v*v_s+0.5) );
			if( v>v_max &&( (j>=z[2]&&j<=z[3])||(j>=z[4]&&j<=z[5])) )	{
				v_max = v;		pos_2=v_pos;
			}
		}
	}
	}
	if( pos_2==-1 )	{
		hIris->s_radius=0;			hIris->s_match=-1.0;
		hIris->s_r=0;				hIris->s_c=0;			
		G_PRINTF( _T("\t!!!Failed:\tNo candidate for sclera!!!\r\n"),hIris->s_r,hIris->s_c,hIris->s_radius,hIris->s_match );
		ret = GIRIST_ERROR_BOUNDARY;		goto QUIT;
	}
	radius_2=pos_2%(R_1-R_0)+R_0;	
	if( radius_1>radius_2 )	{		//radius_2有可能对应于pupil
		pos=pos_2/(R_1-R_0);
		r=GIP_POS2R(pos,N);			c=GIP_POS2C(pos,N);	
		if( (match=GIP_Circle_GMeasure(hObj,r,c,radius_2,xita,g_I,g_rad,mask,GIa_INNER_OBJ))>_ToM_p )	{		
			i=pos_1;		pos_1=pos_2;		pos_2=i;
			z[0]=MAX(radius_2*0.8,radius_2-5)-R_0;		
			if( z[0]<0 )		z[0]=0;
			z[1]=MIN(radius_2*1.2,radius_2+5)-R_0;
		}
	}
	hIris->p_radius=pos_1%(R_1-R_0)+R_0;		
	pos=pos_1/(R_1-R_0);
	hIris->p_r=GIP_POS2R(pos,N);			hIris->p_c=GIP_POS2C(pos,N);	
	hIris->p_match=GIP_Circle_GMeasure(hObj,hIris->p_r,hIris->p_c,hIris->p_radius,xita,g_I,g_rad,mask,GIa_INNER_OBJ );
	hIris->s_radius=pos_2%(R_1-R_0)+R_0;		
	pos=pos_2/(R_1-R_0);
	hIris->s_r=GIP_POS2R(pos,N);			hIris->s_c=GIP_POS2C(pos,N);
	if( hIris->s_radius<=hIris->p_radius*1.1 )
		hIris->s_match=-1.0;
	else
		hIris->s_match=GIP_Circle_GMeasure(hObj,hIris->s_r,hIris->s_c,hIris->s_radius,xita,g_I,g_rad,mask,0 );

	if( hIris->p_match<_ToM_p*1.5 )	{		//INNER_BOUNDARY re-search
		_IRIS_circle_search_( control,hIris,hObj,&l_vote,R_1,R_0,z,IRIS_PUPIL_BOUNDARY );
		if( hIris->p_match<_ToM_p )	{			
			G_PRINTF( _T("\t!!!Warning:\tpupile(%g,%g,%g),p_match=%4.2g\r\n"),hIris->p_r,hIris->p_c,hIris->p_radius,hIris->p_match );
			ret = GIRIST_WARN_INNER_BOUNDARY;			
		}
	}
	if( hIris->s_match<_ToM_s*2 || hIris->s_radius<hIris->p_radius )	{				//OUTER_BOUNDARY re-search
		_IRIS_circle_search_( control,hIris,hObj,&l_vote,R_1,R_0,z,IRIS_SCLERA_BOUNDARY );
		if( hIris->s_match<_ToM_s )	{
			G_PRINTF( _T("\t!!!Warning:\tsclera(%g,%g,%g),s_match=%4.2g\r\n"),hIris->s_r,hIris->s_c,hIris->s_radius,hIris->s_match );
			ret = GIRIST_WARN_OUTER_BOUNDARY;		
		}
	}
	w=hIris->p_radius;		
	if( hIris->s_radius<=hIris->p_radius )	{
		G_PRINTF( _T("\t!!!s_radius<=p_radius:\tp_radius=%4.2g,s_radius=%4.2g\r\n"),hIris->s_radius,hIris->p_radius );
		ret = GIRIST_ERROR_BOUNDARY;	
	}else if( hIris->p_r-w<0||hIris->p_r+w>=M || hIris->p_c-w<0||hIris->p_c+w>=N )	{
		G_PRINTF( _T("\t!!!Pupil is out of Image:\tp_radius=%4.2g,r=%4.2g,c=%4.2g\r\n"),hIris->p_radius,hIris->p_r,hIris->p_c );
		ret = GIRIST_ERROR_BOUNDARY;	
	}
	t = (clock()-start)*1.0/CLOCKS_PER_SEC;

QUIT:
	PACT_ol_clear( &l_vote );
//	GIP_free( mask );
//	GIP_free( g_I );
	GIP_OBJ_clear( hEdge );				GIP_free( hEdge );
	GIP_OBJ_clear( hObj );				GIP_free( hObj );
	GIP_free( vote );			GIP_free( A_ptr );

	return ret;
}

/*
	基于CHT的基准定位，是location模块的基础
	多数情况下比较准确，也有例外：
		1 UBIRIS/2/240	iris与pupil糊在一起		可尝试variation projection
		2 UBIRIS/2/114	瞳孔上有光影
		3 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/30/2009	
*/	
int _IRIS_find_iris_1_( double *control,G_IRIS_ *hIris,int flag )	{
	GIa_HANDLE *hIa=hIris->hIa;
	GIP_OBJECT *hObj_0=hIris->hObj,*hEdge,*hObj;
	int i,j,k,r,c,rX,cX,pos,v_pos,M=hObj_0->M,N=hObj_0->N,MN=M*N,ret,e_nz=0,dgr,dgr_0;
	int *A_ptr,*A_off,*A_R,*temp=hIa->I_buffer,R_0,R_1,xita_deg,*mask=hIris->mask,ldR,grid,pos_1;
	int I_0=hIris->hIa->histo.I_0,I_1=hIris->hIa->histo.I_1;
	GIP_FLOATING *data,*e_data,*g_I,*g_rad;
	float *vote,*v_w,v_max,v,t,w,radius_1,radius_2,s_1=3.5,s_2=1.5,v_s=1000;
	clock_t start=clock( );
	double _ToM_s=control[GIRIS_T_MTCH_S],_ToM_p=control[GIRIS_T_MTCH_P],xita=control[GIRIS_GRAD_BIAS];//PI/18
	double s_p=control[GIRIS_T_R0],q=control[GIRIS_SMOOTH_Q],rad,r_0,r_1;
	PACT_order_list l_vote;

	start=clock( );

	ret=GIP_OK;
	xita_deg=G_RADIEN2DEGREE( xita )*2;
	R_0=(int)(MIN(M,N)*s_p);			R_1=R_0*5;
	_Arc_template_( M,N,R_0,R_1,&A_ptr,&A_off,&A_R,temp,0x0 );
	g_I=hIris->g_I,		g_rad=hIris->g_rad;
	GIP_MEMCLEAR( &l_vote,sizeof(PACT_order_list) );

	hObj=GIP_OBJ_create( M,N,hObj_0->data,0x0 );
	data = hObj->data;
	hEdge = GIP_OBJ_create( M,N,data,0x0 );
	e_data = hEdge->data;
//voting			
	vote = GIP_calloc( sizeof(float),(MN+1)*(R_1-R_0) );		v_w=vote+MN*(R_1-R_0);
	for( i = 0; i < R_1-R_0; i++ )	{
		v_w[i] = (float)( atan2(1.0,(i+R_0) ));
	}
	e_nz = 0;
	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		if( mask[pos]!=1 )			continue;
//		if( mask[pos]<0 )			continue;
		e_data[pos]=0.0;
//		w = (float)(1.0-(data[pos]-I_0)/(I_1-I_0))/R_0;			//偏向暗边界
		w = (float)(I_1-data[pos])/(I_1-I_0);
		ASSERT( w>= 0.0 );
		e_nz++;
		rad = g_rad[pos];
		dgr_0 = G_RADIEN2DEGREE( rad-xita );			
		for( j = dgr_0; j < dgr_0+xita_deg; j++ )	{
			dgr = j%360;
			for( k=A_ptr[dgr];	k<A_ptr[dgr+1];	k++ )	{
				v_pos = pos+A_off[k];
				if( v_pos<0 || v_pos>=MN )	continue;
				v_pos = v_pos*(R_1-R_0)+A_R[k];
//				vote[v_pos]+=v_w[A_R[k]];		
				vote[v_pos]+=v_w[A_R[k]]+w/(A_R[k]+R_0);		
//				vote[v_pos]+=1.0+w;		
			}
		}
	}
	}
	GIP_save_IplImage_d( hEdge->M,hEdge->N,_T(".\\trace\\CHT_edge.jpg"),hEdge->data,hEdge->N,0x0 );
	if( e_nz==0 )	{
		G_PRINTF( _T("\t!!!e_nz=%d!!!\r\n"),e_nz );
		ret = GIRIST_ERROR_BOUNDARY;
		goto QUIT;
	}

	ldR=R_1-R_0;
	v_max=0.0;				pos_1=-1;
	if( 0 )	{
		for( r = 1; r < M-1; r++ )	{	
		for( c = 1; c < N-1; c++ )	{
			i=GIP_RC2POS( r,c,N );
			for( j=0;	j < R_1-R_0; j++ )	{
				v_pos = i*(R_1-R_0)+j;
		//		v = vote[v_pos];
				v = (vote[v_pos]+vote[v_pos-ldR]+vote[v_pos+ldR])/3;	//+vote[v_pos-ldR*N]+vote[v_pos+ldR*N])/5;
				if( v>v_max )	{
					v_max = v;		pos_1=v_pos;
				}
			}
		}
		}
		radius_1=(float)(pos_1%(R_1-R_0)+R_0);	
		pos=pos_1/(R_1-R_0);
		rX=GIP_POS2R(pos,N);						cX=GIP_POS2C(pos,N);	
		grid = MAX( (int)(radius_1*0.05),2 );
	}else	{
		for( r = 1; r < M-1; r++ )	{	
		for( c = 1; c < N-1; c++ )	{
			i=GIP_RC2POS( r,c,N );
			v = 0.0;
			for( j=0;	j < R_1-R_0; j++ )	{
				v_pos = i*(R_1-R_0)+j;
				v += vote[v_pos];
			}
			if( v>v_max )	{
				v_max=v;				rX=r;						cX=c;	
			}		
		}
		}		
		grid = 2;	
	}
	_IRIS_circle_refine_( control,hIris,R_0,R_1,rX,cX,grid,IRIS_PUPIL_BOUNDARY );
 	if( hIris->p_radius<MIN(M,N)*0.1 )	{		//针对类似UBIRIS中的特小pupil
		s_1=6.5;		s_2=1.5;
	}
	r_0=hIris->p_radius*s_2,		r_1=hIris->p_radius*s_1;
	_IRIS_circle_refine_( control,hIris,r_0,r_1,hIris->p_r,hIris->p_c,grid,IRIS_SCLERA_BOUNDARY );
	t = (clock()-start)*1.0/CLOCKS_PER_SEC;

QUIT:
	w=hIris->p_radius;		
	if( hIris->p_match<_ToM_p )	{			
//		G_PRINTF( _T("\t!!!Warning:\tpupile(%g,%g,%g),p_match=%4.2g\r\n"),hIris->p_r,hIris->p_c,hIris->p_radius,hIris->p_match );
		ret = GIRIST_WARN_INNER_BOUNDARY;			
	}
	if( hIris->s_match<_ToM_s )	{
//		G_PRINTF( _T("\t!!!Warning:\tsclera(%g,%g,%g),s_match=%4.2g\r\n"),hIris->s_r,hIris->s_c,hIris->s_radius,hIris->s_match );
		ret = GIRIST_WARN_OUTER_BOUNDARY;		
	}
	if( hIris->s_radius<=hIris->p_radius )	{
//		G_PRINTF( _T("\t!!!s_radius<=p_radius:\tp_radius=%4.2g,s_radius=%4.2g\r\n"),hIris->s_radius,hIris->p_radius );
		ret = GIRIST_ERROR_BOUNDARY;	
	}else if( hIris->p_r-w<0||hIris->p_r+w>=M || hIris->p_c-w<0||hIris->p_c+w>=N )	{
//		G_PRINTF( _T("\t!!!Pupil is out of Image:\tp_radius=%4.2g,r=%4.2g,c=%4.2g\r\n"),hIris->p_radius,hIris->p_r,hIris->p_c );
		ret = GIRIST_ERROR_BOUNDARY;	
	}
	PACT_ol_clear( &l_vote );
	GIP_OBJ_clear( hEdge );				GIP_free( hEdge );
	GIP_OBJ_clear( hObj );				GIP_free( hObj );
	GIP_free( vote );			GIP_free( A_ptr );

	return ret;
}

/*
	hIris已有初步结果，进一步细化
	注意：	
		ret_0==GIP_OK时，grid总体上偏小，基于原结果较准的假设

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/17/2009	
*/	
int _IRIS_find_iris_refine_( double *control,G_IRIS_ *hIris,int low,int high,int method,int ret_0 )	{
	GIP_OBJECT  *hObj=hIris->hObj;
	int M=hObj->M,N=hObj->N,MN=M*N,ret,i,pos,r,c,nz=0,*mask=hIris->mask,isIn,p_grid,s_grid,isReScan=0;
	float p_radius=hIris->p_radius,s_radius=hIris->s_radius,s_p,s_s,w;
//	double _ToM_s=control[GIRIS_T_MTCH_S],_ToM_p=control[GIRIS_T_MTCH_P];
	double 	xita=control[GIRIS_GRAD_BIAS],_ToM_s=control[GIRIS_T_MTCH_S],_ToM_p=control[GIRIS_T_MTCH_P];
	double R_0,R_1;
	GIP_FLOATING *g_I=hIris->g_I,*g_rad=hIris->g_rad;

	if( ret_0==GIRIST_ERROR_BOUNDARY)	{
		isReScan=1;
	}else if( hIris->p_match<_ToM_p*0.6 || hIris->s_match<_ToM_s/4 )	{
		isReScan=1;
	}else	{
		s_p=hIris->p_match;
		hIris->p_match = GIP_Circle_GMeasure(hObj,hIris->p_r,hIris->p_c,hIris->p_radius,xita,g_I,g_rad,GIP_NULL,GIa_INNER_OBJ );
		s_s=hIris->s_match;
		hIris->s_match=GIP_Circle_GMeasure(hObj,hIris->s_r,hIris->s_c,hIris->s_radius,xita,g_I,g_rad,GIP_NULL,0x0 );

		p_grid = hIris->p_match>_ToM_p ? 4 : hIris->p_match>_ToM_p*0.6 ? 8 : 10;
		s_grid = hIris->s_match>_ToM_s ? 5 : hIris->s_match>_ToM_s*0.6 ? 8 : 10;
		if( 1 )	{		//!!testing!!
			GIP_OBJECT *hO=GIP_OBJ_create( M,N,hObj->data,0x0 );;
			GIP_Circle( hO,hIris->p_r,hIris->p_c,p_radius-p_grid,GIP_IMAG_ORIGIN );
			GIP_Circle( hO,hIris->p_r,hIris->p_c,p_radius+p_grid,GIP_IMAG_ORIGIN );
			GIP_Circle( hO,hIris->s_r,hIris->s_c,s_radius-s_grid,GIP_IMAG_ORIGIN );
			GIP_Circle( hO,hIris->s_r,hIris->s_c,s_radius+s_grid,GIP_IMAG_ORIGIN );
//			GIP_Circle( hO,hIris->p_r,hIris->p_c,s_radius*0.8,GIP_IMAG_ORIGIN );
//			GIP_Circle( hO,hIris->p_r,hIris->p_c,s_radius*1.4,GIP_IMAG_ORIGIN );
			GIP_save_IplImage_d( M,N,GIP_path_(_T("refine_0"),g_nzSample),hO->data,N,0x0 );
			GIP_OBJ_clear( hO );		GIP_free( hO );
		}		
		

		if( s_p >_ToM_p )	{		//估计值比较可靠
			R_0=MAX( p_radius-p_grid,MIN(M,N)*control[GIRIS_T_R0] );
			R_1=MAX( p_radius+p_grid,R_0+p_grid/2 );
			_IRIS_circle_refine_( control,hIris,R_0,R_1,hIris->p_r,hIris->p_c,p_grid,IRIS_PUPIL_BOUNDARY );
		}
		if( s_s >_ToM_s )	{		//估计值比较可靠
			_IRIS_circle_refine_( control,hIris,s_radius-s_grid,s_radius+s_grid,hIris->s_r,hIris->s_c,s_grid,IRIS_SCLERA_BOUNDARY );
		}		
		if( hIris->p_match<_ToM_p || hIris->s_match<_ToM_s )	{
			isReScan=0;
		}
	}

	if( isReScan==1 )	{
	//	ret = _IRIS_find_iris_2_( control,hIris,low,high,GIP_M_ALL );
		ret = _IRIS_find_iris_1_( control,hIris,GIP_M_ALL );
		g_nF_3++;
	}	else
		ret = GIP_OK;

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
	利用iris的内外边界已确定

	v0.1	cys
		2/18/2009
*/
void _IRIS_normalize_0_( G_IRIS_ *hIris,int method )	{
	GIP_OBJECT *hObj=hIris->hObj,*hRing;
	int i,r,c,M=hObj->M,N=hObj->N,MN=M*N,pos,I_0,I_1,nz=0,*mask;
	GIP_FLOATING *data=GIP_NULL;
	double w,R_0=hIris->p_radius,R_1=hIris->s_radius,*g_I,*g_rad;
	GIa_HANDLE *hIa=hIris->hIa;

	hRing = GIP_OBJ_create( M,N,hObj->data,0x0 );
	data=hRing->data;

	R_1 = R_0+(R_1-R_0)*0.7;
	I_0 = 256;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		w = sqrt( (r-hIris->s_r)*(r-hIris->s_r)+(c-hIris->s_c)*(c-hIris->s_c) );
		if( w>R_1 )
			continue;
		w = sqrt( (r-hIris->p_r)*(r-hIris->p_r)+(c-hIris->p_c)*(c-hIris->p_c) );
		if( w<hIris->p_radius )
			continue;
		pos= GIP_RC2POS( r,c,N );	
		I_0 = MIN( I_0,data[pos] );
	}
	}	
	
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos= GIP_RC2POS( r,c,N );	
		w = sqrt( (r-hIris->s_r)*(r-hIris->s_r)+(c-hIris->s_c)*(c-hIris->s_c) );
		if( w>R_1 )	{
			data[pos] = I_0;		continue;
		}
		w = sqrt( (r-hIris->p_r)*(r-hIris->p_r)+(c-hIris->p_c)*(c-hIris->p_c) );
		if( w<hIris->p_radius )	{
			data[pos] = I_0;		continue;
		}
		nz++;
	}
	}
	GIP_Histo_Equal( hRing,256,0x0,hIa->I_buffer,0x0 );

	g_I=hIa->D_buffer;			g_rad=g_I+MN;
	mask=hIa->I_buffer;
	I_0=1;			I_1=50;
	nz = GIP_canny_( hRing,I_1,I_0,5,20,mask,g_I,g_rad,0x0 );
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos= GIP_RC2POS( r,c,N );	
		if( mask[pos]==1 )		data[pos]=255.0;
	}
	}


//	hough_line_detect( );
//	_IRIS_find_eyelid_( hRing,M,N,hNormal->data,GIP_NULL,GIP_NULL,GIP_NULL,IRIS_NOISE_EYELID,0x0 );
	
	GIP_save_IplImage_d( hRing->M,hRing->N,GIP_path_(_T("iris_ring"),g_nzSample),hRing->data,hRing->N,0x0 );
	GIP_OBJ_clear( hRing );			GIP_free( hRing );
}

double _IRIS_angle_( G_IRIS_ *hIris,int flag )	{
	GIP_OBJECT *hObj=hIris->hObj;
	int r,c,M=hObj->M,N=hObj->N,pos,w_nz=0;
	GIP_FLOATING *data=hObj->data;
	double p_r,p_c,s_r,s_c,w_r,w_c,R_1,R_0,d_p,d_s,ang;
	
	p_r=hIris->p_r,		p_c=hIris->p_c;
	s_r=hIris->s_r,		s_c=hIris->s_c;
	R_0 = hIris->p_radius;
	R_1 = R_0+(hIris->s_radius-hIris->p_radius)*0.6;

	w_r=0.0,			w_c=0.0;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		if( r<p_r )
			continue;
		pos= GIP_RC2POS( r,c,N );	
		d_p = sqrt( (p_r-r)*(p_r-r)+(p_c-c)*(p_c-c) );
		d_s = sqrt( (s_r-r)*(s_r-r)+(s_c-c)*(s_c-c) );
		if( d_p<=R_0 || d_s>= R_1 )
			continue;
		w_nz++;
		w_r += (r-p_r)*data[pos]*data[pos];
		w_c += (c-p_c)*data[pos]*data[pos];		
	}
	}

	w_r /= w_nz;			w_c /= w_nz;
	ang = atan2( w_r,w_c );
	ang = G_RADIEN2DEGREE( ang );

	return ang;

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
	int i,r,c,M=hNormal->M,N=hNormal->N,pos_o,pos,M_0,I_0;
	GIP_FLOATING *data,*N_data=hNormal->data,r_o,c_o;
	double sR,sX,r_p,c_p,r_s,c_s,a,xita,d_sum;
	GIa_HANDLE *hIa=hIris->hIa;

//	_IRIS_normalize_0_( hIris,0x0 );
//	 _IRIS_angle_( hIris,0x0 );

	M_0=80;
	I_0 = hIris->hIa->histo.I_0;
//	xita = atan2( hIris->p_r-hIris->s_r,hIris->p_c-hIris->s_c );
	d_sum = 0.0;
	data = GIP_alloc( sizeof(GIP_FLOATING)*M_0*N );
	for( c = 0; c < N; c++ )	{
		sX = 2*PI*c/N;
		r_p = hIris->p_r+hIris->p_radius*sin( sX );
		c_p = hIris->p_c+hIris->p_radius*cos( sX );
		ASSERT( r_p>=0 && r_p<hObj->M && c_p>=0 && c_p<hObj->N );
		r_s = hIris->s_r+hIris->s_radius*sin( sX );
		c_s = hIris->s_c+hIris->s_radius*cos( sX );
		d_sum += sqrt( (r_p-r_s)*(r_p-r_s)+(c_p-c_s)*(c_p-c_s) );
		for( r = 0; r < M_0; r++ )	{		
			sR = r*1.0/M_0;
			pos= GIP_RC2POS( r,c,N );	
			r_o = r_p+(r_s-r_p)*sR;
			c_o = c_p+(c_s-c_p)*sR;
			if( r_o<=0 || r_o>=hObj->M-1 || c_o<=0 || c_o>=hObj->N-1 )		//需作标记
				data[pos] = I_0;
			else	{
		//		if( (c_o-hIris->ul_h)*(c_o-hIris->ul_h)>hIris->ul_a*(r_o-hIris->ul_k) && r<M)
		//			hIris->nzNoise++;
		/*		a = 0.0;
				for( i = 0; i < 4; i++ )	{
					r_o = (int)(r_p+(r_s-r_p)*sR+r_off_[i]);
					c_o = (int)(c_p+(c_s-c_p)*sR+c_off_[i]);
					pos_o = GIP_RC2POS( r_o,c_o,hObj->N );	
					a +=  hObj->data[pos_o];
				}
				data[pos] = G_DOUBLE2INT( a/4.0 );*/
				a = GIP_sub_pixel( hObj->M,hObj->N,hObj->N,hObj->data,r_o,c_o,0x0 );
				data[pos] = G_DOUBLE2INT( a );
				ASSERT( data[pos]>=0 && data[pos]<=255 );
	//			data[pos] = hObj->data[GIP_RC2POS(r_o,c_o,hObj->N)];
			}
			if( 1 )				
				if( r==M && c%3==0)		data[pos]=0.0;			
		}
	}
//	hIris->sample = M*N/nzSample;			//采样率
	hIris->ecc = d_sum/N/(hIris->s_radius-hIris->p_radius);
	if( 0 )	{			//	!!testing!!	
		GIP_save_IplImage_d( M_0,N,GIP_path_(_T("iris_normal"),g_nzSample),data,N,0x0 );
	}
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos= GIP_RC2POS( r,c,N );		
		N_data[pos] = data[pos];
	}
	}
	GIP_free( data );

//	GIP_save_IplImage_d( hNormal->M,hNormal->N,".\\trace\\iris_normal_0.jpg",hNormal->data,hNormal->N,0x0 );

//	GIP_Smooth( hNormal,GIP_SMOOTH_MEDIAN,hIa->D_buffer,0.0,0x0 );
//	GIP_Smooth( hNormal,GIP_SMOOTH_ANISO,hIa->D_buffer,0x0 );
	GIP_Histo_Equal( hNormal,256,0x0,hIa->I_buffer,0x0 );
//	GIP_save_IplImage_d( hNormal->M,hNormal->N,".\\trace\\iris_normal_1.jpg",hNormal->data,hNormal->N,0x0 );
//	GIP_Contrast_Enhance( hNormal,256,0,0 );
//	_IRIS_find_eyelid_( hIris,M,N,hNormal->data,GIP_NULL,GIP_NULL,GIP_NULL,IRIS_NOISE_EYELID,0x0 );

	if( 1 )	{			//	!!testing!!	
		GIP_save_IplImage_d( hNormal->M,hNormal->N,GIP_path_(_T("iris_roi"),g_nzSample),hNormal->data,hNormal->N,0x0 );
	}

	GIa_HANDLE_clear( hIa );
	GIa_HANDLE_init( hIa,hNormal,0,0x0 );
//	ASSERT( hIa->contrast==255 );

//	GIP_mean_deviation( hNormal->M,hNormal->N,hNormal->N,hNormal->data,&sR,&sX );
}

/*
	[REF]: A Study of Segmentation and Normalization for Iris Recognition Systems 
	无效：
		1 实测偏心率非常小，接近于同心圆。参见2_22_ecc.info
		2 xita(按 Global Minimum Contour-Based Method得到的起始角)总是为0
		3 应适合于active contour边界

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/22/2009	
*/
void _IRIS_contour_normalize_( G_IRIS_ *hIris,int method )	{
	GIP_OBJECT *hNormal=&(hIris->normal),*hObj=hIris->hObj;
	int i,no,no_min,r,c,M=hNormal->M,N=hNormal->N,pos_o,pos,M_0,I_0;
	GIP_FLOATING *data,*N_data=hNormal->data,r_o,c_o;
	double sR,sX,r_p,c_p,r_s,c_s,a,xita,dr,dc,d_min,d,offcenter;
	GIa_HANDLE *hIa=hIris->hIa;

	M_0=80;
	I_0 = hIris->hIa->histo.I_0;
/*	r_p=GIP_alloc( sizeof(double)*N*4 );		c_p=r_p+N;
	r_s=c_p+N;									c_s=r_s+N;
	for( i = 0; i < N; i++ )	{
		sX = 2*PI*i/N;
		r_p[i] = hIris->p_r+hIris->p_radius*sin( sX );
		c_p[i] = hIris->p_c+hIris->p_radius*cos( sX );
		r_s[i] = hIris->s_r+hIris->s_radius*sin( sX );
		c_s[i] = hIris->s_c+hIris->s_radius*cos( sX );
		ASSERT( r_p[i]>=0 && r_p[i]<hObj->M && c_p[i]>=0 && c_p[i]<hObj->N );
	}*/
	no_min = -1;		d_min=DBL_MAX;
	for( pos=-10; pos < 10; pos++ )	{
		xita = pos*0.1/10;
		d = 0.0;
		for( i = 0; i < N; i++ )	{
			sX = 2*PI*i/N;
			r_p = hIris->p_r+hIris->p_radius*sin( sX+xita );
			c_p = hIris->p_c+hIris->p_radius*cos( sX+xita );
			r_s = hIris->s_r+hIris->s_radius*sin( sX );
			c_s = hIris->s_c+hIris->s_radius*cos( sX );
			dr = r_p-r_s;			dc=c_p-c_s;
			d += sqrt( dr*dr+dc*dc );
		}
		if( d < d_min )	{
			d_min = d;		no_min=pos;
		}
	}
	ASSERT( no_min==0 );
	hIris->ecc = d_min/N/(hIris->s_radius-hIris->p_radius);
//	hIris->x = (hIris->ecc-1.0)*100.0;

	data = GIP_alloc( sizeof(GIP_FLOATING)*M_0*N );
	xita = 0.0;
	for( c = 0; c < N; c++ )	{
		sX = 2*PI*c/N;
		for( r = 0; r < M_0; r++ )	{		
			sR = r*1.0/M_0;
			pos= GIP_RC2POS( r,c,N );	
			r_p = hIris->p_r+hIris->p_radius*sin( sX+xita );
			c_p = hIris->p_c+hIris->p_radius*cos( sX+xita );
			r_s = hIris->s_r+hIris->s_radius*sin( sX );
			c_s = hIris->s_c+hIris->s_radius*cos( sX );
			r_o = r_p+(r_s-r_p)*sR;
			c_o = c_p+(c_s-c_p)*sR;
			if( r_o<=0 || r_o>=hObj->M-1 || c_o<=0 || c_o>=hObj->N-1 )		//需作标记
				data[pos] = I_0;
			else	{
				a = GIP_sub_pixel( hObj->M,hObj->N,hObj->N,hObj->data,r_o,c_o,0x0 );
				data[pos] = G_DOUBLE2INT( a );
				ASSERT( data[pos]>=0 && data[pos]<=255 );
			}
			if( 1 )				
				if( r==M && c%3==0)		data[pos]=0.0;		
		}
	}
//	GIP_free( r_p );
	if( 0 )	{			//	!!testing!!	
		GIP_save_IplImage_d( M_0,N,GIP_path_(_T("iris_normal"),g_nzSample),data,N,0x0 );
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
//	_IRIS_find_eyelid_( hIris,M,N,hNormal->data,GIP_NULL,GIP_NULL,GIP_NULL,IRIS_NOISE_EYELID,0x0 );

	if( 1 )	{			//	!!testing!!	
		GIP_save_IplImage_d( hNormal->M,hNormal->N,GIP_path_(_T("iris_roi"),g_nzSample),hNormal->data,hNormal->N,0x0 );
	}

	GIa_HANDLE_clear( hIa );
	GIa_HANDLE_init( hIa,hNormal,0,0x0 );
//	ASSERT( hIa->contrast==255 );

//	GIP_mean_deviation( hNormal->M,hNormal->N,hNormal->N,hNormal->data,&sR,&sX );

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/16/2008	
*/
void _IRIS_Feature_MALI_( G_IRIS_ *hIris,G_FVector_ *hFV,int method )	{
	int V_len,g_M,g_N,nz=0,r,c,grid=8,pos,M,N;
	GIP_OBJECT *hROI,mY_1,mY_2;
	double mean,devia,*dv;

	hROI=&(hIris->normal);			
	M=hROI->M*0.75;				N=hROI->N;
//	ASSERT( M==M_normal && N==N_normal );
//	GIP_filter_ft( M,N,&mY_1,hROI,gabor_1.hMat,0x0 );	
	GIP_Convol_( &mY_1,hROI,gabor_1.hMat,0x0 );	
	GIP_Convol_( &mY_2,hROI,gabor_2.hMat,0x0 );	

//	ASSERT( hFV->V_len==M_normal/grid*N_normal/grid*4 );
//	hIris->V = GIP_alloc( sizeof(double)*hIris->V_len );		
	dv = (double*)hFV->V;	
	for( r = 0; r < M; r+=grid )	{		
	for( c = 0; c < N; c+=grid )	{
		pos= GIP_RC2POS( r,c,hROI->N );
		GIP_mean_deviation( grid,grid,hROI->N,mY_1.data+pos,&mean,&devia );
		dv[nz++]=mean;			dv[nz++]=devia;
		GIP_mean_deviation( grid,grid,hROI->N,mY_2.data+pos,&mean,&devia );
		dv[nz++]=mean;			dv[nz++]=devia;
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
void _binary_zero_crossings_( unsigned char *pvc,int v_len,double *pv,int pos_1,int pos_2,double *temp )	{
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
//			G_PRINTF( "%lf ",w[k*N+n] );
		}	
//		G_PRINTF( "\r\n" );
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
	int N=step,r,i,pos_1,pos_2;
	double *zc,*temp;

	D_ouput( p_M*p_N,pv,GIP_path_(_T("DCT"),g_nzSample) );

	zc = GIP_calloc( (p_M+1)*N,sizeof(double) );
	temp = zc+p_M*N;
	for( r = 0; r <p_M+1; r++ )	{
		pos_1=r*N;				
		pos_2=r*p_N;
		for( i = 0; i < N; i++ )	{
			zc[pos_1+i] = pv[pos_2+i]-pv[pos_2+p_N+i];
		}
	}
	D_ouput( p_M*N,zc,GIP_path_(_T("DCT"),g_nzSample) );

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
	GIP_save_IplImage_d( p_M,N,GIP_path_(_T("zero_crossing"),g_nzSample),zc,N,0x0 );

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
	unsigned char *codelets,*pvc,*V; 
	GIP_OBJECT *hBand;
	GIP_FLOATING *data,*B_data;
	double *patch,*pv,*d_temp,a,*dv;		//mean,devia,

	M=hROI->M;				N=hROI->N;
	ldu=N;					//ldu始终不变
//	ASSERT( M==M_normal && N==N_normal );
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
//	GIP_save_IplImage_d( hBand->M,hBand->N,".\\trace\\iris_band_.bmp",hBand->data,hBand->N,0x0 );3

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
			patch[i] = a/step_M;
		}
		_DCT_12_( patch,d_temp );					//1D-DCT
	}
	}
	GIP_OBJ_clear( hBand );					GIP_free( hBand );
//	__testing_zero_crossing( p_M,p_N*p_len,p_len,pv,0x0 );
	switch( hFV->type )	{
	case G_FV_ROI:
		dv = (double*)(hFV->V);
//		GIP_fourier_spectrum( hROI->M,hROI->N,hROI,dv,0x0 );
		GIP_MEMCOPY( dv,hROI->data,sizeof(double)*hROI->M*hROI->N );
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
		D_ouput( hFV->V_len,dv,GIP_path_( _T("dv"),g_nzSample) );
//		GIP_MEMCOPY( hFV->V,pv,sizeof(double)*p_M*p_N*p_len );
		break;
	default:
		codelets = GIP_alloc( sizeof(unsigned char)*p_M*p_N );
	//binary code from zero crossings of adjacent patch vectors
		for( i = 0; i < p_M; i++ )	{
		for( j = 0; j < p_N; j++ )	{
			pvc = codelets + i*p_N+j;
			*pvc=0x0;
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
	MASEK's 1-D gabor filer
	[REF]: Recognition of Human Iris Patterns for Biometric Identification	
	[REF]: Iris Recognition System based on Log-Gabor Coding

	Copyright 2008-present, Grusoft
	v0.1	cys
		2/24/2009	
*/
void _IRIS_Feature_MASEK_( GIP_OBJECT *hROI,G_FVector_ *hFV,int method,double *D_buffer )	{
	int r,c,byte,bit,M,N,ldu,unit,t;
	unsigned char *codelets,*V,phase; 
	GIP_FLOATING *data;
	double w,*re,*im,lenda,*feature;	
	GIP_GABOR *hGabor;

	M=hROI->M;				N=hROI->N;
	data=hROI->data;
	ldu=N;	
	unit=hFV->unit_size;

	hGabor = GIP_alloc( sizeof(GIP_GABOR) );
	lenda=32.0,				
	GIP_LOG_GABOR_init( hGabor,1,N,lenda,0.0,0.0,GIP_GABOR_FIELD,GIP_GABOR_1D|GIP_GABOR_FREQ );

	feature=GIP_alloc( sizeof(double)*M*N );
	codelets = GIP_calloc( unit,hFV->V_len );
	re=D_buffer;		im=re+N;
	for( r=0; r<M; r++ )	{	
		GIP_GABOR_filter( hGabor,1,N,data+r*N,re,im,0x0 );	
		for( c=0; c<N; c++ )	{
	//		re=hGabor->re[c];		im=hGabor->im[c];
	//		if( c>=N/2)		{
	//			im[c]=0.0;		re[c]=0.0;
	//		}
			w = atan2( im[c],re[c] );
			phase = (0<w&&w<=PI/2) ? 0xC0 : (PI/2<w&&w<=PI) ? 0x40 : (-PI<w&&w<=-PI*0.5) ? 0 : 0x80;
			feature[GIP_RC2POS( r,c,N )]=phase;
			byte = (r*N+c)*2/8;			bit=(r*N+c)*2%8;
			t = phase;
			phase = t>>bit;
			codelets[byte] |= phase;
		}		
	}

	V =	hFV->V;	
	GIP_MEMCOPY( V,codelets,unit*hFV->V_len );
	
	if( 0 )	{
		GIP_save_IplImage_d( M,N,GIP_path_( _T("feature"),g_nzSample ),feature,N,0x0 );
	}
	GIP_free( feature );

	GIP_free( codelets );
	GIP_GABOR_clear( hGabor );
	GIP_free( hGabor );
}

/*
[**]
	顺平打了三次电话，每次都说没事，真奇怪。。。
[**]

	only lef half part of ROI
	D_buffer[2*N]

	Copyright 2008-present, Grusoft
	v0.1	cys
		3/16/2009	
*/
void _IRIS_Feature_CYS_( double *control,GIP_OBJECT *hROI,G_FVector_ *hFV,int method,double *D_buffer )	{
	int r,c,byte,bit,M,N,ldu,unit,t,pos_1,pos_2;
	unsigned char *codelets,*V,phase; 
	GIP_FLOATING *data;
	double w,*re,*im,lenda,*feature;	
	GIP_GABOR *hGabor;

	M=hROI->M;				N=hROI->N;
	data=hROI->data;
	ldu=N;	
	unit=hFV->unit_size;

	hGabor = GIP_alloc( sizeof(GIP_GABOR) );
	lenda=control[GIRIS_LG_WAVELEN];			//32.0,				
	GIP_LOG_GABOR_init( hGabor,1,N,lenda,0.0,0.0,GIP_GABOR_FIELD,GIP_GABOR_1D|GIP_GABOR_FREQ );

	feature=GIP_alloc( sizeof(double)*M*N );
	codelets = GIP_calloc( unit,hFV->V_len );
	re=D_buffer;				im=re+N;
	for( r=0; r<M; r++ )	{	
		GIP_GABOR_filter( hGabor,1,N,data+r*N,re,im,0x0 );	
		for( c=0; c<N; c++ )	{
	//		re=hGabor->re[c];		im=hGabor->im[c];
			w = atan2( im[c],re[c] );
			phase = (0<w&&w<=PI/2) ? 0xC0 : (PI/2<w&&w<=PI) ? 0x40 : (-PI<w&&w<=-PI*0.5) ? 0 : 0x80;
			feature[GIP_RC2POS( r,c,N )]=phase;
			byte = (r*N+c)*2/8;			bit=(r*N+c)*2%8;
			t = phase;
			phase = t>>bit;
			codelets[byte] |= phase;
		}		
	}
/*	 
	for( r=0; r<M; r++ )	{	
		GIP_GABOR_filter( hGabor,1,N,data+r*N,re,im,0x0 );	
		for( c=0; c<N; c++ )	{
			pos_1=GIP_RC2POS( r,c,N );		
			pos_2= r==M-1 ? pos_1-N : pos_1+N;
			//pos_2= c==N-1 ? pos_1-1 : pos_1+1;
			byte = (r*N+c)/8;			bit=(r*N+c)%8;
			t = feature[pos_1]==feature[pos_2] ? 0x00 : 0x80;
			phase = t>>bit;
			codelets[byte] |= phase;
		}		
	}
*/
	V =	hFV->V;	
	GIP_MEMCOPY( V,codelets,unit*hFV->V_len );
	
	if( 1 )	{
		GIP_save_IplImage_d( M,N,GIP_path_( _T("feature"),g_nzSample ),feature,N,0x0 );
	}
	GIP_free( feature );

	GIP_free( codelets );
	GIP_GABOR_clear( hGabor );
	GIP_free( hGabor );
}

/*
	D_buffer[2*N]

	Copyright 2008-present, Grusoft
	v0.1	cys
		4/10/2009	
*/
void _IRIS_Feature_CYS_1_( double *control,GIP_OBJECT *hROI,G_FVector_ *hFV,int method,double *D_buffer )	{
	int r,c,byte,bit,M,N,ldu,unit,t,pos_1,pos_2,M_step,r_1;
	unsigned char *codelets,*V,phase; 
	GIP_FLOATING *data;
	double w,*re,*im,lenda,*feature,*temp;	
	GIP_GABOR *hGabor;

	M=hROI->M;				N=hROI->N;
	M_step=hFV->M_step;
	data=hROI->data;
	ldu=N;	
	unit=hFV->unit_size;

	hGabor = GIP_alloc( sizeof(GIP_GABOR) );
	lenda=control[GIRIS_LG_WAVELEN];			
	GIP_LOG_GABOR_init( hGabor,1,N,lenda,0.0,0.0,GIP_GABOR_FIELD,GIP_GABOR_1D|GIP_GABOR_FREQ );

	feature=GIP_alloc( sizeof(double)*M*N*2 );
	temp = feature+M*N;
	codelets = GIP_calloc( unit,hFV->V_len );
	re=D_buffer;				im=re+N;
	for( r=0; r<M/M_step; r++ )	{	
		if( M_step==0 )	{
			for( c=0; c<N; c++ )	{
				temp[c]=data[GIP_RC2POS(r,c,N)];
			}
		}	else	{
			for( c=0; c<N; c++ )	{
				temp[c]=0.0;
				for( r_1=r*M_step;r_1<(r+1)*M_step;r_1++ )	{
					temp[c] += data[GIP_RC2POS(r_1,c,N)];
				}
				temp[c]/=M_step;
			}
		}
//		GIP_GABOR_filter( hGabor,1,N,data+r*M_step*N,re,im,0x0 );	
		GIP_GABOR_filter( hGabor,1,N,temp,re,im,0x0 );	
		for( c=0; c<N; c++ )	{
	//		re=hGabor->re[c];		im=hGabor->im[c];
			w = atan2( im[c],re[c] );
			phase = (0<w&&w<=PI/2) ? 0xC0 : (PI/2<w&&w<=PI) ? 0x40 : (-PI<w&&w<=-PI*0.5) ? 0 : 0x80;
			feature[GIP_RC2POS( r,c,N )]=phase;
			byte = (r*N+c)*2/8;			bit=(r*N+c)*2%8;
			t = phase;
			phase = t>>bit;
			ASSERT( byte<=hFV->V_len );
			codelets[byte] |= phase;
		}		
	}

	V =	hFV->V;	
	GIP_MEMCOPY( V,codelets,unit*hFV->V_len );
	
	if( 1 )	{
		GIP_save_IplImage_d( M,N,GIP_path_( _T("feature"),g_nzSample ),feature,N,0x0 );
	}
	GIP_free( feature );

	GIP_free( codelets );
	GIP_GABOR_clear( hGabor );
	GIP_free( hGabor );
}

/*
	利用光斑边界急速衰减的特点

	1 权重以data[pos]为主，即光斑接近于最亮
[**]
	今天查出颈椎不好，没想到老的这么快。。。
	希望能再写四到五年。。。。
[**]
	I_buffer[3*M*N]

	[ISSUE-NEEDED]
	1 存在pupil之外有spot的例子：MMU\1\right\aevar3

	v0.1	cys
		3/21/2009
*/
int _IRIS_spot_remove( double *control,G_IRIS_ *hIris,int *I_buffer,double *D_buffer,int flag )	{
	GIP_OBJECT *hObj=hIris->hObj;
	int M=hObj->M,N=hObj->N,I_0,I_1,high,low,i,j,r,c,r1,c1,pos,pos_1,nSpot=0,R=2;
	float w[4],a[4],a_min,val,w_sum;
	double *g_I=D_buffer,gr,gc,g_max,dist;
	int *stack,*list,*map,seed,noR,nz,nzB,nz_1,top,neibor[]={-1,1,-N,N};
	GIP_FLOATING *data=hObj->data,area_0,wb;

	GIP_save_IplImage_d( M,N,GIP_path_( _T("spot_before"),g_nzSample ),data,N,0x0 );

	I_0=hIris->hIa->histo.I_0;				I_1=hIris->hIa->histo.I_1;
	area_0 = M*N*control[GIRIS_SPOT_AREA];			///100*3;
	if( area_0<1 )
		return;
	high = MAX(I_1-15,I_1*0.9);			low=I_0+(I_1-I_0)*0.3;
//	I_buffer=GIP_alloc( sizeof(int)*M*N*3 );
	stack=I_buffer;				
	list=stack+M*N;				map=list+M*N;
	for( i = 0; i < M*N; i++ )		{
		g_I[i]=0.0;
		map[i]=-1;
	}

	g_max = 0.0;
	nz = ( 2*R+1)*(2*R+1)-1;
	for( r=R; r<M-R; r++ )	{	//加权系数以data[pos]为主	
	for( c=R; c<N-R; c++ )	{
		pos=GIP_RC2POS( r,c,N );	
		if( r==69 && c==69 )
			i=0;	
		w_sum = 0.0;
		for( r1=r-R;	r1<=r+R; r1++ )		{
		for( c1=c-R;	c1<=c+R; c1++ )		{
			if( r1==r && c1==c )
				continue;
			pos_1=GIP_RC2POS( r1,c1,N );	
			dist = sqrt( (r1-r)*(r1-r)+(c1-c)*(c1-c) );
			w_sum += (data[pos]-data[pos_1])/dist;
		}
		}
		g_I[pos] = w_sum/nz+data[pos]-high;
//		gr = GIP_Mask_5_( pos,M,N,N,data,G_MASK_SPOT_R,0x0 );
//		gc = GIP_Mask_5_( pos,M,N,N,data,G_MASK_SPOT_C,0x0 );
//		g_I[pos] = (sqrt( gr*gr+gc*gc ))+data[pos]-high;
//		g_I[pos] = (log( gr*gr+gc*gc ))*data[pos];
		g_max = MAX( g_max,g_I[pos] );
	}
	}
	
	noR=0;
	for( r=1; r<M-1; r++ )	{	
	for( c=1; c<N-1; c++ )	{
		pos=GIP_RC2POS( r,c,N );
		if( data[pos]<data[pos-1] || data[pos]<data[pos+1] || data[pos]<data[pos-N] || data[pos]<data[pos+N] )
			continue;
		if( r==69 && c==69 )
			i=0;
//		gr = GIP_Mask_5_( pos,M,N,N,data,G_MASK_FIVE_R,0x0 );
//		gc = GIP_Mask_5_( pos,M,N,N,data,G_MASK_FIVE_C,0x0 );
//		g_I[pos] = (sqrt( gr*gr+gc*gc ))+data[pos];
		if( g_I[pos]<g_max*0.7 || data[pos]<high )			continue;
		if( map[pos] != -1 )		continue;
		seed = pos;
		top=0;			nz=0;			nz_1=0;
		stack[top]=seed;
		wb = 0.0;		nzB=0;
		map[pos] = ++noR;
		while( top>=0 )	{
			pos = stack[top--];	
			list[nz] = pos;
			nz++;		
			for( i = GIP_NEIBOR_LEFT; i <= GIP_NEIBOR_UP; i++ )	{
				pos_1 = pos + neibor[i];
				r1 = GIP_POS2R( pos_1,N );			c1 = GIP_POS2C( pos_1,N );
				if( r1<1 || r1>=M-1 || c1<1 || c1>=N-1)		continue;
				ASSERT( pos_1>=0 && pos_1<M*N ) ;
				if( map[pos_1]!=-1 )	
					continue;
				if( data[pos_1]>=high )			nz_1++;
				if( data[pos_1]>data[pos] )	{
					wb+=data[pos];		nzB++;
					continue;
				}
				stack[++top]=pos_1;
				map[pos_1]=noR;
				
			}
		}		
		ASSERT( nz<=M*N && nzB>0 /*&& nzB<nz*/ );
		wb /= nzB;
		if( nz>=area_0 || wb>low /*|| nz_1<2*/ )	{	//面积小，边界黑
			for( i = 0; i < nz; i++ )	
			{	pos = list[i];			map[pos]=-1;		}
			continue;
		}
		nSpot++;
		for( i = nz-1; i >= 0; i-- )	{		//线性插值显著改进了spot填充的效果
			pos = list[i];
			r1 = GIP_POS2R( pos,N );			c1 = GIP_POS2C( pos,N );
			if( r1==63 && c1==46 )
				j=0;
			w_sum=0.0;			a_min=FLT_MAX;
			for( j = GIP_NEIBOR_LEFT; j <= GIP_NEIBOR_UP; j++ )	{
				pos_1 = pos;
				while( map[pos_1]==noR )	{
					pos_1 += neibor[j];
				}
				w[j] = (pos_1-pos)/neibor[j];		ASSERT( w[j]>0 );
				a[j]=data[pos_1];		a_min=MIN(a_min,a[j]);
				w[j]=1.0/w[j];	w_sum += w[j];
			}
			val = 0.0;
			for( j = GIP_NEIBOR_LEFT; j <= GIP_NEIBOR_UP; j++ )	
				val += a[j]*w[j]/w_sum;
			data[pos] = MIN( data[pos],a_min );		
			//data[pos] = MIN( val,data[pos] );
			ASSERT( data[pos]>=0.0 && data[pos]<=255.0 );
		}
		for( i = 0; i < nz; i++ )	{
			pos = list[i];			map[pos]=-1;
		}
	}		
	}

	GIP_save_IplImage_d( M,N,GIP_path_( _T("spot_remove"),g_nzSample ),data,N,0x0 );

	return nSpot;
//	GIP_free( I_buffer );

}

/*
	注意：
		不修改hIris的数据

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/1/2008	
*/
void 	_IRIS_ouput_( double *control,G_IRIS_ *hIris,_TCHAR* sPath,int flag )	{
	GIP_OBJECT *hObj=hIris->hObj,*hOut;
	int M=hObj->M,N=hObj->N;

	hOut=GIP_OBJ_create( M,N,hObj->data,0x0 );
	GIP_Circle( hOut,hIris->p_r,hIris->p_c,hIris->p_radius,GIP_IMAG_ORIGIN );
	GIP_Circle( hOut,hIris->s_r,hIris->s_c,hIris->s_radius,GIP_IMAG_ORIGIN );
	GIP_save_IplImage_d( M,N,sPath,hOut->data,N,0x0 );

	GIP_OBJ_clear( hOut );
	GIP_free( hOut );
}

/*
[**]
	小马和豪哥，兰博和教练之间的关系：师徒 父子般的感情 及同志
[**]	

	注意:
		1 hIMG不被改变
		2 直接将find_iris_2_的返回值传给调用函数处理

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/1/2008	
	v0.2	cys
		3/21/2009
		采用多重搜索
*/
G_IRIS_ *G_IRIS_process( double *control,GIP_IMAGE *hIMG,int model,int flag,int *ret )	{
	float scale,w;
	G_IRIS_ *hIris;
	GIP_OBJECT *hObj,obj_0,obj_1;
//	GIP_IMAGE *hIMG_core=hIMG;
	int M,N,isTrace=1,I_0,I_1,low,high,n_limit,isScale=0,nSpot;
	double _ToI_0=control[GIRIS_T_INTN_0],_ToI_1=control[GIRIS_T_INTN_1],q_0=control[GIRIS_SMOOTH_Q];
	clock_t start;

	*ret=GIP_OK;

//	n_limit=5120;	
	n_limit=15000;	
	GIP_import_IplImage( hIMG,&(obj_0.M),&(obj_0.N),&(obj_0.data),0x0 );
	M=hIMG->height;				N=hIMG->width; 
	if( hIMG->height*hIMG->width>=n_limit )	{
		scale=sqrt(n_limit*1.0/hIMG->height/hIMG->width);
		M=G_DOUBLE2INT(hIMG->height*scale);		
		N=G_DOUBLE2INT(hIMG->width*scale);
	//	M &= -(1<<2);					N &= -(1<<2);
	}
	if( M!=hIMG->height || N!=hIMG->width )	{
		isScale=1;
		GIP_OBJ_init( &(obj_1),M,N,0x0 );
		GIP_OBJ_resize( &obj_0,&(obj_1),5,0 );
//		hIMG_core = GIP_IMG_resize( hIMG,M,N,0x0 );
#ifdef _GIRIST_DEVELOP_
//		G_PRINTF( _T("RESIZE:\tM:%d->%d,N:%d->%d\r\n"),hIMG->height,M,hIMG->width,N );	
#endif
	}	else	{
		GIP_OBJ_init( &(obj_1),M,N,0x0 );
		GIP_MEMCOPY( obj_1.data,obj_0.data,sizeof(GIP_FLOATING)*M*N );
//		hIMG_core = hIMG;
	}
//	GIP_import_IplImage( hIMG_core,&(obj_1.M),&(obj_1.N),&(obj_1.data),0x0 );
/*
*/	hIris=GIP_alloc( sizeof(G_IRIS_) );
	G_IRIS_init_( control,&(obj_1),hIris );
	hObj=hIris->hObj;
	nSpot = _IRIS_spot_remove( control,hIris,hIris->hIa->I_buffer,hIris->hIa->D_buffer,0x0 );
	I_0=hIris->hIa->histo.I_0;				I_1=hIris->hIa->histo.I_1;
	low=(int)(I_0+(I_1-I_0)*_ToI_0);		high=(int)(I_0+(I_1-I_0)*_ToI_1);
	if( isScale==1 )	
		control[GIRIS_SMOOTH_Q]=0.8;
	_IRIS_transform( control,hIris,GIP_M_SMOOTH );
//	*ret = _IRIS_find_iris_2_( control,hIris,low,high,0x0 );			g_nF_1++;
	*ret = _IRIS_find_iris_1_( control,hIris,0x0 );			g_nF_1++;
	if( isTrace==1 )		//	!!testing!!
		_IRIS_ouput_( control,hIris,GIP_path_(_T("CHT"),g_nzSample),0x0 );
/*	if( *ret !=GIP_OK )	{
		GIP_MEMCOPY( hObj->data,obj_1.data,sizeof(GIP_FLOATING)*M*N );
		_IRIS_transform( control,hIris,GIP_M_HISTO|GIP_M_SMOOTH );
		*ret = _IRIS_find_iris_1_( control,hIris,0x0 );			g_nF_2++;
	//	*ret = _IRIS_find_iris_2_( control,hIris,low,high,0 );				g_nF_2++;
		if( isTrace==1 )				//	!!testing!!
			_IRIS_ouput_( control,hIris,GIP_path_(_T("CHT"),g_nzSample),0x0 );
	}*/

	if( isScale==1 )	{
		G_IRIS_ t_iris;

		GIP_MEMCOPY( &t_iris,hIris,sizeof(G_IRIS_) );
//		cvReleaseImage( &hIMG_core );
		G_IRIS_clear( hIris );
		G_IRIS_init_( control,&(obj_0),hIris );
		hObj = hIris->hObj;
//		high=(int)(I_0+(I_1-I_0)*(_ToI_1+0.2) );
//		_IRIS_spot_remove( control,hIris,hIris->hIa->I_buffer,hIris->hIa->D_buffer,0x0 );
	start=clock( );
		_IRIS_transform( control,hIris,GIP_M_SMOOTH );
	hIris->x = (clock( )-start)*1.0/CLOCKS_PER_SEC;
		hIris->p_r=t_iris.p_r/scale;		hIris->p_c=t_iris.p_c/scale;		
		hIris->p_radius=t_iris.p_radius/scale;			hIris->p_match=t_iris.p_match;		
		hIris->s_r=t_iris.s_r/scale;		hIris->s_c=t_iris.s_c/scale;		
		hIris->s_radius=t_iris.s_radius/scale;			hIris->s_match=t_iris.s_match;			
		if( 0 )			
			_IRIS_ouput_( control,hIris,GIP_path_(_T("CHT"),g_nzSample),0x0 );
		control[GIRIS_SMOOTH_Q]=q_0;
		*ret = _IRIS_find_iris_refine_( control,hIris,low,high,0x0,*ret );		
		//	GIP_Parabolic( hEdge,hIris->ul_a,hIris->ul_h,hIris->ul_k,GIP_IMAG_ORIGIN );
		if( 1 )				//	!!testing!!
			_IRIS_ouput_( control,hIris,GIP_path_(_T("CHT"),g_nzSample),0x0 );

	}
	M=hObj->M;					N=hObj->N; 
	GIP_MEMCOPY( hObj->data,obj_0.data,sizeof(GIP_FLOATING)*hObj->M*hObj->N );		//复原数据
	w=hIris->p_radius;		
	if( *ret==GIRIST_ERROR_BOUNDARY )	{			//based on circular hash transform
		goto GIP_exit;		
	}else	{
		if( hIris->p_match<control[GIRIS_T_MTCH_P]/(nSpot+1) )	{			
			G_PRINTF( _T("\t!!!Warning:\tpupile(%g,%g,%g),p_match=%4.2g\r\n"),hIris->p_r,hIris->p_c,hIris->p_radius,hIris->p_match );
			*ret = GIRIST_WARN_INNER_BOUNDARY;			
		}
		if( hIris->s_match<control[GIRIS_T_MTCH_S] )	{
			G_PRINTF( _T("\t!!!Warning:\tsclera(%g,%g,%g),s_match=%4.2g\r\n"),hIris->s_r,hIris->s_c,hIris->s_radius,hIris->s_match );
			*ret = GIRIST_WARN_OUTER_BOUNDARY;		
		}
		if( hIris->s_radius<=hIris->p_radius )	{
			G_PRINTF( _T("\t!!!s_radius<=p_radius:\tp_radius=%4.2g,s_radius=%4.2g\r\n"),hIris->s_radius,hIris->p_radius );
			*ret = GIRIST_ERROR_BOUNDARY;	
			goto GIP_exit;		
		}else if( hIris->p_r-w<0||hIris->p_r+w>=M || hIris->p_c-w<0||hIris->p_c+w>=N )	{
			G_PRINTF( _T("\t!!!Pupil is out of Image:\tp_radius=%4.2g,r=%4.2g,c=%4.2g\r\n"),hIris->p_radius,hIris->p_r,hIris->p_c );
			*ret = GIRIST_ERROR_BOUNDARY;	
			goto GIP_exit;		
		}	
	}
	_IRIS_normalize_( hIris,0x0 );
//	_IRIS_contour_normalize_( hIris,0x0 );

GIP_exit:
	GIP_OBJ_clear( &obj_0 );
	GIP_OBJ_clear( &obj_1 );
	return hIris;
}

/*
	根据hEntry的status确定iris的参数已确定，并normalize

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/8/2009	
*/
G_IRIS_ *G_IRIS_on_entry( double *control,GIP_LIB_entry *hEntry,int model,int flag,int *ret )	{
	G_IRIS_ *hIris=GIP_NULL;
	GIP_OBJECT *hObj,obj_0;
	int M,isTrace=1,I_0,I_1,low,high,n_limit=90000,iscolor=0x0;
	double _ToI_0=control[GIRIS_T_INTN_0],_ToI_1=control[GIRIS_T_INTN_1],s_N=1.0/48/512;
	GIP_IMAGE *hIMG_0;
	GIa_HANDLE *hIa = GIP_NULL;

	*ret=GIP_ERROR_internal_error;

	hIMG_0=GIP_LoadImage( hEntry->sFilePath,iscolor,0x0 );		ASSERT( hIMG_0!=0x0 );
	GIP_import_IplImage( hIMG_0,&(obj_0.M),&(obj_0.N),&(obj_0.data),0x0 );
	if( hEntry->status==G_ENTRY_OK )	{
		hIris=GIP_alloc( sizeof(G_IRIS_) );
		G_IRIS_init_( control,&(obj_0),hIris );
		hIris->p_r=hEntry->p_r;			hIris->p_c=hEntry->p_c;				hIris->p_radius=hEntry->p_radius;
		hIris->s_r=hEntry->s_r;			hIris->s_c=hEntry->s_c;				hIris->s_radius=hEntry->s_radius;
		hIris->p_match=hEntry->p_match;			hIris->s_match=hEntry->s_match;		

		_IRIS_normalize_( hIris,0x0 );
//	_IRIS_contour_normalize_( hIris,0x0 );
	}else	if( 1 )	{
		hIris = G_IRIS_process( control,hIMG_0,0x0,flag,ret );
		hEntry->p_r=hIris->p_r;		hEntry->p_c=hIris->p_c;		hEntry->p_radius=hIris->p_radius;		
		hEntry->s_r=hIris->s_r;		hEntry->s_c=hIris->s_c;		hEntry->s_radius=hIris->s_radius;		
		hEntry->p_match=hIris->p_match;						hEntry->s_match=hIris->s_match;	
	}

	hObj=hIris->hObj;
	I_0=hIris->hIa->histo.I_0;				I_1=hIris->hIa->histo.I_1;
	low=(int)(I_0+(I_1-I_0)*_ToI_0);		high=(int)(I_0+(I_1-I_0)*_ToI_1);
	hIa = hIris->hIa;	

//	G_PRINTF( _T("%-4d %3g %3g %4.2g %4.2g\t %-5.4g %-5.3g %-5.2g %-4.2g\r\n"),1,hIris->p_radius,hIris->s_radius,
//			hIris->p_match,hIris->s_match,hIa->mean,hIa->deviation,hIris->nzNoise*s_N,hIris->x );	
	*ret=GIP_OK;

GIP_exit:
	GIP_OBJ_clear( &obj_0 );

	cvReleaseImage( &hIMG_0 );

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
void G_IRIS_save_training( GIP_LIB *hLib,G_FVector_ *fv,TCHAR *sLibPath,char *sample,int flag )	{
	int nEntry=hLib->nEntry,nSample=hLib->nSample,N;
	FILE *fp;

	fp=_tfopen( sLibPath,_T("wb") );
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

	fp=_tfopen( sLibPath,_T("rb") );
	ASSERT( fp!=0 );

	fread( fv,sizeof(G_FVector_),1,fp );	
	fv->V= GIP_alloc( fv->unit_size*fv->V_len );
	fread( hLib,sizeof(GIP_LIB),1,fp );				
	hLib->fv_dat=GIP_NULL;	
/*	hLib->W_ptr=GIP_NULL;				hLib->W_no=GIP_NULL;	
	hLib->W_ptr=GIP_alloc( sizeof(int)*(hLib->nClass+1) );		
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
[**]
	tgbus被封，郁闷。。。
[**]	

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/25/2009	
*/
void _IRIS_Feature_Vector( double *control,GIP_OBJECT *hROI,G_FVector_ *hFV,int flag )	{
	int N=hROI->N;
	double *D_buffer;
//	double s=1.0/CLOCKS_PER_SEC;
//	clock_t start_0;

	D_buffer=GIP_alloc( sizeof(double)*N*2 );
//	start_0=clock( );
	switch( hFV->type )	{
	case G_FV_MASEK:
		_IRIS_Feature_MASEK_( hROI,hFV,0x0,D_buffer ); 
		break;
	case G_FV_CYS:
		_IRIS_Feature_CYS_( control,hROI,hFV,0x0,D_buffer ); 
		break;
	case G_FV_CYS_1:
		_IRIS_Feature_CYS_1_( control,hROI,hFV,0x0,D_buffer ); 
		break;
	case G_FV_DCT_CODE:
	case G_FV_ROI:
	case G_FV_COEFFICIENT:
		_IRIS_Feature_DCT_( GIP_NULL,hROI,hFV,0x0 ); 
		break;
	case G_FV_MALI:
		_IRIS_Feature_MALI_( GIP_NULL,hFV,0x0 ); 
		break;
	default:
		ASSERT( 0 );
		break;
	}
	GIP_free( D_buffer );

//	s = (clock( )-start_0)*s;
}

/*
[**]
	圣诞快乐,嘻嘻
[**]	

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/24/2008	
	v0.2	cys
		2/27/2009		由control控制流程	

void G_IRIS_program( char *sTrainPath,char *sLibPath,int flag )	{
	int nSample=0,N,i,nRet,iscolor=0x0,nOK=0,ret,cls_old=-1,mod;		
	int M_normal,N_normal;
	G_IRIS_ *hIris;
	GIP_LIB lib;
	G_FVector_ fv,*hfv=&fv;
	GIP_IMAGE *src_image,*img;
	GIP_LIB_entry *hEntry;
	GIa_HANDLE* hIa;
	clock_t start=clock( ),start_0;
	double s=1.0/CLOCKS_PER_SEC,*w_sort,s_N=1.0/48/512,Mp_0=100,Ms_0=100;
	unsigned char *sample;

	GIP_init_heapinfo( );

	mod = (int)(control[GIRIS_MODE_]);
	M_normal=(int)(control[GIRIS_M_NORMAL]),		N_normal=(int)(control[GIRIS_N_NORMAL]);

	start_0=clock( );
	if( mod==GIRIS_MODE_VERIFY )	{
	//	G_IRIS_load_training( &lib,&fv,"..\\iris_training_1.data",&sample,0x0 );
		GIP_LIB_load( &lib,hfv,sLibPath, 0x0 );
	}else if( mod==GIRIS_MODE_TRAINING )		{
//		N=M_normal*N_normal/16;
		G_FVector_init( &fv,M_normal,N_normal,(int)(control[GIRIS_FV_SIFT]),0x0 );
		nRet = GIP_LIB_init( &lib,sTrainPath,&fv,control,0x0 );			ASSERT(nRet==0x0);
//		lib.fv_n_0 = N;
		nSample=lib.nEntry;
		N = hfv->V_len;
	//	nSample=1;
		sample=lib.fv_dat;
//		sample=(unsigned char*)GIP_alloc( hfv->unit_size*nSample*N );		
		w_sort=GIP_alloc( sizeof(double)*nSample );
		G_PRINTF( _T(" no  R_0 R_1 p_mt s_mt \t mean devia noise  x  \r\n") );
		for( i = 0; i < nSample; i++ )	{
			g_nzSample = nOK;
			hEntry = lib.arrEntry+i;
			if( hEntry->cls!=cls_old )	{
				G_PRINTF( _T("%s\r\n"),hEntry->sFilePath );
				cls_old = hEntry->cls;
			}
			start=clock( );
			src_image=GIP_LoadImage(hEntry->sFilePath,iscolor,0x0);		ASSERT( src_image!=0x0 );
			img = cvCloneImage( src_image );
//			hObj = GIP_IMG2OBJ( src_image,0x0 );
//			GIP_MASK_testing( hObj,0x0 );
//			GIP_OBJ_clear( hObj );		GIP_free( hObj );
			hIris = G_IRIS_process( img,src_image,0x0,flag,&ret );
			hIa = hIris->hIa;
			if( ret==GIP_OK )	{
				hEntry->status=1;
				_IRIS_Feature_Vector( hIris,hfv,0x0 );
				memcpy( sample+nOK*N*hfv->unit_size,hfv->V,hfv->unit_size*N );
				w_sort[nOK] = hIris->s_match;
				Mp_0=MIN(Mp_0,hIris->p_match);		Ms_0=MIN(Ms_0,hIris->s_match);
				nOK++;
			}
			G_PRINTF( _T("%-4d %3g %3g %4.2g %4.2g\t %-5.4g %-5.3g %-5.2g %-4.2g\r\n"),i+1,hIris->p_radius,hIris->s_radius,
				hIris->p_match,hIris->s_match,hIa->mean,hIa->deviation,hIris->nzNoise*s_N,hIris->x );	
			G_IRIS_clear( hIris );
			GIP_free( hIris );		hIris=GIP_NULL;
			cvReleaseImage( &src_image );
			cvReleaseImage( &img );
//			G_PRINTF( "I_%d:\ttime=%g\r\n",i+1,(clock( )-start)*s );	
//			break;
		}
		lib.nSample = nOK;
		G_PRINTF( _T("\r\nFINISH:fail=%g(%d),time=%g,Mp_0=%5.3g,Ms_0=%5.3g\r\n"),(nSample-nOK)*1.0/nSample,nSample-nOK,(clock( )-start_0)*s,
			Mp_0,Ms_0 );
	//	G_IRIS_save_training( &lib,&fv,"..\\iris_training_1.data",sample,0x0 );
	//	GIP_LIB_reducing( &lib,sample, 0x0 );
	//	_sort_lib_( &lib,w_sort,GIP_NULL,0x0 );
		GIP_free( w_sort );
	//	GIP_LIB_scan( &lib,0x0 );
		GIP_LIB_save( &lib,sLibPath, 0x0 );
	}
	if( 1 )
		GIP_LIB_discriminate( &lib,1,0x0 );

	G_FVector_clear( &fv );
	GIP_LIB_clear( &lib );
//	GIP_free( sample );	
//	GIP_GABOR_clear( &gabor_1 );			GIP_GABOR_clear( &gabor_2 );

	GIP_exit_heapinfo( );

}
*/

/*
	from G_IRIS_program;根据GiristApp的要求修改

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/2/2009	
		3/17/2009	分离出G_IRIS_coding
*/
int  G_IRIS_location( GIP_LIB *hLib,int flag )	{
	GU_CTRL_DATA *hUserCtrl=g_hUserCtrl;
	int nSample=0,N,i,iscolor=0x0,nOK=0,ret,cls_old=-1,mod;	
	int M_normal,N_normal;
	G_IRIS_ *hIris;
	G_FVector_ *hfv=hLib->hfv;
	GIP_IMAGE *src_image,*img;
	GIP_LIB_entry *hEntry;
	GIa_HANDLE* hIa;
	clock_t start=clock( ),start_0;
	double s=1.0/CLOCKS_PER_SEC,*w_sort,s_N=1.0/48/512,Mp_0=100,Ms_0=100,*control;
//	unsigned char *sample;

	ASSERT( hUserCtrl!=GIP_NULL && hLib!=GIP_NULL && hfv!=NULL );
//	GU_SetUserCtrl( hUserCtrl,_T("Begin...\r\n"),-1,0x0 );
	control = hLib->user_control;
	mod = (int)(control[GIRIS_MODE_]);
	M_normal=(int)(control[GIRIS_M_NORMAL]),		N_normal=(int)(control[GIRIS_N_NORMAL]);
	g_nF_1=0,	g_nF_2=0,	g_nF_3=0;

	start_0=clock( );
//	G_FVector_init( hfv,M_normal,N_normal,0x0 );
	nSample=hLib->nEntry;
	N = hfv->V_len;
//	nSample=1;
//	sample=hLib->fv_dat;
	w_sort=GIP_alloc( sizeof(double)*nSample );
	GU_SetUserCtrl( hUserCtrl,GIP_NULL,0,0x0 );
	G_PRINTF( _T(" no  R_0 R_1 p_mt s_mt \t mean devia noise  x  \r\n") );
	for( i = 0; i < nSample; i++ )	{
		if( hUserCtrl->bTerminate==TRUE )
			break;
		g_nzSample = nOK;
		hEntry = hLib->arrEntry+i;
		if( hEntry->cls!=cls_old )	{
			G_PRINTF( _T("%s\r\n"),hEntry->sFilePath );
			cls_old = hEntry->cls;
		}
		start=clock( );
		src_image=GIP_LoadImage( hEntry->sFilePath,iscolor,0x0 );		ASSERT( src_image!=0x0 );
		img = cvCloneImage( src_image );
//			hObj = GIP_IMG2OBJ( src_image,0x0 );
//			GIP_MASK_testing( hObj,0x0 );
//			GIP_OBJ_clear( hObj );		GIP_free( hObj );
		hIris = G_IRIS_process( control,img,0x0,0x0,&ret );
		hIa = hIris->hIa;
		switch( ret )	{
		case GIP_OK: 
		case GIRIST_WARN_INNER_BOUNDARY:
		case GIRIST_WARN_OUTER_BOUNDARY:
			hEntry->status=ret==GIP_OK ? G_ENTRY_OK : G_ENTRY_WARNING;
			hEntry->p_c=hIris->p_c;		hEntry->p_r=hIris->p_r;			hEntry->p_radius=hIris->p_radius;
			hEntry->s_c=hIris->s_c;		hEntry->s_r=hIris->s_r;			hEntry->s_radius=hIris->s_radius;
			hEntry->x = _IRIS_angle_( hIris,0x0 );
			hEntry->p_match=hIris->p_match;		hEntry->s_match=hIris->s_match;	
//			_IRIS_Feature_Vector( hIris,hfv,0x0 );
//			memcpy( sample+i*N*hfv->unit_size,hfv->V,hfv->unit_size*N );
			w_sort[nOK] = hIris->s_match;
			Mp_0=MIN(Mp_0,hIris->p_match);		Ms_0=MIN(Ms_0,hIris->s_match);
			nOK++;
			break;
		case GIRIST_ERROR_BOUNDARY: 
			hEntry->status=G_ENTRY_FAILED;
			break;
		default:
			break;
		}

		G_PRINTF( _T("%-4d %3g %3g %4.2g %4.2g\t %-5.4g %-5.3g %-5.2g %-4.2g\r\n"),i,hIris->p_radius,hIris->s_radius,
			hIris->p_match,hIris->s_match,hIa->mean,hIa->deviation,hIris->nzNoise*s_N,hIris->x );	
		G_IRIS_clear( hIris );
		GIP_free( hIris );		hIris=GIP_NULL;
		cvReleaseImage( &src_image );
		cvReleaseImage( &img );
		GU_SetUserCtrl( hUserCtrl,GIP_NULL,i+1,0x0 );
	}
	hLib->nSample = nOK;
	if( i==nSample )	{	//正常结束
		G_PRINTF( _T("\r\nFINISH:fail=%g(%d),time=%g,Mp_0=%5.3g,Ms_0=%5.3g\r\n"),
			(nSample-nOK)*1.0/nSample,nSample-nOK,(clock( )-start_0)*s,Mp_0,Ms_0 );
		G_PRINTF( _T("F_1=%d,F_2=%d,F_3=%d \r\n"),g_nF_1,g_nF_2,g_nF_3 );
	}
	hLib->L_spd = (clock( )-start_0)*s/hLib->nSample;
//	GIP_LIB_reducing( &lib,sample, 0x0 );
//	_sort_lib_( &lib,w_sort,GIP_NULL,0x0 );
	GIP_free( w_sort );

	ret=GIP_OK;
	return ret;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/17/2009	
*/
int  G_IRIS_coding( GIP_LIB *hLib,int flag )	{
	GU_CTRL_DATA *hUserCtrl=g_hUserCtrl;
	int nEntry=hLib->nEntry,N,i,iscolor=0x0,nOK=0,ret,cls_old=-1,mod;	
	G_FVector_ *hfv=hLib->hfv;
	GIP_LIB_entry *hEntry;
	clock_t start=clock( ),start_0;
	double s=1.0/CLOCKS_PER_SEC,*control= hLib->user_control;
	G_U8 *sample;
	G_IRIS_ *hIris=GIP_NULL;

	ASSERT( hUserCtrl!=GIP_NULL && hLib!=GIP_NULL && hfv!=NULL );
//	GU_SetUserCtrl( hUserCtrl,_T("Begin...\r\n"),-1,0x0 );
	mod = (int)(control[GIRIS_MODE_]);

	start_0=clock( );
	N = hfv->V_len;
	if( hLib->fv_dat!=GIP_NULL )	GIP_free( hLib->fv_dat );
	hLib->fv_dat=(G_U8*)GIP_alloc( hfv->unit_size*hfv->V_len*hLib->nEntry );		
	sample=hLib->fv_dat;
	GU_SetUserCtrl( hUserCtrl,GIP_NULL,0,0x0 );
	for( i = 0; i < nEntry; i++ )	{
		if( hUserCtrl->bTerminate==TRUE )
			break;
		g_nzSample = nOK;
		hEntry = hLib->arrEntry+i;
		if( hEntry->status==G_ENTRY_OK )	{
			hIris = G_IRIS_on_entry( control,hEntry,0x0,flag,&ret );
			ASSERT( hIris!=GIP_NULL );
			if( hIris!=GIP_NULL )	{
				_IRIS_Feature_Vector( control,&(hIris->normal),hfv,0x0 );
				memcpy( sample+i*N*hfv->unit_size,hfv->V,hfv->unit_size*N );
				G_IRIS_clear( hIris );			GIP_free( hIris );
				hIris=GIP_NULL;
			}else	{
			}
		}

		GU_SetUserCtrl( hUserCtrl,GIP_NULL,i+1,0x0 );
	}

	s = (clock( )-start_0)*s;

	ret=GIP_OK;

	return ret;
}

/*
	基于GIP_LIB输出信息到文件
	与G_IRIS_build相比，尽量减少计算 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/13/2009	
*/
int G_IRIS_LIB_info( GIP_LIB *hLib,TCHAR* sFilePath,int flag )	{
	int nEntry=hLib->nEntry,nSample=hLib->nSample,i,nRet,cls_old=-1;	
	int M_normal,N_normal,mod;
//	G_IRIS_ *hIris;
//	G_FVector_ *hfv=hLib->hfv;
//	GIP_IMAGE *src_image;
	GIP_LIB_entry *hEntry;
	double s=1.0/CLOCKS_PER_SEC,*w_sort,s_N=1.0/48/512,Mp_0=100,Ms_0=100,*control=hLib->user_control;
	FILE *fp=_tfopen( sFilePath,_T("w") ); 

	ASSERT( hLib!=GIP_NULL );
	mod = (int)(control[GIRIS_MODE_]);
	M_normal=(int)(control[GIRIS_M_NORMAL]),		N_normal=(int)(control[GIRIS_N_NORMAL]);
	_ftprintf( fp,_T("mode=%d,M_normal=%d,N_normal=%d\r\n"),mod,M_normal,N_normal );
	_ftprintf( fp,_T("\r\n\"%s\"\r\n"),hLib->sPath );
	_ftprintf( fp,_T("nEntry=%d,nClass=%d\r\n"),nEntry,hLib->nClass );
	_ftprintf( fp,_T("nOK=%d,nFail=%d\r\n"),nSample,nEntry-nSample );

	w_sort=GIP_alloc( sizeof(double)*nSample );
	_ftprintf( fp,_T("\r\n no  p_r p_c R_p s_r s_c R_s \r\n") );
	for( i = 0; i < nEntry; i++ )	{
		hEntry = hLib->arrEntry+i;
		if( hEntry->cls!=cls_old )	{
			_ftprintf( fp, _T("%s\r\n"),hEntry->sFilePath );
			cls_old = hEntry->cls;
		}
	//	src_image=GIP_LoadImage( hEntry->sFilePath,iscolor,0x0 );		ASSERT( src_image!=0x0 );
		if( hEntry->status==G_ENTRY_OK )	{
		}else	{
		//	hIris = G_IRIS_process( img,src_image,0x0,0x0,&ret );
		}

		_ftprintf( fp,_T("%-4d %3g %3g %3g %3g %3g %3g \r\n"),i+1,
			hEntry->p_r,hEntry->p_c,hEntry->p_radius,hEntry->s_r,hEntry->s_c,hEntry->s_radius
		 );	
	//	cvReleaseImage( &src_image );
	}
	fclose( fp );
//	_sort_lib_( &lib,w_sort,GIP_NULL,0x0 );
	GIP_free( w_sort );

	nRet=GIP_OK;
	return nRet;
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
void _IRIS_matching_( G_IRIS_ *hIris,G_FVector_ *hFv,GIP_LIB *hLib,unsigned char *sample,int flag )	{
	int i,nSample,no_min,V_len=hFv->V_len,bpSF=hFv->bpSF;
	int d_type=GIP_DISTANCE_EUCLID;
	unsigned char *V=hFv->V,*temp;
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
	int nSample=0,iscolor=0x0,ret;	//CV_LOAD_IMAGE_GRAYSCALE=0;;	
	G_IRIS_ *hIris;
	G_FVector_ fv;
	GIP_LIB lib;
	GIP_IMAGE *src_image,*img;
	clock_t start=clock( );
	double s=1.0/CLOCKS_PER_SEC;
	unsigned char *sample;

	GIP_init_heapinfo( );

	G_IRIS_load_training( &lib,&fv,"..\\iris_training_1.data",&sample,0x0 );

	src_image=GIP_LoadImage( sFilePath,iscolor,0x0 );		ASSERT( src_image!=0x0 );
	img = cvCloneImage( src_image );
//	hIris = G_IRIS_process( img,src_image,0x0,flag,&ret );
	if( ret==GIP_OK )	{
		//	_IRIS_Feature_Vector( hIris,&fv,0x0 );
		_IRIS_matching_( hIris,&fv,&lib,sample,0x0 );
	}

	G_IRIS_clear( hIris );
	GIP_free( hIris );		hIris=GIP_NULL;
	cvReleaseImage( &src_image );
	cvReleaseImage( &img );
	G_PRINTF( _T("M_:\ttime=%g\r\n"),(clock( )-start)*s );	

	GIP_LIB_clear( &lib );
	GIP_free( sample );	
	G_FVector_clear( &fv );			//G_FVector_clear( &fv_lib );

	GIP_exit_heapinfo( );

}

/*
	load control parameter

	v0.1	cys
		2/27/2009		
*/
int G_IRIS_init_control( TCHAR *sControlPath,double *control,int flag )	{
	int len,no,ret=GIRIST_ERROR_INIT_CONTROL;
	FILE *fp;
//	double a;
	TCHAR sLine[MAX_LINE],seps[]=_T(" ,\t\n="),*token;

	GIP_MEMCLEAR( control,sizeof(double)*GIP_USER_CONTROL );
//一般固定为(48,512)
	control[GIRIS_M_NORMAL]=48;		
	control[GIRIS_N_NORMAL]=512;

	if( sControlPath==GIP_NULL  ||						//缺省设置
		(fp=_tfopen( sControlPath,_T("r") ))==NULL )	{		
		control[GIRIS_MODE_]=1.0;
//		control[GIRIS_T_GRAD]=3.0;
		control[GIRIS_T_GRAD]=5.0;
		control[GIRIS_T_INTN_0]=0.0;
		control[GIRIS_T_INTN_1]=0.8;
//		control[GIRIS_T_INTN_1]=0.6;
		control[GIRIS_T_MTCH_S]=0.25;
		control[GIRIS_T_MTCH_P]=0.6;
		control[GIRIS_FV_SIFT]=6;
		control[GIRIS_LG_WAVELEN]=32;
		control[GIRIS_GRAD_BIAS]=PI/18;
		control[GIRIS_SPOT_AREA]=0.03;
		control[GIRIS_T_R0]=0.05;
		control[GIRIS_SMOOTH_Q]=1.0;

		if( sControlPath==GIP_NULL )	{
			ret=GIP_OK;		goto GIP_exit;
		}else
			goto _FIO_EXIT_;
	} 

	while( fReadLine( fp,sLine ) != 0 )	{
		len=(int)(_tcslen( sLine ));
		if( len==0 )			continue;
		if( sLine[0]==_T('c') )			
		{	G_PRINTF( _T("%s\r\n") ,sLine );			continue;			}
		if( sLine[0]==_T('!')  )		continue;			//注释行		
		no = -1;
		token = _tcstok( sLine, seps );
		if( _tcscmp (token,_T("GIRIS_MODE_") )==0x0 )		{
			no = GIRIS_MODE_;
		}else if( _tcscmp (token,_T("GIRIS_T_GRAD") )==0x0 )		{
			no = GIRIS_T_GRAD;
		}else if( _tcscmp (token,_T("GIRIS_T_INTN_0") )==0x0 )		{
			no = GIRIS_T_INTN_0;
		}else if( _tcscmp (token,_T("GIRIS_T_INTN_1") )==0x0 )		{
			no = GIRIS_T_INTN_1;
		}else if( _tcscmp (token,_T("GIRIS_T_MTCH_S") )==0x0 )		{
			no = GIRIS_T_MTCH_S;
		}else if( _tcscmp (token,_T("GIRIS_T_MTCH_P") )==0x0 )		{
			no = GIRIS_T_MTCH_P;
		}else if( _tcscmp (token,_T("GIRIS_FV_SIFT") )==0x0 )		{
			no = GIRIS_FV_SIFT;
		}else if( _tcscmp (token,_T("GIRIS_LG_WAVELEN") )==0x0 )		{
			no = GIRIS_LG_WAVELEN;
		}

		if( no != -1 )	{
			token = _tcstok( NULL, seps );
			if( _stscanf( token,_T("%lf"),control+no )!=1 )
			{	goto _FIO_EXIT_;	}
		}
	}
	ret=GIP_OK;
_FIO_EXIT_:
	if( ret==0 )	{
		if( fclose( fp )!=0 )
			ret = GIRIST_ERROR_INIT_CONTROL;
//		G_PRINTF( );				//输出控制参数
	}
//	if( ret!=0 )
//		GIP_ERROR( _T("G_IRIS_init_control"),ret, _T("") );
	
GIP_exit:
	return ret;
}
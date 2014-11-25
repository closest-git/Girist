#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/gip_util.h"
#include "../util/gip_ss.h"
#include "../util/gip_imag_transform.h"
#include "../util/Compact_DataStruc.h"
#include "GIa_basic.h"

/*
	temp[2*nz]

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/4/2006	
*/	
void GIa_Histo_init( GIP_OBJECT *hObj,GIa_HISTOGRAM *histo,int type,int *temp,int flag )	{
	int i,L,nz,M=hObj->M,N=hObj->N,MN=M*N,r,ldu=hObj->N;
	int neibors[8],n,*map,*stack,pos,pos_1,top;
	GIP_FLOATING *data=hObj->data;
	GIa_Histo_Entry *hEntry;
	double s;

	memset( histo,0x0,sizeof(GIa_HISTOGRAM) );
	histo->I_limit=0xff;
	L = histo->I_limit+1;

	histo->entrys=GIP_alloc( sizeof(GIa_Histo_Entry)*L );	
	memset( histo->entrys,0x0,sizeof(GIa_Histo_Entry)*L );	
	for( i = 0; i < L; i++ )	{
		hEntry = histo->entrys+i;
		hEntry->level=i;
	}	
	
	map = temp;			stack=temp+MN;
	for( i = 0; i < MN; i++ )	{
		r = GIP_INTENSITY( data[i] );		//(int)(data[i]+0.5);
		hEntry = histo->entrys+r;		ASSERT( hEntry->level==data[i] );
		hEntry->r += GIP_POS2R( i,ldu );	hEntry->c += GIP_POS2C( i,ldu );
		hEntry->nz++;
		map[i] = r;
	}

	for( i = 0; i < MN; i++ )	{
		r = map[i];
		if( r==-1 )		continue;
		hEntry = histo->entrys+r;		ASSERT( hEntry->level==r );
		top=0;
		stack[top]=i;		
		while( top>=0 )	{
			pos = stack[top--];	
			nz = GIP_Neibors( pos,neibors,M,N,ldu,0x0 );
			for( n = 0; n < nz; n++ )	{
				pos_1 = neibors[n];
				if( map[pos_1]!=r )	continue;
				map[pos_1]=-1;
				stack[++top]=pos_1;				
			}
		}
		hEntry->nBlob++;
	}
	histo->nPeak = 0;
	histo->I_0 = -1;			histo->I_1 = -1;
	for( i = 0; i < L; i++ )	{
		hEntry = histo->entrys+i;
		nz = hEntry->nz;
		if( nz==0 )		continue;
		if( histo->I_0 == -1 )		histo->I_0=i;
		if( i>histo->I_1 )			histo->I_1=i;
		s = nz*1.0/MN;
		if( s > histo->r_peak )
		{	histo->I_peak=i;		histo->r_peak=s;	}
		if( (i==0 || nz>histo->entrys[i-1].nz) && (i==L-1 || nz>histo->entrys[i+1].nz) )	{
			histo->nPeak++;
		}
		hEntry->r /= hEntry->nz;		hEntry->c /= hEntry->nz;
	}

	if( 1 )	{
		double a = 0.0,b=0.0;
		for( i = 0; i < L; i++ )	{
			hEntry = histo->entrys+i;
			a += hEntry->level*hEntry->nz;
		}
		for( i = 0; i < MN; i++ )	b += data[i];
		ASSERT( a==b );
	}
	if( 1 )	{
		int j;
		GIP_OBJECT obj,*hObj=&obj;
		GIP_OBJ_init( hObj,L,100,0x0 );
		pos = 0;		
		n=histo->entrys[histo->I_peak].nz;
		for( i = 0; i < L; i++ )	{
			hEntry = histo->entrys+i;
			nz = (int)(100.0*hEntry->nz/n+0.5);
			for( j = 0; j < nz; j++ )	hObj->data[pos++] = 0;
			for( j = nz; j < 100; j++ )	hObj->data[pos++] = 255;
		}
//		GIP_save_IplImage_d( hObj->M,hObj->N,GIP_path_(_T("histo"),L),hObj->data,hObj->N,0x0 );
		GIP_OBJ_clear( hObj );
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/4/2006	
*/	
void GIa_Histo_clear( GIa_HISTOGRAM *histo )	{
	if( histo->entrys!= GIP_NULL )	{
		GIP_free( histo->entrys );		histo->entrys=GIP_NULL;
	}
}

/*
	该算法总是收敛到直方图的波谷？
	因此适用于双峰直方图，即背景和目标区分明显

	效率 ****+
	健壮 ***+
	稳定 ****

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/5/2006	
*/	
void GIP_GlobalThresh_2( GIP_OBJECT *hObj,double *threshold,int flag )	{
	int i,M=hObj->M,N=hObj->N,MN=M*N,n_1,n_2;
	GIP_FLOATING *data = hObj->data;
	double u_1,u_2,a,thresh,th_1,epsi;

	thresh = 0.0;		//背景与对象的面积相近时，可取thresh为平均值 
	for( i = 0; i < MN; i++ )	thresh+= data[i];
	thresh /= MN;
	epsi = fabs(thresh*0.01);
	while( 1 )	{
		u_1=0.0;		u_2=0.0;
		n_1=0;			n_2=0;
		for( i = 0; i < MN; i++ )	{
			a = data[i];
//			a = GIP_LAPLAS_( i,M_p,N_p,N_p,hObj->data,0x0 );
			if( a>thresh )
			{	u_1 += a;		n_1++;	}
			else
			{	u_2 += a;		n_2++;	}
		}
		th_1 = (u_1/n_1+u_2/n_2)/2.0;
		if( fabs(th_1-thresh)<=epsi )
			break;
		thresh = th_1;
	}

	*threshold= thresh;
}

/*
	基于阈值迭代总能收敛到(直方图)波谷的假设

	效率 ***+
	健壮 -			仅针对iris寻找pupil的情况，需要完善各种情况
	稳定 *			基于阈值迭代总能收敛到(直方图)波谷的假设

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/9/2006	
*/	
int GIP_GlobalThresh_x( GIP_OBJECT *hObj,double I_0,double I_1,int *nTh,double *arrThresh,int flag )	{
	int i,M=hObj->M,N=hObj->N,MN=M*N,n_1,n_2,nz,ret=0,nLevel=0;
	GIP_FLOATING *data = hObj->data;
	double u_1,u_2,a,thresh,th_1,epsi;

	while( 1 )	{
		thresh = 0.0;	nz=0;	//背景与对象的面积相近时，可取thresh为平均值 
		for( i = 0; i < MN; i++ )	{
			a = data[i];		
			if( a < I_0 || a > I_1 )		continue;
			thresh+= data[i];		nz++;
		}
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

	*nTh=nLevel;
	return ret;
}


/*
	get basic attributes

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/9/2006	
*/	
void GIa_Histo_attrib( GIa_HISTOGRAM *histo,double *attrib,int flag )	{
	int i,I_0=histo->I_0,I_1=histo->I_1,n;
	double a,s,m;
	GIa_Histo_Entry *hEntry;

	a = 0.0;	n=0;
	for( i = I_0; i <= I_1; i++ )	{
		hEntry=histo->entrys+i;
		a += hEntry->level*hEntry->nz;
		n += hEntry->nz;
	}
	m = attrib[GIA_ATTRIB_MEAN] = a/n;
	a= 0.0;
	for(  i = I_0; i <= I_1; i++ )	{
		hEntry=histo->entrys+i;
		s = hEntry->level-m;
		a += s*s*hEntry->nz;;
	}
	attrib[GIA_ATTRIB_DEVIA] = sqrt(a/n);

	attrib[GIA_ATTRIB_CONTRAST] = histo->I_1-histo->I_0;
}

/*
	Circular Hough Transform

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/9/2008	
*/
double GIa_CHT( GIP_OBJECT *hObj,int nz,int* arrPos,double *R,double *c_r,double *c_c,int flag )	{
	int i,r,c,pos,M=hObj->M,N=hObj->N,*hough,hough_limit,h_max;
	int r_0,r_1,c_0,c_1,hc_0,hc_1,hr_0,hr_1,hR_1,hR_0,x_R,x_c,x_r,h_r,h_c,h_R;
	double *data=hObj->data,a,err=0.0;
	
	ASSERT( M>2 && N>2 && nz>8 );

	r_0=M*2;	r_1=-1;
	c_0=N*2;	c_1=-1;
	for( i = 0; i < nz; i++ )	{		
		pos = arrPos[i];
		r=GIP_POS2R(pos,N);		c=GIP_POS2C(pos,N);
		r_0=MIN( r_0,r);		r_1=MAX( r_1,r);		
		c_0=MIN( c_0,c);		c_1=MAX( c_1,c);		
	}
//	ASSERT( r_0>0 && r_1<M-1 && c_0>0 && c_1<N-1 );
/*	if( BIT_TEST( flag,GIa_CHT_CENTER_FIX_) )	{
		hR_1=MIN(r_1-r_0,c_1-c_0)/2+1;			hR_0=hR_1/2-1;
		hr_0 = *c_r;							hr_1 = *c_r+1;
		hc_0 = *c_c;							hc_1 = *c_c+1;
	}else*/	{
		hR_1=MIN(r_1-r_0,c_1-c_0)/2+1;			hR_0=hR_1/2-1;
		hr_0 = (r_1+r_0)/2-hR_0;				hr_1 = (r_1+r_0)/2+hR_0;
		hc_0 = (c_1+c_0)/2-hR_0;				hc_1 = (c_1+c_0)/2+hR_0;
	}

	x_R=hR_1-hR_0+1;		x_r=(hr_1-hr_0+1);			x_c=(hc_1-hc_0+1);
	hough_limit=x_R*x_r*x_c;
	hough=GIP_alloc( sizeof(int)*hough_limit );
	memset( hough,0x0,sizeof(int)*hough_limit );

	for( i = 0; i < nz; i++ )	{		
		pos = arrPos[i];
		r=GIP_POS2R(pos,N);		c=GIP_POS2C(pos,N);
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
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/12/2008	
*/
int _CoVar_Win_( int r,int c,int ldu,GIP_FLOATING *u,int W_M,int W_N,double *A,double *x,double *y,int flag )	{
	int i,s,t,nz=0,pos;
	double x_mean,y_mean;

	for( s = r;s < r+W_M; s++ )	{
	for( t = c;t < c+W_N; t++ )	{
		pos=GIP_RC2POS( s,t,ldu ); 
		if( u[pos]!=0.0)	continue;		
		x[nz]=s;			y[nz]=t;		//0.0代表边界上的点
		nz++;
	}
	}
	if( nz==0 )
		return nz;

	x_mean=0.0;		y_mean=0.0;
	for( i = 0; i < nz; i++ )	{
		x_mean+=x[i];		y_mean+=y[i];	
	}
	x_mean/=nz;				y_mean/=nz;
	A[0]=A[1]=A[2]=A[3]=0.0;
	for( i = 0; i < nz; i++ )	{
		A[0] += (x[i]-x_mean)*(x[i]-x_mean);
		A[3] += (y[i]-y_mean)*(y[i]-y_mean);
		A[1] += (x[i]-x_mean)*(y[i]-y_mean);
	}
	A[0]/=nz;		A[1]/=nz;		A[3]/=nz;
	A[2]=A[1];

	return nz;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/13/2008	
*/
int _CoVar_Edge_( int E_0,int E_1,int *E_pos,int ldu,int W_size,double *A,double *x,double *y,int flag )	{
	int i,k,s,t,nz=0,pos;
	double x_mean,y_mean,a,b;

	nz = E_1-E_0;
	if( nz==0 || nz<W_size )
		return nz;

	for( i = 0;	i < nz; i++ )	{
		pos=E_pos[i+E_0]; 
		x[i]=GIP_POS2R(pos,ldu);			y[i]=GIP_POS2C(pos,ldu);		//0.0代表边界上的点
	}
	
	for( k = 0; k<nz-W_size; k+=W_size )	{
		x_mean=0.0;		y_mean=0.0;
		for( i = k; i < k+W_size; i++ )	{
			x_mean+=x[i];		y_mean+=y[i];	
		}
		x_mean/=W_size;				y_mean/=W_size;
		A[0]=A[1]=A[2]=A[3]=0.0;
		if( 1 )	{
			for( i = k; i < k+W_size; i++ )	{
				A[0] += (x[i]-x_mean)*(x[i]-x_mean);
				A[3] += (y[i]-y_mean)*(y[i]-y_mean);
				A[1] += (x[i]-x_mean)*(y[i]-y_mean);
			}
			A[0]/=W_size;		A[1]/=W_size;		A[3]/=W_size;
		}else	{
			for( i = k; i < k+W_size; i++ )	{
				A[0] +=x[i]*x[i];		A[1] +=x[i]*y[i];		A[3] +=y[i]*y[i];
			}
			A[0]=A[0]/W_size-x_mean*x_mean;		
			A[1]=A[1]/W_size-x_mean*y_mean;		
			A[3]=A[3]/W_size-y_mean*y_mean;
		}
		A[2]=A[1];
		a = sqrt( (A[3]-A[0])*(A[3]-A[0])+4.0*A[1]*A[1] );
		b = (A[3]+A[0]);
	}

	return nz;
}


/*
	Elliptical Features Extraction Using EigenValues of Covariance Matrices
	buffer[W_rank*W_rank*3]

	效率 **
		需对每个Edge点执行窗口操作		
	健壮 *
		由于离散为像素，1）园的轮廓并不是旋转不变的 2）在小窗口里近似为直线
		需要半径已知,以确定合适的窗口尺寸。
	稳定 *

	注意:
		1 hEdge已是EdgeMap

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/12/2008	
*/
void GIa_Circle_Detect_Eigen( GIP_OBJECT *hEdge,int type,double *D_buffer,int *I_buffer,int flag )	{
	int i,j,r,c,s,t,pos,M=hEdge->M,N=hEdge->N,MN=M*N,W_rank=60,nz;
	int *E_ptr,*E_pos,nEdge,len;
	GIP_FLOATING *data=hEdge->data;
	double *A,*x,*y,lend_1,lend_2,a,b;
	GIP_OBJECT *hNew;

	hNew=GIP_OBJ_create( M,N,GIP_NULL,0x0 );
	for( i = 0;	i < MN; i++ )	hNew->data[i]=1.0;
//	GIP_Circle( hEdge,M/2,N/2,30.0,GIP_IMAG_EDGE );

	E_ptr=I_buffer;		E_pos=E_ptr+MN;	
	nEdge=GIP_Edge2List( hEdge,E_ptr,E_pos,0x0,0x0 );
/*	for( i = 0; i < nEdge; i++ )	{
		for( j = E_ptr[i]; j < E_ptr[i+1]; j++ )	{
			hNew->data[E_pos[j]] = 0.0;
		}
	}*/
	A=D_buffer;		
	x=D_buffer+4;			y=D_buffer+MN;
	for( i = 0; i < nEdge; i++ )	{
		len = E_ptr[i+1]-E_ptr[i];
		_CoVar_Edge_( E_ptr[i],E_ptr[i+1],E_pos,N,W_rank,A,x,y,0x0 );
		a = sqrt( (A[3]-A[0])*(A[3]-A[0])+4.0*A[1]*A[1] );
		b = (A[3]+A[0]);
		if( a<b*0.1 || len>300 )	{		//its a circle
			for( j = E_ptr[i]; j < E_ptr[i+1]; j++ )	{
				hNew->data[E_pos[j]] = 0.0;
			}
		}
	}
/*	for( r = 0; r+W_rank < M; r++ )	{		
	for( c = 0; c+W_rank < N; c++ )	{
		pos=GIP_RC2POS( r,c,N ); 
		nz = _CoVar_Win_( r,c,N,data,W_rank,W_rank,A,x,y,0x0 );
		if( nz==0 )
			continue;
		a = sqrt( (A[3]-A[0])*(A[3]-A[0])+4.0*A[1]*A[1] );
		b = (A[3]+A[0]);
		if( a<b*0.01 )	{		//its a circle
			for( s = r;s < r+W_rank; s++ )	{
			for( t = c;t < c+W_rank; t++ )	{
				pos=GIP_RC2POS( s,t,N ); 
				hNew->data[pos]=data[pos];
			}
			}
		}
	}
	}	
	GIP_save_IplImage_d( hNew->M,hNew->N,".\\trace\\eigen.jpg",hNew->data,hNew->N,0x0 );
*/	
	GIP_OBJ_clear( hNew );
}

/*
	Circle Detection Using Gradient Pair Vectors
	buffer[W_rank*W_rank*3]

	效率 ***
		与CHoughTransform相比，半径能确定		
	健壮 *
		1 基于梯度，即使采用SOBEL等算子，误差仍很大（尤其对于较大半径的园）
		2 基于circular object比背景暗的假设
	稳定 *

	ISSUE-NEEDED:
		1 t_pair需要动态数组替代，绝大多数情况下 nzPair<MN
		2 int *R_ptr,int*R_center,*R_sum

	注意：

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/14/2008	
*/
void GIa_Circle_Detect_GPV( GIP_OBJECT *hObj,int th_1,int type,double *D_buffer,int *I_buffer,int flag )	{
	int i,j,k,r,s,t,c,r1,c1,r2,c2,pos,M=hObj->M,N=hObj->N,MN=M*N,ldu=N,id;
	int grid=360,window=10,R,nz_1,nz_2,R_limit=MIN(M,N)*0.4,*R_ptr,*R_center,*R_sum,*R_cur;
	int nzPair,*t_pair,*gv_ptr,*gv_ind,no;
	GIP_FLOATING *data=hObj->data;
	double uc,ur,g,angle,xita=2*PI/grid,g_thrsh=5.0;
	GIP_OBJECT *hNew;
	PACT_order_list lGV;

	hNew=GIP_OBJ_create( M,N,GIP_NULL,GIP_OBJ_2ZERO );	
	PACT_ol_init( &lGV,grid+1,MN,0x0 );
//	GIP_Circle( hObj,M/2,N/2,100.0,GIP_IMAG_EDGE );

	for( r = 1; r < M-1; r++ )	{		
	for( c = 1; c < N-1; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
//		if( data[pos]>0.0 )		continue;
//		if( data[pos]>th_1 || data[pos-1]>th_1 || data[pos+1]>th_1 
//			|| data[pos-ldu]>th_1 || data[pos+ldu]>th_1 )		
//			continue;
		if( data[pos]>data[pos+ldu] && data[pos]>data[pos-ldu] &&
			data[pos]>data[pos+1] && data[pos]>data[pos-1] )		
			continue;
//		uc = (data[pos+1]-data[pos-1])/2.0;
//		ur = (data[pos+ldu]-data[pos-ldu])/2.0;
		ur = GIP_Mask_3_( pos,M,N,ldu,data,G_MASK_SOBEL_R,0x0 );		
		uc = GIP_Mask_3_( pos,M,N,ldu,data,G_MASK_SOBEL_C,0x0 );	
		if( ur==0.0 && uc==0.0 )		continue;
		g = sqrt( uc*uc+ur*ur );
		if( g < g_thrsh )	continue;
//		hNew->data[pos]=1.0;
		angle = atan2( ur,uc );
		if( angle<0.0 )		angle=2*PI+angle;
		id=(int)(angle/(xita)+0.5);			ASSERT( id>=0 && id <=grid );
		PACT_ol_insert( &lGV,pos,id );
	}
	};
//	GIP_save_IplImage_d( hNew->M,hNew->N,".\\trace\\GPV.jpg",hNew->data,hNew->N,0x0 );
	
	nzPair=0;
	t_pair=GIP_alloc( sizeof(int)*MN );
	R_ptr=GIP_alloc( sizeof(int)*(R_limit+1)*2 );			R_cur=R_ptr+R_limit+1;
	for( i = 0; i < R_limit+1; i++ )	R_ptr[i]=0;
	gv_ptr=I_buffer;			gv_ind=I_buffer+lGV.compact+1;
	PACT_ol_2ccs( &lGV,gv_ptr,gv_ind,0x0 );
	for( s = 0; s < grid/2; s++ )	{
		nz_1 = gv_ptr[s+1]-gv_ptr[s];
		if( nz_1==0 )		continue;
		for( no = -window/2; no < window/2; no++ )	{
			t=(s+no+grid/2)%grid;
			nz_2 = gv_ptr[t+1]-gv_ptr[t];
			if( nz_2==0 )		continue;
			for( i = gv_ptr[s]; i < gv_ptr[s+1]; i++ )	{			//matching
				r1=GIP_POS2R(gv_ind[i],ldu);			c1=GIP_POS2C(gv_ind[i],ldu);
				for( j = gv_ptr[t]; j < gv_ptr[t+1]; j++ )	{			//matching
					r2=GIP_POS2R(gv_ind[j],ldu);		c2=GIP_POS2C(gv_ind[j],ldu);
					if( fabs(data[gv_ind[i]]-data[gv_ind[j]])>0 )		
						continue;
					angle = atan2( r1-r2,c1-c2);
					if( angle<0.0 )		angle=2*PI+angle;
					id=(int)(angle/(xita)+0.5);
					if( abs(id-s)>window/2 && abs(id-s-grid/2)>window/2 )		
						continue;
					r=(int)((r1+r2)/2.0+0.5),		c=(int)((c1+c2)/2.0+0.5);
					R=(int)(sqrt( (r1-r2)*(r1-r2)+(c1-c2)*(c1-c2) )/2+0.5);
					if( R>=R_limit || R<=1 )		
						continue;
					pos=GIP_RC2POS( r,c,N );
					R_ptr[R] ++;			
					hNew->data[pos] += 1.0;///R;
					t_pair[2*nzPair]=gv_ind[i];			t_pair[2*nzPair+1]=gv_ind[j];
					nzPair++;
				}
			}
		}
	}
//	GIP_save_IplImage_d( hNew->M,hNew->N,".\\trace\\GPV.jpg",hNew->data,hNew->N,0x0 );
	ASSERT(nzPair<MN );
	R_center=GIP_alloc( sizeof(int)*(nzPair*2) );		
	R_sum=R_center+nzPair;			R_sum[0]=0;
	R_ptr[R_limit]=nzPair;
	for( i = R_limit-1; i >= 0; i-- )	{
		R_ptr[i] = R_ptr[i+1]-R_ptr[i];			R_cur[i]=R_ptr[i];
	}
	ASSERT( R_ptr[0]==0 );
	for( i = 0; i < nzPair; i++ )	{
		r1=GIP_POS2R(t_pair[2*i],ldu);				c1=GIP_POS2C(t_pair[2*i],ldu);
		r2=GIP_POS2R(t_pair[2*i+1],ldu);			c2=GIP_POS2C(t_pair[2*i+1],ldu);
		r=(int)((r1+r2)/2.0+0.5),		c=(int)((c1+c2)/2.0+0.5);
		R=(int)( sqrt( (r1-r2)*(r1-r2)+(c1-c2)*(c1-c2) )/2+0.5 );
		pos=GIP_RC2POS( r,c,N );
		no=-1;
		for( j = R_ptr[R]; j < R_cur[R]; j++ )	{
			if( R_center[j]==pos )
			{	no=j;		break;	}
		}
		if( no==-1 )	{
			R_center[R_cur[R]]=pos;
			R_sum[R_cur[R]]=1;
			R_cur[R]++;			
		}else
			R_sum[no]++;
	}
	r1=0;		no=0;
	for( i = 0; i < R_limit; i++ )	{
		ASSERT( R_cur[i]<=R_ptr[i+1] );
	//	if( i>=73  )	continue;
		for( j = R_ptr[i]; j < R_cur[i]; j++ )	{
			if( R_sum[j]>R_sum[no] )
			{	R=i;		no = j;			}
		}
	}
	r1=GIP_POS2R(R_center[no],ldu);				c1=GIP_POS2C(R_center[no],ldu);
	GIP_Circle( hObj,r1,c1,R,GIP_IMAG_ORIGIN );	//draw circle at R_center[j] with Radius i
//	GIP_save_IplImage_d( hObj->M,hObj->N,".\\trace\\GPV_1.jpg",hObj->data,hObj->N,0x0 );

	GIP_free( R_ptr );			GIP_free( R_center );			GIP_free( t_pair );
	PACT_ol_clear( &lGV );
	GIP_OBJ_clear( hNew );
}
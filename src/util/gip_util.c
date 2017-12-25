#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "GIP_util.h"

double GIP_INFO[GIP_INFO_ITEM];
double GIP_CONTROL[GIP_CONTROL_ITEM];

GIP_HEAP_INFO heap_info;

/*
	v0.1	cys
		9/2/2004
*/
void GIP_init_heapinfo( )	{
	heap_info.alloc = 0;
	heap_info.nAlloc = 0;
	heap_info.nFree = 0;
	heap_info.old = 0;
	heap_info.peak = 0.0;

//	DECLARE_CS( );
}

/*
	v0.1	cys
		9/2/2004
*/
void GIP_exit_heapinfo( )	{
	int nStay=heap_info.nAlloc-heap_info.nFree;
	if( nStay!=0 )	{
		G_PRINTF( _T("%d block unfreed, unfreed memory is %d"),nStay,heap_info.alloc );
	}
	ASSERT( nStay == 0  );
	ASSERT( heap_info.alloc == 0 );

//	DELETE_CS( );
}

/*
	v0.1	cys
		7/1/2005
*/
double GIP_peak_heap( )	{	
	return heap_info.peak;
}

/*
	return the heap difference between two check point
	v0.1	cys
		9/3/2004
*/
size_t HEAP_diff( )	{
	size_t diff = heap_info.alloc-heap_info.old;
	heap_info.old = heap_info.alloc;
	return diff;
}
/*
	v0.1	cys
		9/2/2004
*/
void* GIP_alloc( size_t size )	{
	void* ptr = GIP_NULL;
	size_t m_alloc;				//6/8/2007

//	ENTER_CS( );
	ptr = malloc( size );
	if( ptr==0x0 )
		goto END;
//	ptr = VirtualAlloc( NULL,size,MEM_RESERVE|MEM_COMMIT,PAGE_READWRITE );
	m_alloc = _msize( ptr );
	heap_info.alloc += m_alloc;
	ASSERT( ptr!=GIP_NULL && m_alloc>=size );
//	heap_info.temp[heap_info.nAlloc] = ptr;
	heap_info.nAlloc++;
	heap_info.peak = MAX( heap_info.peak,heap_info.alloc );
//	LEAVE_CS( );
END:
	return ptr;
}

/*
	v0.1	cys
		9/2/2004
*/
void* GIP_calloc(  size_t num,size_t size )	{
	void* ptr = GIP_NULL;
//	ENTER_CS( );
	ptr = calloc( num,size );
	ASSERT( ptr!=GIP_NULL );
	heap_info.alloc += _msize( ptr );
//	heap_info.temp[heap_info.nAlloc] = ptr;
	heap_info.peak = MAX( heap_info.peak,heap_info.alloc );
	heap_info.nAlloc++;
//	LEAVE_CS( );
	return ptr;
}
/*
	v0.1	cys
		6/13/2005
*/
void* GIP_realloc( void* old,size_t size )	{
	void* ptr = GIP_NULL;
//	ENTER_CS( );
	heap_info.alloc += old==GIP_NULL ? size:size-_msize( old );
	ptr = realloc( old,size );
	ASSERT( ptr!=GIP_NULL );
	if( old==GIP_NULL )heap_info.nAlloc++;
//	heap_info.temp[heap_info.nAlloc] = ptr;
//	LEAVE_CS( );
	return ptr;

}

/*
	v0.1	cys
		9/2/2004
*/
void GIP_free( void* memblock  )	{
	size_t size;

//	ENTER_CS( );
	size = _msize( memblock );
	heap_info.alloc -= size;
	free( memblock );
/*	VirtualFree( memblock,0,MEM_DECOMMIT );	*/
	heap_info.nFree++;
//	LEAVE_CS( );
}

/*
	Int Vector输出

	v0.1	cys
		11/8/2004 void I_ouput
*/
void I_ouput( int len,int* val,const char* sPath )	{
	int i;
	FILE* fp=_tfopen( sPath,"w" );

	G_PRINTF( "{I_ouput} is called to print file <<<%s>>> \r\n",sPath  );

	fprintf( fp,"%d\r\n",len );
	for( i = 0; i < len; i++ )	{
		fprintf( fp,"%d\t%d\r\n",i+1,val[i] );
	}
	fclose( fp );
}

/*
	Read Int Vector 

	v0.1	cys
		5/23/2005 void I_input
*/
void I_input( int len,int* val,const char* sPath )	{
	int i,no;

	FILE* fp=_tfopen( sPath,"r" );
	fscanf( fp,"%d",&no );
	for( i = 0; i < len; i++ )	{
		fscanf( fp,"%d\t%d\r\n",&no,val+i );
	}
	fclose( fp );
}

/*
	v0.1	cys
		5/25/2009 
*/
void D_input( int len,double* val,const TCHAR* sPath )	{
	int i,len_0,no;
	double a;
	FILE* fp=_tfopen( sPath,"r" );

	_ftscanf( fp,_T("%d\r\n"),&len_0 );
	ASSERT( len_0>=len );
	for( i = 0; i < len; i++ )	{
		_ftscanf( fp,_T("%d\t%lf\r\n"),&no,&(a) );
		val[i] = a;
	}
	fclose( fp );
}

/*
	v0.1	cys
		5/23/2004 
*/
void D_ouput( int len,double* val,const TCHAR* sPath )	{
	int i;
	FILE* fp=_tfopen( sPath,"w" );

	G_PRINTF( _T("{D_ouput} is called to print file <<<%s>>> \r\n"),sPath  );

	_ftprintf( fp,_T("%d\r\n"),len );
	for( i = 0; i < len; i++ )	{
		_ftprintf( fp,_T("%d\t%.15g\r\n"),i+1,val[i] );
	}
	fclose( fp );
}

/*
	输出将val对应的distribution

	v0.1	cys
		8/6/2004 
*/
void distri_ouput( int len,double *w,double* val,const TCHAR* sPath,int flag )	{
	int i;
	double a;
	FILE* fp=_tfopen( sPath,"w" );

	G_PRINTF( _T("{distri_ouput} is called to print file <<<%s>>> \r\n"),sPath  );

	_ftprintf( fp,_T("%d\r\n"),len );
	a = 0.0;
	for( i = 0; i < len; i++ )	{
//		a = (i*1.0)/(len-1);
		a += w[i];
		_ftprintf( fp,_T("%.6g\t%.15g\r\n"),a,val[i] );
	}
	fclose( fp );
}

typedef struct {
	int* m_stkInt,maxInt,topInt,maxDouble,topDouble;
	double* m_stkDouble;
}GIP_BUFFER;
static GIP_BUFFER grus_Buffer;

/*
	v0.1	cys
		1/31/2005
*/
void GIP_BUFFER_Init( int maxI,int maxD )	{
	grus_Buffer.maxInt = maxI;				grus_Buffer.topInt = 0;
	grus_Buffer.maxDouble = maxD;			grus_Buffer.topDouble = 0;
	grus_Buffer.m_stkInt = GIP_alloc( sizeof(int)*maxI );
	grus_Buffer.m_stkDouble = GIP_alloc( sizeof(double)*maxD );
}

/*
	v0.1	cys
		1/31/2005
*/
void GIP_BUFFER_Free( )	{
	GIP_free( grus_Buffer.m_stkInt );		grus_Buffer.m_stkInt=GIP_NULL;
	GIP_free( grus_Buffer.m_stkDouble );	grus_Buffer.m_stkDouble=GIP_NULL;
}

/*
	v0.1	cys
		1/31/2005
*/
int* GIP_BI_alloc( int len )	{
	int*m_stkInt = grus_Buffer.m_stkInt;
	ASSERT( grus_Buffer.topInt+len <= grus_Buffer.maxInt );
	grus_Buffer.topInt += len;
	return m_stkInt+grus_Buffer.topInt-len;
}

/*
	v0.1	cys
		1/31/2005
*/
double* GIP_BD_alloc( int len )	{
	double*m_stkDouble = grus_Buffer.m_stkDouble;
	ASSERT( grus_Buffer.maxDouble+len <= grus_Buffer.maxDouble );
	grus_Buffer.maxDouble += len;
	return m_stkDouble+grus_Buffer.maxDouble-len;
}

/*
	v0.1	cys
		1/31/2005
*/
void GIP_BI_free( int len )	{
	ASSERT( grus_Buffer.topInt-len >= 0  );
	grus_Buffer.topInt -= len;
}

/*
	v0.1	cys
		1/31/2005
*/
void GIP_BD_free( int len )	{
	ASSERT( grus_Buffer.topDouble-len >= 0  );
	grus_Buffer.topDouble -= len	;
}

int grus_stp,*grus_signal,grus_signal_len;

/*
	v0.1	cys
		3/8/2005
*/
void GIP_STAMP_Init( int order )	{
	grus_stp = 0;
	grus_signal_len = order;
	grus_signal=GIP_alloc(sizeof(int)*grus_signal_len);
//	cblas_iset( order,grus_stp,grus_signal,1 );
	memset( grus_signal,0x0,sizeof(int)*grus_signal_len );
}

/*
	v0.1	cys
		3/8/2005
*/
void GIP_STAMP_Clear( )	{
	GIP_free(grus_signal);
}

/*
	v0.1	cys
		3/8/2005
*/
void GIP_STAMP_Next( int need )	{
	if( grus_stp+need > GIP_INT32_MAX )	{
		grus_stp = 0;
//		cblas_iset( grus_signal_len,grus_stp,grus_signal,1 );
		memset( grus_signal,0x0,sizeof(int)*grus_signal_len );
	}
	grus_stp+=need;
}

/*
	v0.1	cys
		3/8/2005
*/
void GIP_STAMP_Need( int need )	{
	if( grus_stp+need > GIP_INT32_MAX )	{
		grus_stp = 0;
//		cblas_iset( grus_signal_len,grus_stp,grus_signal,1 );
		memset( grus_signal,0x0,sizeof(int)*grus_signal_len );
	}
}

/*
	v0.1	cys
		3/8/2005

_inline void GIP_STAMP_Set( int ele ){
	grus_signal[ele] = grus_stp;
}
*/
/*
	v0.1	cys
		3/8/2005

_inline int GIP_STAMP_Test( int ele ){
	return grus_signal[ele] == grus_stp;
}
*/
/*
	v0.1	cys
		9/9/2004 void NO_VEC( )	
*/
int NO_VEC( int r,int len,int * vec )	{
	int i,no=-1;
	for( i = 0; i < len; i++ )	{
		if( vec[i] == r )
		{	no = i;		break;	}
	}
	return no;
}

/*
	v0.1	cys
		2/15/2004	VERIFY_GRAPH_2	
*/
int VERIFY_GRAPH_2( int order,int* ptr,int* ind )	{
	int i,j,v,v_pre;
	for( i = 0; i < order; i++ )	{
		v_pre = -1;
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			v = ind[j];
			ASSERT( v >= 0 && v < order );
//			ASSERT( v > v_pre );
			v_pre = v;
		}
	}

	return 0;
}

/*
	v0.1	cys
		3/4/2005	VERIFY_VEC_SAME	
*/
void VERIFY_VEC_SAME( int len_1,int* vec_1,int len_2,int* vec_2 )	{
	int i;
	ASSERT( len_1 == len_2 );
	for( i = 0; i < len_1; i++ )	{
		ASSERT( vec_1[i] == vec_2[i] );
	}
}

/*
	matrix is ordered!!!

	v0.1	cys
		9/21/2004
	v0.2	cys
		1/6/2005 Col_Merge_MA( int head,int len,int* ind,int* ma )
*/
void Col_Merge_MA( int head,int len,int* ind,int* ma )	{
	int i,r,cur,pre;
	for( i = 0; i <len; i++ )	{
		r = ind[i];
		ASSERT( r >= head );
		cur = head;
		while( r > cur )	{
			pre = cur;
			cur = ma[cur];
		}
		if( r != cur )	{	//insert r before m
			ma[pre] = r;
			ma[r] = cur;
		}
		head = r;	//ind is ordered,so inc head
	}
}

/*
	打开文件，写入必要的初始信息 

	v0.1	cys
		11/24/2007
*/
void grus_dump_init( void **hDump,char *sTitle )	{
	double *info=GIP_INFO,*control=GIP_CONTROL;
	int dump=(int)(control[GIP_DUMP]);
	FILE *fp=0x0;
	char sPath[GIP_MAX_LINE],*strR="\r\n";
	
//	if( !BIT_IS( dump,GIP_DUMP_ALL) )
//		goto END;
	_stprintf( sPath,"%s.info",GIP_APP_NAME );

	fp=_tfopen( sPath,"w" );
	fprintf( fp,"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\" ?>%s"
				"<GSP_TEST_INFO>%s",strR,strR );
	fprintf( fp,  
		"<DECLARE>\r\n\t<APP>%s</APP>%s"
		"\t<COPYRIGHT>Copyright (c) 2005-present by YingShi Chen.	All Rights Reserved.</COPYRIGHT>%s"
		"\t<MAIL>gsp@grusoft.com</MAIL>%s</DECLARE>%s%s",sTitle,strR,strR,strR,strR,strR
	);
//END:
	ASSERT( fp != 0x0 );
	*hDump = fp;

	return;
}

/*
	打开文件，写入必要的初始信息 

	v0.1	cys
		11/24/2007
*/
void grus_dump_finish( void *hDump )	{
	FILE *fp=(FILE *)hDump;
	char *strR="\r\n";
	
	if( fp != 0x0 )	{
		fprintf( fp,"%s</GSP_TEST_INFO>%s",strR,strR );
		fclose( fp );
	}
}

/*
	v0.1	cys
		12/23/2007
*/
void GIP_Get_Control_Info( double **control,double **info,int flag )	{
	if( control!=0x0 )
		*control=GIP_CONTROL;
	if( info!=0x0 )
		*info=GIP_INFO;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	

void GIP_IMAGE_init( int M,int N,double *val,GIP_IMAGE *hImag )	{
	hImag->M = M;		hImag->N = N;		
	if( hImag->val != GIP_NULL )
		GIP_free(hImag->val);		
	hImag->val=GIP_alloc( sizeof(double)*M*N );
	memcpy( hImag->val,val,sizeof(double)*M*N );
}
*/
/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/26/2008	

void GIP_IMAGE_clear( GIP_IMAGE *hImag )	{
	if( hImag->val != GIP_NULL )
	{	GIP_free(hImag->val);		hImag->val=GIP_NULL;	}
}
*/

/*
	[ISSUE-NEEDED]
		file_name需要转化为TCHAR

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/9/2008	
*/
void GIP_Error( int status, const TCHAR* head,const TCHAR* err_msg, const TCHAR* file_name, int line )		{
/*    if( code == GIP_OK )
        cvSetErrStatus( code );
    else
    {
        CvContext* context = icvGetContext();

        if( code != CV_StsBackTrace && code != CV_StsAutoTrace )
        {
            char* message = context->err_msg;
            context->err_code = code;

            _tcscpy( message, err_msg );
            context->err_ctx.file = file_name;
            context->err_ctx.line = line;
        }

        if( context->err_mode != CV_ErrModeSilent )
        {
            int terminate = context->error_callback( code, func_name, err_msg,
                                                    file_name, line, context->userdata );
            if( terminate )
            {
#if !defined WIN32 && !defined WIN64
                assert(0); // for post-mortem analysis with GDB
#endif
                exit(-abs(terminate));
            }
        }
    }	
	;*/
//	G_PRINTF( _T("%s(%d):\tERROR(%d):%s."),file_name,line,status,err_msg );
	G_PRINTF( _T("(%d):\tERROR(%d):%s."),line,status,err_msg );
}

/*
	return |Va-Vb|	(|| is Euclidean norm)

	v0.1	cys
		5/23/2008
*/
double  GIP_DOFF2( int dim,double *Va,double *Vb )	{
	int i;
	double off=0.0,a;
	for( i = 0; i < dim; i++ )	{
		a = Va[i]-Vb[i];
		off += a*a;
	}
	off = sqrt(off);

	return off;
}

/*
	return position of max |Va-Vb|	

	v0.1	cys
		5/27/2008
*/
int  GIP_MAXOFF_POS( int dim,double *Va,double *Vb,int flag )	{
	int i,pos=0;
	double off=0.0,max_off=0.0;
	for( i = 0; i < dim; i++ )	{
		off = fabs(Va[i]-Vb[i]);
		if( off > max_off )	{
			max_off=off;		pos=i;
		}
	}
	off = sqrt(off);

	return pos;
}

/*
	user interactive control

	v0.1	cys
		5/29/2008
*/
void GIP_interactive( int M,int N,double *u,int ldu,int no,int *t_step,int flag )	{
/*	int nKey;
	char sTitle[80];

	_stprintf( sTitle,".\\trace\\u_%d.bmp",no );		
	GIP_save_IplImage_d( M,N,sTitle,u,ldu,0x0 );
//			_stprintf( sTitle,"u_%d.bmp",tvm_step );		
//			GIP_show_IplImage_d( M,N,sTitle,(tvm.u+tvm.shift),tvm.ldu,GIP_OPENCV_CREATEWND );
//			_stprintf( sTitle,".\\trace\\interface_%d.bmp",step );		
	nKey = cvWaitKey( );
	switch( nKey )	{
		case 27:
			G_PRINTF( "TV MODEL:	fixed point iteration break out\r\n"	);
			GIP_EXIT;
		case '+':
			*t_step *= 10;
			break;
		case '-':
			*t_step /= 10;
			break;
		default:
			break;
	}
//			GIP_save_IplImage_d( M,N,sTitle,(temp+tvm.shift),tvm.ldu,0x0 );
GIP_exit:
	;*/
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/14/2008	
*/
void Equation_output( int M,int N,int dim,int *ptr,int *ind,double *val,double *rhs,const char* sPath )	{
	int MN=M*N,i,j,nnz,r,c,*crs_ptr,*crs_ind,*temp;
	double *crs_val;
	FILE *fp=NULL;
	
	nnz = ptr[dim];
	crs_ptr=GIP_alloc( sizeof(int)*(dim+1+nnz) );		crs_ind=crs_ptr+dim+1;
	crs_val=GIP_alloc( sizeof(double)*nnz );
	temp=GIP_alloc( sizeof(int)*(dim+1) );
	memset( temp,0x0,sizeof(int)*(dim+1) );
	for( i = 0; i < nnz; i++ )	{
		r = ind[i];
		temp[r]++;
	}
	crs_ptr[dim]=nnz;
	for( i = dim-1;	i >= 0; i-- )	{
		crs_ptr[i] = crs_ptr[i+1]-temp[i];
		temp[i] = crs_ptr[i];
	}
	ASSERT( crs_ptr[0] == 0 );

	for( i = 0;	i < dim; i++ )	{
		for( j = ptr[i]; j < ptr[i+1]; j++ )	{
			r = ind[j];
			crs_ind[temp[r]] = i;		crs_val[temp[r]] = val[j];
			temp[r]++;
		}
	}
	for( i = 0;	i < dim; i++ )	ASSERT( temp[i]==crs_ptr[i+1] );

	G_PRINTF( "<<%s>> outputing...\r",sPath);
	fp = _tfopen( sPath,"w" );
	if( fp==NULL )	{
//		G_PRINTF( "<<%s>> outputing failed.\n",sPath);
		GIP_ERROR( "Equation_output",-1, "outputing failed" );
	}
	fprintf( fp,"no\tu\trhs\n" );
	for( i = 0;	i < dim; i++ )	{
		fprintf( fp,"%d\t%g\n",i+1,rhs[i] );
	}
	fprintf( fp,"************************\r\n",dim,nnz );
	fprintf( fp,"%d %d\n",dim,nnz );
	for( i = 0;	i < dim; i++ )	{
		for( j = crs_ptr[i]; j <crs_ptr[i+1]; j++ )	{
			c = crs_ind[j];
			fprintf( fp,"%d\t%d\t%lf\n",i+1,c+1,crs_val[j] );
		}
	}
	fclose( fp );
	G_PRINTF( "<<%s>> outputing finish.\n",sPath);
GIP_exit:		
	GIP_free(crs_ptr);			
	GIP_free(crs_val);
	GIP_free(temp);
}

/*
	modified form Mat_Trans.c in GSP 2.0
	修改:
		1 自动分配内存
		2 支持crs_data=NULL

	v0.1	cys
		8/28/2008	Initial Draft	
*/
void GIP_CRS2CCS( int nRow,int nCol,GIP_FLOATING *crs_data,int* crs_ptr,
			 int *crs_ind,GIP_FLOATING **data,int **colptr,int **rowind,int* temp )	{
	int i,col,row,nnz;

	nnz=crs_ptr[nRow];
	if( crs_data!=GIP_NULL )
		*data=GIP_alloc( sizeof(GIP_FLOATING)*nnz );
	*colptr=GIP_alloc( sizeof(GIP_FLOATING)*(nCol+1) );
	*rowind=GIP_alloc( sizeof(GIP_FLOATING)*nnz );

	memset( temp,0,sizeof(int)*nCol );

	for( i = 0; i < nnz; i++ )	
		temp[crs_ind[i]]++;
	(*colptr)[0] = 0;
	for( i = 1; i < nRow+1; i++ )	{	
		(*colptr)[i] = (*colptr)[i-1]+temp[i-1];
		temp[i-1] = (*colptr)[i-1];
	}
	for( row = 0; row < nRow; row++ )	{
		for( i = crs_ptr[row]; i < crs_ptr[row+1]; i++ )	{
			col = crs_ind[i];
			if( crs_data!=GIP_NULL )
				(*data)[temp[col]] = crs_data[i];
			(*rowind)[temp[col]] = row;
			temp[col]++;
		}
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/27/2008	
*/
int GIP_conmpact_index( int M,int N,int ld,GIP_FLOATING *val )	{
	int index=0,r,c,pos;
	GIP_FLOATING a;

	for( r = 0; r < M; r++ )	{		
		for( c = 0; c < N; c++ )	{
			if( r==0 || c==0 || r==M-1 || c==N-1 )	{
				index++;	continue;
			}
			pos = GIP_RC2POS( r,c,ld );
			a = val[pos];
			if( a == val[pos+1] && a == val[pos-1] && a == val[pos+ld] && a == val[pos-ld] )
				continue;
			index++;							
		}
	}	

	return index;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/27/2008	
*/
void GIP_OBJ_clear( GIP_OBJECT *hObj )	{
	if( hObj->data != GIP_NULL )	
	{	GIP_free( hObj->data );		hObj->data=GIP_NULL; }
}

/*
	RMSE:	residual-mean-squared-error
	SNR:	signal noise ratio

	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/10/2008	
*/
void GIP_RES_quality( int M,int N,int ldu,GIP_FLOATING *u,GIP_FLOATING *u_0,double *rmse,double *snr,double *v_1_,double *v_2_ )	{
	int r,c,pos;
	double a1,a2,u_avg=0.0,v_avg=0.0,off_2;

	u_avg=0.0,	v_avg=0.0;
	off_2=0.0;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
		u_avg += u[pos];
		v_avg += u_0[pos]-u[pos];
		off_2 += (u[pos]-u_0[pos])*(u[pos]-u_0[pos]);
	}
	}
	u_avg /= M*N;
	v_avg /= M*N;
	*rmse = sqrt(off_2/(M*N));

	a1=0.0;			
	a2=0.0;
	for( r = 0; r < M; r++ )	{		
	for( c = 0; c < N; c++ )	{
		pos = GIP_RC2POS( r,c,ldu );
		a1 += (u[pos]-u_avg)*(u[pos]-u_avg);
		a2 += (u_0[pos]-u[pos]-v_avg)*(u_0[pos]-u[pos]-v_avg);
	}
	}

	*snr = a2==0.0 ? 1.0e100 : sqrt(a1)/sqrt(a2);
	*v_1_ = v_avg;
	*v_2_ = *rmse;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		11/11/2008		
*/
void GIP_MODEL_info( char *sModel,GIP_MODEL model,int flag )	{
	switch( model )	{
	case GIP_MODEL_130:
		_stprintf( sModel,"%s","MODEL_130" );;
		break;
	case GIP_SEG_MUMSHAH_s:
		_stprintf( sModel,"%s","MUMSHAH_TV_Hs" );;
		break;
	case GIP_SEG_MUMSHAH_2:
		_stprintf( sModel,"%s","MUMSHAH_TV_H2" );;
		break;
	case GIP_MUMSHAH_1:
		_stprintf( sModel,"%s","MUMSHAH_TV_H1" );;
		break;
	case GIP_MUMSHAH_filter:
		_stprintf( sModel,"%s","MUMSHAH_TV_filter" );
		break;
/*	case GIP_DECOMP_Hs:
		_stprintf( sModel,"%s","DECOMP_TV_Hs" );;
		break;*/
	case GIP_MODEL_CGM:
		_stprintf( sModel,"%s","CGM" );;
		break;
	case GIP_SEG_CV:
		_stprintf( sModel,"%s","SEG_CV" );;
		break;
	default:
		_stprintf( sModel,"%s","!!!!" );
		break;
	}
	;
}

/*
	ISSUE-NEEDED:
		1 以后可用模板代替

	//{-1,1,-ldu,ldu}
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/5/2006	
*/	
int GIP_Neibors( int pos,int *neibors,int M,int N,int ldu,int flag )	{
	int r=GIP_POS2R(pos,ldu),c=GIP_POS2C(pos,ldu),nz=0;

	if( c>0 )	neibors[nz++]=pos-1;
	if( c<N-1 )	neibors[nz++]=pos+1;
	if( r>0 )	neibors[nz++]=pos-ldu;
	if( r<M-1 )	neibors[nz++]=pos+ldu;

	return nz;
}

/*
	ISSUE-NEEDED:
		1 以后可用模板代替

	//{-1,1,-ldu,ldu}
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/5/2006	
*/
double GIP_LAPLAS_( int pos,int M,int N,int ldu,GIP_FLOATING *u,int flag )	{
	int r=GIP_POS2R(pos,ldu),c=GIP_POS2C(pos,ldu);
	double u1,u2,u3,u4;
	u1=c==0 ? u[pos+1] : u[pos-1];
	u2=c==N-1 ? u[pos-1] : u[pos+1];
	u3=r==0 ? u[pos+ldu] : u[pos-ldu];
	u4=r==M-1 ? u[pos-ldu] : u[pos+ldu];

	return u1+u2+u3+u4-u[pos]*4;
}

/*
	ISSUE-NEEDED:
		1 以后可用模板代替

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/12/2006	
*/
static double _MASK_3_w_sum[6]={
	6,6,8,8,8,8
};
static double GIP_MASK_3_3[6][9]={
	{-1,-1,-1,0,0,0,1,1,1},		//PREWITT_X
	{-1,0,1,-1,0,1,-1,0,1},		//PREWITT_Y
	{-1,-2,-1,0,0,0,1,2,1},		//SOBEL_R			
	{-1,0,1,-2,0,2,-1,0,1},		//SOBEL_C			
	{-2,-1,0,-1,0,1,0,1,2},		//SOBEL_45			
	{0,-1,-2,1,0,-1,2,1,0},		//SOBEL_135			
};
double GIP_Mask_3_( int pos,int M,int N,int ldu,GIP_FLOATING *u,GIP_MASK_OPERATOR type,int flag )	{
	int i,r,c,supp=9,pos_1;
	int GIP_SUPP_r[]={-1,-1,-1,0,0,0,1,1,1},GIP_SUPP_c[]={-1,0,1,-1,0,1,-1,0,1};
	double a,*w;								//{-1,1,-ldu,ldu}

	a = 0.0;
	w=GIP_MASK_3_3[type-G_MASK_3_3];
	for( i = 0; i < supp; i++ )	{
		r=GIP_POS2R(pos,ldu),		c=GIP_POS2C(pos,ldu);
		r+=GIP_SUPP_r[i];			c+=GIP_SUPP_c[i];			
		if( r==-1 )	r=1;
		if( r==M )	r=M-2;
		if( c==-1 )	c=1;
		if( c==N )	c=N-2;
		ASSERT( r>=0 && r<M && c>=0 && c<N );
		pos_1=GIP_RC2POS( r,c,ldu );
		a += w[i]*u[pos_1];
	}
	a /= _MASK_3_w_sum[type-G_MASK_3_3]; 

	return a;
}

/*
	[REF]	http://www.cim.mcgill.ca/~image529/TA529/Image529_99/assignments/edge_detection/references/sobel.htm

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/12/2009	
*/
static double _MASK_5_w_sum[4]={
	18,18,18,18
};
static double _MASK_5_w[4][25]={
	{	//5_5_R		
		0,-1,-1,-1, 0,
		0,-2,-2,-2, 0,
		0, 0, 0, 0, 0,
		0, 2, 2, 2, 0,
		0, 1, 1, 1, 0,
	},			
	{	//5_5_c		
		 0, 0,0,0,0,
		-1,-2,0,2,1,
		-1,-2,0,2,1,
		-1,-2,0,2,1,
		 0, 0,0,0,0,
	},			
	{	//SPOT_R		
		0,-2,-2,-2, 0,
		0,-1,-1,-1, 0,
		0, 0, 0, 0, 0,
		0, 1, 1, 1, 0,
		0, 2, 2, 2, 0,
	},			
	{	//SPOT_c		
		 0, 0,0,0,0,
		-2,-1,0,1,2,
		-2,-1,0,1,2,
		-2,-1,0,1,2,
		 0, 0,0,0,0,
	}
};
double GIP_Mask_5_( int pos,int M,int N,int ldu,GIP_FLOATING *u,GIP_MASK_OPERATOR type,int flag )	{
	int i,r,c,pos_1,r_0,c_0,dr,dc;
	double a,*w;								//{-1,1,-ldu,ldu}

	r_0=GIP_POS2R(pos,ldu),		c_0=GIP_POS2C(pos,ldu);
	a = 0.0;
	w=_MASK_5_w[type-G_MASK_5_5];
	i=0;
	for( dr = -2; dr <= 2; dr++ )	{
	for( dc = -2; dc <= 2; dc++ )	{
		r=r_0+dr;			c=c_0+dc;			
		if( r<0 )	r=-r;
		if( r>=M )	r=2*M-1-r;
		if( c<0 )	c=-c;
		if( c>=N )	c=2*N-1-c;
		ASSERT( r>=0 && r<M && c>=0 && c<N );
		pos_1=GIP_RC2POS( r,c,ldu );
		a += w[i++]*u[pos_1];
	}
	}
	a /= _MASK_5_w_sum[type-G_MASK_5_5];

	return a;
}

/*
	for binary image[0,1]
	{-1,1,-ldu,ldu}	{-ldu-1,-ldu+1,ldu-1,ldu+1}
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/13/2006	
*/
int GIP_adjacent_2(int M,int N,int ldu,int *u,int r,int c,int *adj,GIP_ADJACENT_TYPE type  )	{
	int i,nAdj=0,r_1,c_1,pos_1,pos_2,pos_3,N4=4,N8=8,target;	
	int GIP_SUPP_r[]={0,0,-1,1,-1,-1,1,1},GIP_SUPP_c[]={-1,1,0,0,-1,1,-1,1};

	pos_1=GIP_RC2POS( r,c,ldu );
	target=u[pos_1];
	for( i = 0; i < N8; i++ )	{
		r_1=r+GIP_SUPP_r[i];			c_1=c+GIP_SUPP_c[i];
		if( !GIP_RC_VALID(r_1,0,M-1) || !GIP_RC_VALID(c_1,0,N-1) )
			continue;
		if( i>=N4 && type==GIP_ADJACENT_4 )		//GIP_ADJACENT_4
			break;
		pos_1=GIP_RC2POS( r_1,c_1,ldu );
		if( i>N4 && type==GIP_ADJACENT_M )	{	//GIP_ADJACENT_M
			pos_2 = GIP_RC2POS( r,c_1,ldu );
			pos_3 = GIP_RC2POS( r_1,c,ldu );
			if( u[pos_1]==target && u[pos_2]!=target && u[pos_3]!=target)	{
				adj[nAdj++]=pos_1;
			}
		}else	{								//GIP_ADJACENT_8
			if( u[pos_1]==target )
			{	adj[nAdj++]=pos_1;			}
		}
	}	
	
	return nAdj;
}

/*
	to get the property of sub-pixel
	bilinear interpolation
	v0.1	cys
		1/13/2006	
*/		
double GIP_sub_pixel( int M,int N,int ldu,GIP_FLOATING *u,GIP_FLOATING r_x,GIP_FLOATING c_x,int flag )	{
	int r_o,c_o,isR,isC,pos,r_1,c_1,r_2,c_2;
	double a,a_1,a_2;

	r_o = G_DOUBLE2INT(r_x);			c_o = G_DOUBLE2INT(c_x);	
	isR = fabs(r_o-r_x)<FLT_EPSILON;
	isC = fabs(c_o-c_x)<FLT_EPSILON;

	if( isR && isC )	{
		pos=GIP_RC2POS( r_o,c_o,ldu );
		a = u[pos];
	}else if( isR )	{
		c_1=G_DOUBLE2INT(c_x-0.5);		c_2=G_DOUBLE2INT(c_x+0.5);
		ASSERT( c_2==c_1+1 );
		pos=GIP_RC2POS( r_o,c_1,ldu );
		a_1=u[pos];				a_2=u[pos+1];
		a = a_1*(c_2-c_x)+a_2*(c_x-c_1);
	}else if( isC )	{
		r_1=G_DOUBLE2INT(r_x-0.5);		r_2=G_DOUBLE2INT(r_x+0.5);
		ASSERT( r_2==r_1+1 );
		pos=GIP_RC2POS( r_1,c_o,ldu );
		a_1=u[pos];				a_2=u[pos+ldu];
		a = a_1*(r_2-r_x)+a_2*(r_x-r_1);
	}else	{	//bilinear interpolation
		r_1=G_DOUBLE2INT(r_x-0.5);		r_2=G_DOUBLE2INT(r_x+0.5);
		c_1=G_DOUBLE2INT(c_x-0.5);		c_2=G_DOUBLE2INT(c_x+0.5);
		ASSERT( r_2==r_1+1 && c_2==c_1+1 );
		pos=GIP_RC2POS( r_1,c_1,ldu );
		a_1=(u[pos]*(c_2-c_x)+u[pos+1]*(c_x-c_1) );
		pos=GIP_RC2POS( r_2,c_1,ldu );
		a_2=(u[pos]*(c_2-c_x)+u[pos+1]*(c_x-c_1) );

		a = a_1*(r_2-r_x)+a_2*(r_x-r_1);
	}		

	return a;
}

/*

	沿(dr,dc)方向偏移的第len个像素

	注意
	1 按像素偏移，与一般的偏移概念有区别

	v0.1	cys
		2/21/2009
*/
void GIP_pixel_shift_( int r,int c,double dr,double dc,int len,int *r1,int *c1 )	{
	double a;

	ASSERT( !(dr==0.0 && dc==0.0 ) );
	if( fabs(dr)>fabs(dc) )	{
		len = dr>0.0 ? len : -len;
		*r1 = r+len;
		a = len*dc/dr;
		*c1 = c+G_DOUBLE2INT( a );
	}else	{
		len = dc>0.0 ? len : -len;
		*c1 = c+len;
		a = len*dr/dc;
		*r1 = r+G_DOUBLE2INT( a );
	}
}

#define GIP_MAX_TITLE       260
static TCHAR g_sPath[GIP_MAX_TITLE];
/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/22/2009	
*/
TCHAR *GIP_path_( TCHAR *sTitle,int no )	{
	if( no==-1 )	{
		_stprintf( g_sPath,_T(".\\trace\\%s.jpg"),sTitle );
	}else	{
		_stprintf( g_sPath,_T(".\\trace\\%s_%d.jpg"),sTitle,no );
	}

	return g_sPath;
}

/*
	暂采用冒泡排序,arr被改动

	ISSUE-NEEDED:
		需采用类似快速排序的算法

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/29/2009	
*/
void GIP_arr_top( GIP_FLOATING *arr,int len,int nTop,int *arrPos,int flag )	{
	int i,j,*order=GIP_alloc( sizeof(int)*len ),t;
	double a;

	ASSERT( nTop<len );
	for( i = 0; i < len; i++ )	
		order[i]=i;

	for( i = 0; i < len; i++ )	{
	for( j = i+1; j < len; j++ )	{
		if( arr[i]<arr[j] ){
			a=arr[i];		arr[i]=arr[j];		arr[j]=a;
			t=order[i];		order[i]=order[j];	order[j]=t;
		}
	}
	}
	for( i=0; i<nTop; i++ )	{
		arrPos[i] = order[i];
	}
	GIP_free( order );
}

/*
	将arr快速分成两份[low,SEP] <= [SEP+1,high]
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/6/2009	
*/
void GIP_quick_partion( GIP_FLOATING *arr,int low,int high,int SEP,int flag )	{
	int scanUp,scanDown;
	GIP_FLOATING pivot,a;

	if( high-low <= 0 )
		return;
	else if ( high-low==1 )	{
		if( arr[high]<arr[low] )
			GIP_SWAP( arr[high],arr[low],a );
		return;
	}

	pivot = arr[SEP];
	GIP_SWAP( arr[SEP],arr[low],a );
	scanUp = low+1;		scanDown=high;
	do{
		while( scanUp<=scanDown && arr[scanUp]<=pivot )
			scanUp++;
		while(pivot < arr[scanDown] )
			scanDown--;
		if( scanUp<scanDown )
			GIP_SWAP( arr[scanUp],arr[scanDown],a );

	}while( scanUp<scanDown );
	arr[low]=arr[scanDown];
	arr[scanDown]=pivot;
	
	if( scanDown>SEP )
		GIP_quick_partion( arr,low,scanDown-1,SEP,flag );
	else if( scanDown<SEP )	
		GIP_quick_partion( arr,scanDown+1,high,SEP,flag );
	else
		return;

}
#include<windows.h>
#include <memory.h>
#include <FLOAT.h>
#include <stdio.h>
#include <time.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/gip_util.h"
#include "../util/gip_ss.h"
#include "../util/gip_fv.h"
#include "../util/gip_thread.h"
#include "GIP_library.h"

/*
	仅采用MAJOR控制版本，以后可扩充
	3/12/2009
*/
#define _LIB_VERSION_MAJOR	1
#define _LIB_VERSION_MINOR	0

static char BLAS_N[1]={'N'},BLAS_T[1]={'T'},BLAS_C[1]={'C'},BLAS_L[1]={'L'},BLAS_R[1]={'R'},BLAS_U[1]={'U'};
static int inc_1=1; 
static GIP_FLOATING one_=1.0,fuyi_=-1.0,zero_=0.0;		 
static TCHAR drive[_MAX_DRIVE],dir[_MAX_DIR],fname[_MAX_FNAME],ext[_MAX_EXT];

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/19/2009	
*/
G_FILE_TYPE GIP_File_Type( TCHAR *ext )		{
	G_FILE_TYPE type;

	if( _tcscmp (ext,_T(".jpg") )==0x0 )	{
		type = G_IMAGE_JPG;
	}else if( _tcscmp (ext,_T(".bmp") )==0x0 )	{
		type = G_IMAGE_BMP;
	}else if( _tcscmp (ext,_T(".pgm") )==0x0 )	{
		type = G_IMAGE_PGM;
	}else if( _tcscmp (ext,_T(".") )==0x0 )	{
		type = G_FILE_DIRECTORY;
	}else
		type = G_FILE_UNKNOWN;

	return type;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/19/2009	
*/
void GIP_LIB_Append( GIP_LIB *hLib,TCHAR *sPath,int no,int cls,int flag )	{
	GIP_LIB_entry *hEntry,*arrEntry=hLib->arrEntry;

	if( hLib->nEntry==hLib->nMaxEntry )	{
		hLib->nMaxEntry *= 2;
		arrEntry = GIP_alloc( sizeof(GIP_LIB_entry)*hLib->nMaxEntry );
		GIP_MEMCOPY( arrEntry,hLib->arrEntry,sizeof(GIP_LIB_entry)*hLib->nEntry );
		GIP_free( hLib->arrEntry );
		hLib->arrEntry = arrEntry;
	}
/*
		if( (src_image=GIP_LoadImage(sTestPath,iscolor))==0 )    {
			G_PRINTF("Image was not loaded.\n");
			return -1;
		}
*/	hEntry=hLib->arrEntry+hLib->nEntry;
	GIP_MEMCLEAR( hEntry,sizeof(GIP_LIB_entry) );
	ASSERT( _tcslen(sPath)<MAX_PATH_LEN );
	_tcscpy( hEntry->sFilePath,sPath );
//	_tsplitpath( hEntry->sFilePath, drive, dir, fname, ext );
//	hEntry->sTitle=GIP_alloc( sizeof(char)*(_tcslen(fname)+1) );
//	_tcscpy( hEntry->sTitle,fname );
	hEntry->no = no;
	hEntry->cls = cls;

	hLib->nEntry++;
	ASSERT( hLib->nEntry<=hLib->nMaxEntry );	
}

/*
	http://www.diybl.com/course/3_program/c++/cppjs/2008331/107762.html

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/19/2009	
*/
void _List_File( GIP_LIB *hLib,int flag )		{
    HANDLE hSearch;
    WIN32_FIND_DATA data;
	G_FILE_TYPE type;
	int no=0;

    hSearch=FindFirstFile(_T("*"),&data);
    do{
        if(data.dwFileAttributes==FILE_ATTRIBUTE_DIRECTORY
           &&_tcscmp (data.cFileName,_T("."))  && _tcscmp (data.cFileName,_T("..")) )	{
            SetCurrentDirectory(data.cFileName);
            _List_File( hLib,flag );
            SetCurrentDirectory(_T("..") );
        }
        else if(_tcscmp (data.cFileName,_T(".")) && _tcscmp (data.cFileName,_T("..")))	{
 			_tsplitpath( data.cFileName, drive, dir, fname, ext );
            type = GIP_File_Type( ext );
			if( type>G_IMAGE_BEGIN && type<G_IMAGE_END ) {
				TCHAR sLine[MAX_LINE];
				GetCurrentDirectory( MAX_LINE, sLine);
				_stprintf( sLine,_T("%s\\%s"),sLine,data.cFileName );
				GIP_LIB_Append( hLib,sLine,no,hLib->nClass,0x0 );
				no++;
			}
        }
	}while(FindNextFile(hSearch,&data));
	if( no>0 )
 		hLib->nClass++;

    FindClose(hSearch);
}

/*
	注意:
		hfv直接复制给hLib->fv,因此hfv->V指向的内存由GIP_LIB负责释放 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/25/2008	
*/
int GIP_LIB_init( GIP_LIB *hLib,TCHAR *sPath,G_FVector_ *hfv_0,double *control,int flag )	{
	TCHAR sLine[MAX_LINE],sCurPath[MAX_LINE];
	FILE*  fp;
	int i,nMaxEntry=1000,len,model,method,nLoop,isGetParameter=0,M_normal,N_normal;
	G_FILE_TYPE type;
//	GIP_LIB_entry *hEntry;

	memset( hLib,0x0,sizeof(GIP_LIB) );

	ASSERT( control!=GIP_NULL );
	for( i = 0; i < GIP_USER_CONTROL; i++ )
		hLib->user_control[i]=control[i];

	GetCurrentDirectory( MAX_LINE,sCurPath );
	hLib->version = _LIB_VERSION_MAJOR;		
	hLib->D_span = 1000;
	hLib->D_inter=GIP_calloc( sizeof(float),(hLib->D_span+1) );
	hLib->D_intra=GIP_calloc( sizeof(float),(hLib->D_span+1) );

	hLib->arrEntry = GIP_alloc( sizeof(GIP_LIB_entry)*nMaxEntry );
	hLib->nMaxEntry = nMaxEntry;

	_tcscpy( hLib->sPath,sPath );
	if( BIT_TEST( flag,GIP_PATH_DIRECTORY ) )	{
//		G_PRINTF( _T("directory is %s"), sPath );
		SetCurrentDirectory(sPath);
		_List_File( hLib,0x0 );
	}else if( BIT_TEST( flag,GIP_PATH_INFILE ) ){
		fp=_tfopen( sPath,"r" );
		if( fp==NULL )	{
			G_PRINTF( "error at open file: %s", sPath );
			return -1;
		} 

		while( fReadLine( fp,sLine ) != 0 )	{
			len=(int)(_tcslen( sLine ));
			if( len==0 )			continue;
			if( sLine[0]=='c' )			
			{	G_PRINTF( "%s",sLine );			continue;			}
			if( sLine[0]=='!' )		continue;			//注释行		
			if( sLine[0]=='q' && sLine[1]=='u' && sLine[2]=='i' && sLine[3]=='t' )		{
				break;
			}
			if( sLine[0]=='p' )		{	
				if( _stscanf( sLine+1,"%d %d %d",&model,&method,&nLoop )!=3 )
				{	break;		}
				ASSERT( nLoop > 0 );
				isGetParameter=1;
				continue;	
			}
			if( isGetParameter==0 )		
				continue;
			_tsplitpath( sLine, drive, dir, fname, ext );
			type = GIP_File_Type( ext );
			if( type==G_FILE_DIRECTORY )	{
				SetCurrentDirectory(sLine);
				_List_File( hLib,0x0 );
			}else if( type>G_IMAGE_BEGIN && type<G_IMAGE_END ) {
				GIP_LIB_Append( hLib,sLine,0,hLib->nClass,0x0 );
				hLib->nClass++;
			}
		}
		fclose( fp );	
	}else	{
		_tsplitpath( sPath, drive, dir, fname, ext );
		type = GIP_File_Type( ext );
		if( type>G_IMAGE_BEGIN && type<G_IMAGE_END ) {
			GIP_LIB_Append( hLib,sPath,0,hLib->nClass,0x0 );
			hLib->nClass++;
		}
	}

    SetCurrentDirectory(sCurPath);

	hLib->hfv = hfv_0;
	ASSERT( hfv_0!=GIP_NULL );
	hLib->fv_dat = (G_U8*)GIP_alloc( hfv_0->unit_size*hfv_0->V_len*hLib->nEntry );		
//	G_PRINTF( "%s load in.", sLibFile );

	return 0x0;
}

/*
	注意:
		hfv直接复制给hLib->fv,因此hfv->V指向的内存由GIP_LIB负责释放 

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/24/2009	
*/
int GIP_LIB_extract( GIP_LIB *hLib,GIP_LIB *hLib_0,TCHAR *sPath,G_FVector_ *hfv_0,int *mark,int flag )	{
	TCHAR sLine[MAX_LINE],sCurPath[MAX_LINE];
	FILE*  fp;
	double *control=hLib_0->user_control;
	int i,nMaxEntry=hLib_0->nEntry,nEntry,pos,cls;
	G_FILE_TYPE type;
	GIP_LIB_entry *hEntry,*hEntryOld;
	TCHAR *drive_0,*dir_0,*drive_1,*dir_1,fname[_MAX_FNAME],ext[_MAX_EXT];

	drive_0=GIP_alloc( sizeof(TCHAR)*(_MAX_DRIVE+_MAX_DIR)*2 );			dir_0=drive_0+_MAX_DRIVE;
	drive_1=dir_0+_MAX_DIR;			dir_1=drive_1+_MAX_DRIVE;

	memset( hLib,0x0,sizeof(GIP_LIB) );
	GIP_MEMCOPY( hLib,hLib_0,sizeof(GIP_LIB) );

	_tcscpy( hLib->sPath,sPath );
	hLib->D_inter=GIP_alloc( sizeof(float)*(hLib->D_span+1) );
	hLib->D_intra=GIP_alloc( sizeof(float)*(hLib->D_span+1) );
	nEntry=0;
	for( i = 0; i < nMaxEntry; i++ )	{
		if( mark[i]==0 )	continue;
		nEntry++;
	}
	hLib->arrEntry = GIP_alloc( sizeof(GIP_LIB_entry)*nEntry );
	hLib->nMaxEntry = nEntry;

	pos=0;
	for( i = 0; i < nMaxEntry; i++ )	{
		if( mark[i]==0 )	continue;
		hEntry = hLib->arrEntry+pos;
		GIP_MEMCOPY( hEntry,hLib_0->arrEntry+i,sizeof(GIP_LIB_entry) );
		hEntry->no=pos;			hEntry->cls=-1;
		pos++;
	}
	ASSERT( pos==nEntry );
//重新对cls排序
	cls = 0;
	hEntryOld = hLib->arrEntry;
	hEntryOld->cls = 0;
	_wsplitpath_s( hEntryOld->sFilePath, drive_0,_MAX_DRIVE,dir_0,_MAX_DIR,fname,_MAX_FNAME,ext,_MAX_EXT );
	for( i = 1; i < nMaxEntry; i++ )	{
		hEntry = hLib->arrEntry+i;
		_wsplitpath_s( hEntry->sFilePath, drive_1,_MAX_DRIVE,dir_1,_MAX_DIR,fname,_MAX_FNAME,ext,_MAX_EXT );
		if( _tcscmp(drive_0,drive_1)!=0 || _tcscmp(dir_0,dir_1)!=0 )	{
			hEntry->cls = cls++;
		}
		hEntryOld = hEntry;			
		drive_1=drive_0;				dir_1=dir_0;	
	}
	hLib->hfv = hfv_0;
	ASSERT( hfv_0!=GIP_NULL );
	hLib->fv_dat = (G_U8*)GIP_alloc( hfv_0->unit_size*hfv_0->V_len*hLib->nEntry );		

	GIP_free( drive_0 );
//	G_PRINTF( "%s load in.", sLibFile );

	return 0x0;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/25/2008			
		3/9/2009	
			G_FVector_独立于GIP_LIB
*/
void GIP_LIB_clear( GIP_LIB *hLib )	{
	int i;
	GIP_LIB_entry *hEntry;
	if( hLib == GIP_NULL )	return;

//	G_FVector_clear( &(hLib->fv) );			//G_FVector_独立于GIP_LIB

	if( hLib->arrEntry != GIP_NULL )	{
		for( i = 0; i < hLib->nEntry; i++ )	{
			hEntry=hLib->arrEntry+i;
//			GIP_free( hEntry->sTitle );			hEntry->sTitle=GIP_NULL;
		}
		GIP_free( hLib->arrEntry );				hLib->arrEntry=GIP_NULL; 
	}
	hLib->nEntry = -1;
/*
	if( hLib->W_ptr != GIP_NULL )	{
		GIP_free( hLib->W_ptr );			hLib->W_ptr=GIP_NULL;
	}
	if( hLib->W_no != GIP_NULL )	{
		GIP_free( hLib->W_no );			hLib->W_no=GIP_NULL;
	}*/
	if( hLib->fv_dat != GIP_NULL )	{
		GIP_free( hLib->fv_dat );			hLib->fv_dat=GIP_NULL;
	}
	if( hLib->D_inter != GIP_NULL )	{
		GIP_free( hLib->D_inter );			hLib->D_inter=GIP_NULL;
	}
	if( hLib->D_intra != GIP_NULL )	{
		GIP_free( hLib->D_intra );			hLib->D_intra=GIP_NULL;
	}

	GIP_MEMCLEAR( hLib,sizeof(GIP_LIB) );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/20/2009	
*/
int GIP_LIB_scan( GIP_LIB *hLib,int flag )	{
	int i,nEntry=hLib->nEntry,cls=-1,nz=0,W_dim;
//	int *W_ptr=hLib->W_ptr,*W_no=hLib->W_no;
	int *W_ptr,*W_no;
	GIP_LIB_entry *hEntry;

	if( W_ptr!=GIP_NULL )
	{	GIP_free(W_ptr);		W_ptr=GIP_NULL;		}
	if( W_no!=GIP_NULL )
	{	GIP_free(W_no);			W_no=GIP_NULL;		}
	W_ptr=GIP_alloc( sizeof(int)*(hLib->nClass+1) );		//hLib->W_ptr=W_ptr;
	W_no=GIP_alloc( sizeof(int)*nEntry );					//hLib->W_no=W_no;

	W_ptr[0]=0;
	W_dim=0;
	for( i = 0; i < nEntry; i++ )	{
		hEntry = hLib->arrEntry+i;
		if( hEntry->status != 1 )	
			continue;		//处理失败
		if( cls==-1 )		cls=hEntry->cls;
		if( hEntry->cls != cls )	{
			W_dim++;		
			W_ptr[W_dim]=nz;
			cls = hEntry->cls;
		}
		W_no[nz] = i;
		nz++;
	}
	if( nz>W_ptr[W_dim] )	{
		W_dim++;		
		W_ptr[W_dim]=nz;
	}
	ASSERT( nz == hLib->nSample && W_dim<=hLib->nClass );

//	hLib->W_dim = W_dim;

	return 0x0;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/26/2008	
*/
void GIP_LIB_reducing( GIP_LIB *hLib,GIP_FLOATING *sample_0, int flag )	{
	int m,nSample = hLib->nSample,N,i,j,u_block;
//	int *W_ptr=hLib->W_ptr,*W_no=hLib->W_no;
	int *W_ptr,*W_no;
	double *W_pca,*W_lda,*mean,*u_mean,*u_s,*u;
	GIP_FLOATING *sample,*S_b,*S_w;

//	GIP_LIB_scan( hLib,0x0 );

	N = hLib->fv_n_0;
	m = hLib->nClass;
	mean=GIP_alloc( sizeof(GIP_FLOATING)*N );
	W_pca=GIP_alloc( sizeof(GIP_FLOATING)*N*m );
//	hLib->W=GIP_alloc( sizeof(GIP_FLOATING)*N*m );
//	GIP_CalcPCA( nSample,N,m,sample,mean,W_pca );		//get PCA matrix
//reduce dimension
	BLAS_GEMM_0( BLAS_T,BLAS_N,m,nSample,N,one_,W_pca,N,sample_0,N,zero_,sample,m );		//dgemm( sample_0,W_pca,sample )
//	get S_b,S_w
	u=GIP_alloc( sizeof(GIP_FLOATING)*m*(m+2) );		
	u_mean=u+m*m;			u_s=u_mean+m;
	for( i = 0; i < m; i++ )	{	//get mean for each class
		u_block = W_ptr[i+1]-W_ptr[i];
		for( j = 0; j < u_block; j++ )	u_s[j]=1.0/u_block;
		BLAS_GEMV_0( BLAS_N,N,u_block,one_,sample+W_ptr[i]*m,m,u_block,inc_1,one_,u+i*m,inc_1 );				
	}
	for( j = 0; j < m; j++ )	u_s[j]=1.0/m;
	BLAS_GEMV_0( BLAS_N,m,m,one_,u,m,u_s,inc_1,one_,u_mean,inc_1 );				
	for( i = 0; i < m; i++ )	{	
		//SB += (W_ptr[i+1]-W_ptr[i])*(ui-u_mean)*(ui-u_mean)'
		for( j = W_ptr[i]; j < W_ptr[i+1]; j++ )	{
		//SW += (u-ui)*(u-ui)'
		}
	}
//	get W_lda(general eigenvalue)
//	dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info);
//	dgemm( W_pca,W_lda,hLib->W )
//	dgemm( sample,hLib->W,hLib->fv );
	GIP_free( W_pca );				GIP_free( u );

}

/*
	基于用户需要，允许空LIB(nEntry=0)保存

	教训：
		1 采用"w"选项，写入无格式数据，导致stream中包含多个FEOF，读取失败。改为"wb"选项。	2/13/2009

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/20/2009	
	v0.2	cys	arrEntry全部存储，不再记录W_*
		3/6/2009	
*/
int GIP_LIB_save( GIP_LIB *hLib,TCHAR *sLibPath, int flag )	{
	int nEntry,N,ret=0,D_span=hLib->D_span;
//	int W_dim=hLib->W_dim,*W_ptr=hLib->W_ptr,*W_no=hLib->W_no;
	GIP_LIB_entry *arrEntry=hLib->arrEntry;
	FILE *fp;
	char *fv_dat=hLib->fv_dat;
	long pos;
	G_FVector_ *hfv=hLib->hfv;

	nEntry=hLib->nEntry;
//	W_nz=W_ptr[W_dim];
	fp=_tfopen( sLibPath,_T("wb") );
	if( fp==0x0 )
	{	ret=-1;		goto _FIO_EXIT_;		}
	if( fwrite( hLib->user_control,sizeof(double),GIP_USER_CONTROL,fp )!=GIP_USER_CONTROL )
	{	ret=-10;		goto _FIO_EXIT_;	}
	if( fwrite( hLib->sRoot,sizeof(TCHAR),MAX_PATH_LEN,fp )!=MAX_PATH_LEN )
	{	ret=-11;		goto _FIO_EXIT_;	}
	GIP_FWRITE_1( hLib->version,fp );

	GIP_FWRITE_1( hLib->CRR,fp );
	GIP_FWRITE_1( hLib->L_spd,fp );
	GIP_FWRITE_1( hLib->M_spd,fp );

	GIP_FWRITE_1( hLib->nEntry,fp );
	GIP_FWRITE_1( hLib->nSample,fp );
	GIP_FWRITE_1( hLib->nMaxEntry,fp );
	GIP_FWRITE_1( hLib->nClass,fp );
/*	GIP_FWRITE_1( hLib->W_dim,fp );
	if( fwrite( W_ptr,sizeof(int),W_dim+1,fp )!=W_dim+1 ) 
	{	ret=-3;		goto _FIO_EXIT_;	}
	if( fwrite( W_no,sizeof(int),W_nz,fp )!=W_nz )
	{	ret=-4;		goto _FIO_EXIT_;	}*/
	if( fwrite( arrEntry,sizeof(GIP_LIB_entry),nEntry,fp )!=nEntry )
	{	ret=-5;		goto _FIO_EXIT_;	}

	N=hfv->V_len;
	GIP_FWRITE_1( hfv->type,fp );
	GIP_FWRITE_1( hfv->flag,fp );
	GIP_FWRITE_1( hfv->V_len,fp );
	GIP_FWRITE_1( hfv->sf_M,fp );
	GIP_FWRITE_1( hfv->sf_N,fp );
	GIP_FWRITE_1( hfv->M_step,fp );
	GIP_FWRITE_1( hfv->N_step,fp );
	GIP_FWRITE_1( hfv->bpSF,fp );
	GIP_FWRITE_1( hfv->nzSF,fp );
	GIP_FWRITE_1( hfv->unit_size,fp );
	GIP_FWRITE_1( hfv->nShift,fp );

	pos = ftell( fp );
	if( fwrite( fv_dat,hfv->unit_size,N*nEntry,fp ) != N*nEntry )
	{	ret=-6;		goto _FIO_EXIT_;	}
	pos = ftell( fp );
	GIP_FWRITE_1( D_span,fp );
	if( fwrite( hLib->D_inter,sizeof(float),(D_span+1),fp )!=D_span+1 )
	{	ret=-8;		goto _FIO_EXIT_;	}
	if( fwrite( hLib->D_intra,sizeof(float),(D_span+1),fp )!=D_span+1 )
	{	ret=-8;		goto _FIO_EXIT_;	}

_FIO_EXIT_:
	if( fp != 0x0 && fclose( fp )!= 0 )
	{	ret=-7;			}
	if( ret==0 )	{
		G_PRINTF( _T("SAVE %s:nEntry=%d"),sLibPath,nEntry );
	}else	{
		ret = fp!=0x0 ? ferror( fp ) : ret;
		G_PRINTF( _T("\t!!!Failed to save %s. err=%d"),sLibPath,ret );
	}

	if( ret==0 && 0 )	{		//!!testing!!
		GIP_LIB *hLib;
		hLib=GIP_alloc( sizeof(GIP_LIB) );
		GIP_LIB_load( hLib,hfv,sLibPath,0x0 );
		GIP_LIB_clear( hLib );			GIP_free( hLib );
	}

	return ret;
}

/*
	注意：
		1 必须先释放hLib与hfv
		2 load后相当于GIP_LIB_init与G_FVector_init
		3 存在失败的可能，必须返回错误代码

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/20/2009	
	v0.2	cys	arrEntry全部存储，不再记录W_*
		3/6/2009	
		3/9/2009		
			G_FVector_独立于GIP_LIB	
*/
int GIP_LIB_load( GIP_LIB *hLib,G_FVector_ *hfv,TCHAR *sLibPath, int flag )	{
	int N,ret=0,nz,D_len;
	long pos;
	FILE *fp;

	fp=_tfopen( sLibPath,_T("rb") );
	if( fp==0x0 )
	{	ret=-1;		goto _FIO_EXIT_;	}
	
	GIP_LIB_clear( hLib );	
	G_FVector_clear( hfv );
	if( fread( hLib->user_control,sizeof(double),GIP_USER_CONTROL,fp )!=GIP_USER_CONTROL )
	{	ret=-10;		goto _FIO_EXIT_;	}
	if( fread( hLib->sRoot,sizeof(TCHAR),MAX_PATH_LEN,fp )!=MAX_PATH_LEN )
	{	ret=-11;		goto _FIO_EXIT_;	}
	GIP_FREAD_1( hLib->version,fp );
	if( hLib->version<=0 || hLib->version>_LIB_VERSION_MAJOR )	{	//版本错误
		ret=GIP_LIB_VERSION_ERR;		goto _FIO_EXIT_;	
	}else if( hLib->version!=_LIB_VERSION_MAJOR  )	{				//旧版本
		ASSERT( FALSE );		
	}
	GIP_FREAD_1( hLib->CRR,fp );
	GIP_FREAD_1( hLib->L_spd,fp );
	GIP_FREAD_1( hLib->M_spd,fp );
	GIP_FREAD_1( hLib->nEntry,fp );
	GIP_FREAD_1( hLib->nSample,fp );
	GIP_FREAD_1( hLib->nMaxEntry,fp );
	GIP_FREAD_1( hLib->nClass,fp );
	if( hLib->CRR<0.0 || hLib->CRR>1.0 || hLib->nEntry<0 || hLib->nClass>hLib->nEntry )	
	{	ret=-11;		goto _FIO_EXIT_;	}
/*	GIP_FREAD_1( hLib->W_dim,fp );
	W_dim=hLib->W_dim;				ASSERT( W_dim>0 );
	hLib->W_ptr=GIP_alloc( sizeof(int)*(W_dim+1) );
	if( fread( hLib->W_ptr,sizeof(int),W_dim+1,fp )!=W_dim+1 ) 
	{	ret=-3;		goto _FIO_EXIT_;	}
	W_nz=hLib->W_ptr[W_dim];
	ASSERT( W_nz>=0 && W_nz<=hLib->nEntry );
	hLib->W_no=GIP_alloc( sizeof(int)*W_nz );
	if( fread( hLib->W_no,sizeof(int),W_nz,fp )!=W_nz )  
	{	ret=-4;		goto _FIO_EXIT_;	}
*/
	hLib->arrEntry=GIP_alloc( sizeof(GIP_LIB_entry)*hLib->nEntry );
	if( fread( hLib->arrEntry,sizeof(GIP_LIB_entry),hLib->nEntry,fp )!=hLib->nEntry )  
	{	ret=-5;		goto _FIO_EXIT_;	}

//	hfv=GIP_alloc( sizeof(G_FVector_) );			
	hLib->hfv = hfv;

	GIP_FREAD_1( hfv->type,fp );
	if( hfv->type<0 || hfv->type>=G_FV_TYPE_END )	
	{	ret=-7;		goto _FIO_EXIT_;	}
	GIP_FREAD_1( hfv->flag,fp );
	GIP_FREAD_1( hfv->V_len,fp );
	GIP_FREAD_1( hfv->sf_M,fp );
	GIP_FREAD_1( hfv->sf_N,fp );
	GIP_FREAD_1( hfv->M_step,fp );
	GIP_FREAD_1( hfv->N_step,fp );
	GIP_FREAD_1( hfv->bpSF,fp );
	GIP_FREAD_1( hfv->nzSF,fp );
	GIP_FREAD_1( hfv->unit_size,fp );
	GIP_FREAD_1( hfv->nShift,fp );

	hfv->V = GIP_alloc( hfv->unit_size*hfv->V_len );
	N=hfv->V_len;
	hLib->fv_dat=(char*)GIP_alloc( hfv->unit_size*N*hLib->nEntry );
	pos = ftell( fp );
	nz = fread( hLib->fv_dat,hfv->unit_size,N*hLib->nEntry,fp );
	if( nz !=N*hLib->nEntry )  	{
		ret = ferror( fp );
		ret = feof( fp );
		ret=-6;		goto _FIO_EXIT_;	
	}
	if( 0 )	{
		hLib->D_span = 1000;
		hLib->D_inter=GIP_alloc( sizeof(float)*(hLib->D_span+1) );
		hLib->D_intra=GIP_alloc( sizeof(float)*(hLib->D_span+1) );
	}else	{
		GIP_FREAD_1( hLib->D_span,fp );
		if( hLib->D_span<0 || hLib->D_span>=1000000 )	
		{	ret=-9;		goto _FIO_EXIT_;	}
		D_len = hLib->D_span+1;
		hLib->D_inter=GIP_alloc( sizeof(float)*D_len );
		hLib->D_intra=GIP_alloc( sizeof(float)*D_len );
		if( fread( hLib->D_inter,sizeof(float),D_len,fp )!=D_len )
		{	ret=-8;		goto _FIO_EXIT_;	}
		if( fread( hLib->D_intra,sizeof(float),D_len,fp )!=D_len )
		{	ret=-8;		goto _FIO_EXIT_;	}
	}

_FIO_EXIT_:
	if( fp != 0x0 && fclose( fp )!= 0 )
	{	ret=-7;			}
	if( ret==0 )	{
		_tcscpy( hLib->sPath,sLibPath );
		G_PRINTF( _T("LOAD %s"),sLibPath );	
	}else	{
		G_PRINTF( _T("\t!!!Failed to load %s. err=%d"),sLibPath,ret );
	}	

	return ret;
}

/*
[**]
	看狄龙之边缘岁月
	“作兄弟的，有今生，无来世”   嘿嘿。。。。
[**]

	注意：
	1 基于\MMU iris database\44\.的测试，nShift的最佳值为12

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/7/2009	
	v0.2	cys		返回distri(distribution具有重要的意义)
		3/10/2009	
		4/14/2009	增加dist_4，用于返回dist信息
*/
int GIP_LIB_discriminate_0( GIP_LIB *hLib,G_DCRIMI_ *hDcrimi,double *dist_4,int flag )	{
	GU_CTRL_DATA *hUserCtrl=g_hUserCtrl;
	int nEntry=hLib->nEntry,nCls=hLib->nClass,V_len,*C_ptr;
	int i,j,k,no_e0,fail_1,fail_2,D_span=hLib->D_span,D_no,nShift,no;
	int mode=(int)(hLib->user_control[GIRIS_MODE_]);
	G_FVector_ *hfv=hLib->hfv;
	G_FV_TYPE type=hfv->type;
//	int W_dim=hLib->W_dim,*W_ptr=hLib->W_ptr,*W_no=hLib->W_no,W_nz;
	double *inter_0,*inter_1,*intra_0,*intra_1,d,d_min,d_max,s,err_1,err_2,s_D,*w_sort,*dis;
	float *D_inter=hLib->D_inter,*D_intra=hLib->D_intra;
	G_U8 *temp,*mask;
	clock_t start=clock( );
	GIP_LIB_entry *hEntry,*arrEntry=hLib->arrEntry;
	
	nShift = hfv->nShift;
	nShift = G_DOUBLE2INT(nShift*1.0/180*512/2);
	V_len = hfv->V_len;
	if( hfv->V!=GIP_NULL )		GIP_free(hfv->V);
	hfv->V=GIP_alloc( sizeof(G_U8)*V_len*(2*nShift+1) );
	if( hfv->mask!=GIP_NULL )		GIP_free(hfv->mask);
	hfv->mask=GIP_calloc( sizeof(G_U8),V_len*(2*nShift+1) );
//	W_nz = W_ptr[W_dim];
	C_ptr = GIP_alloc( sizeof(int)*(nCls+1) );
	temp = GIP_alloc( MAX(hfv->unit_size*V_len,nEntry*sizeof(int)) );
	mask = GIP_calloc( hfv->unit_size,V_len*nEntry );
	if( dist_4==GIP_NULL )	{
		inter_0 = GIP_alloc( sizeof(double)*nEntry*4 );
	}else	{
		inter_0 = dist_4;
	}
	inter_1=inter_0+nEntry;	intra_0=inter_1+nEntry;	intra_1=intra_0+nEntry;

	fail_1=0,				fail_2=0;
	for( i = 0; i < nEntry; i++ )	{
		intra_0[i]=0.0;		intra_1[i]=0.0;
		inter_0[i]=0.0;		inter_1[i]=0.0;
	}
	for( i = 0; i <= D_span; i++ )	{
		D_inter[i]=0.0;		D_intra[i]=0.0;
	}
	w_sort=GIP_alloc( sizeof(double)*nEntry*2 );		dis=w_sort+nEntry;
	for( i = 0; i < nEntry; i++ )		w_sort[i]=0.0;
	s_D=D_span/1.0;		//适用于HAMMING DISTANCE
	GU_SetUserCtrl( hUserCtrl,GIP_NULL,0,0x0 );
	G_PRINTF( _T(" title\t  cls  no   intra_0 intra_1 inter_0 inter_1\t  off") );
	C_ptr[0]=0;
	for( i = 0; i < nCls; i++ )	{
		j=C_ptr[i]+1;
		while( j<nEntry && arrEntry[j].cls==i )	j++;
		C_ptr[i+1] = j;
	//	if( W_dim>69 && i!=69 )	 continue;
		for( j=C_ptr[i]; j < C_ptr[i+1]; j++ )	{
			if( hUserCtrl!=GIP_NULL && hUserCtrl->bTerminate==TRUE )
				break;
			no = j;
			hEntry = hLib->arrEntry+no;
			if( hEntry->status!=G_ENTRY_OK )	{
				fail_1++;			fail_2++;
				continue;
			}
		//	G_fv_shift_( nShift,512/4,2,V_len,hLib->fv_dat+no*V_len,hfv->V,0x0 );
			G_fv_shift_( hfv->type,nShift,V_len,hLib->fv_dat+no*V_len,hfv->V,0x0 );
			G_hamming_dist_( nEntry,dis,hLib->fv_dat,mask,hfv,GIP_NULL,0x0 );	
			if( C_ptr[i+1]-C_ptr[i]>1 )		{
				d_min=DBL_MAX;			d_max=0.0;
				for( k = C_ptr[i]; k < C_ptr[i+1]; k++ )	{		//intra_class distance
					if( k == j || arrEntry[k].status!=G_ENTRY_OK )		continue;
//					d = G_fv_distance( hLib->fv_dat,mask,hfv,no,k,temp );
					d = dis[k];
					if( d>d_max )		d_max = d;
					if( d<d_min )		d_min = d;
					D_no=G_DOUBLE2INT( d*s_D );
					D_intra[D_no]++;
				}
				//	D_no=G_DOUBLE2INT( d_min*s_D );
				//	D_intra[D_no]++;
				intra_0[j]=d_min;		intra_1[j]=d_max;
			}
//			w_sort[no]=intra_0[j]==0.0 ? 0.0 : intra_1[j]/intra_0[j];
//			G_PRINTF( _T("no=%d,w=%g,path=%s"),no,w_sort[no],hEntry->sFilePath );
			no_e0=-1;
			if( nCls>1 )	{
				d_min=DBL_MAX;			d_max=0.0;
				for( k = 0; k < nEntry; k++ )	{					//inter_class distance
					if( k>=C_ptr[i] && k<C_ptr[i+1] )		continue;
					if( arrEntry[k].status!=G_ENTRY_OK )		continue;
//					d = G_fv_distance( hLib->fv_dat,mask,hfv,no,k,temp );
					d = dis[k];
					D_no=G_DOUBLE2INT( d*s_D );
					D_inter[D_no]++;
					if( d>d_max )		d_max = d;
					if( d<d_min )		
					{	no_e0=k;	d_min = d;	}
				}
				inter_0[j]=d_min;		inter_1[j]=d_max;
			}
			w_sort[no]=inter_0[j]==0.0 ? 0.0 : intra_0[j]/inter_0[j];
			w_sort[no]=inter_0[j];
			hEntry->x = w_sort[no];
			s=intra_1[j]==0.0 ? -1: inter_0[j]/intra_1[j];
			if( intra_0[j]>=inter_0[j] )		fail_1++;
			if( intra_1[j]>=inter_0[j] )		fail_2++;
			G_PRINTF( _T("  I_%d\t%4d %3d %7.3g %7.3g %7.3g(%d) %7.3g\t%5.3g"),j,i,j-C_ptr[i],
				intra_0[j],intra_1[j],inter_0[j],no_e0,inter_1[j],s );
			GU_SetUserCtrl( hUserCtrl,GIP_NULL,j+1,0x0 );
		}
	}
#ifdef _GIRIST_DEVELOP_
	GIP_LIB_sort_( hLib,w_sort,(int*)(temp),0x0 );
#endif
	err_1 = fail_1*1.0/nEntry;				err_2 = fail_2*1.0/nEntry;
	G_DCRIMI_set( (void*)hLib,hDcrimi,D_span,D_inter,D_intra,err_1,0x0 );
	hLib->CRR = hDcrimi->CRR;

	GIP_free( C_ptr );
	GIP_free( w_sort );
	GIP_free( mask );
	GIP_free( temp );
	if( dist_4==GIP_NULL )	
		GIP_free( inter_0 );

	s = (clock( )-start)*1.0/CLOCKS_PER_SEC;
	hLib->M_spd = s==0.0 ? 0 : (hDcrimi->nz_a+hDcrimi->nz_r)/s;

	G_PRINTF( _T("  total=%d,time=%g,err_1=%g(%d),err_2=%g(%d)"),nEntry,s,err_1,fail_1,err_2,fail_2 );

	return 0x0;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/15/2009	
*/
int _sym_pos_( int r,int c,int ld ){
	int pos = r<=c ? (ld+ld-r+1)*(r)/2+(c-r) : (ld+ld-c+1)*(c)/2+(r-c);

	return pos;
}
/*
	基于GIP_LIB_discriminate_0，采用sym distance

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/15/2009	
*/
int GIP_LIB_discriminate( GIP_LIB *hLib,G_DCRIMI_ *hDcrimi,double *dist_4,int flag )	{
	GU_CTRL_DATA *hUserCtrl=g_hUserCtrl;
	int nEntry=hLib->nEntry,nCls=hLib->nClass,V_len,*C_ptr;
	int i,j,k,no_e0,fail_1,fail_2,D_span=hLib->D_span,D_no,nShift,no,pos;
	int mode=(int)(hLib->user_control[GIRIS_MODE_]);
	G_FVector_ *hfv=hLib->hfv;
	G_FV_TYPE type=hfv->type;
//	int W_dim=hLib->W_dim,*W_ptr=hLib->W_ptr,*W_no=hLib->W_no,W_nz;
	double *inter_0,*inter_1,*intra_0,*intra_1,d,d_min,d_max,s,err_1,err_2,s_D,*w_sort,*dis;
	float *D_inter=hLib->D_inter,*D_intra=hLib->D_intra,*dis_sym;
	G_U8 *temp,*mask;
	clock_t start=clock( );
	GIP_LIB_entry *hEntry,*arrEntry=hLib->arrEntry;
	
	nShift = hfv->nShift;
	nShift = G_DOUBLE2INT(nShift*1.0/180*512/2);
	V_len = hfv->V_len;
	if( hfv->V!=GIP_NULL )		GIP_free(hfv->V);
	hfv->V=GIP_alloc( sizeof(G_U8)*V_len*(2*nShift+1) );
	if( hfv->mask!=GIP_NULL )		GIP_free(hfv->mask);
	hfv->mask=GIP_calloc( sizeof(G_U8),V_len*(2*nShift+1) );
//	W_nz = W_ptr[W_dim];
	C_ptr = GIP_alloc( sizeof(int)*(nCls+1) );
	temp = GIP_alloc( MAX(hfv->unit_size*V_len,nEntry*sizeof(int)) );
	mask = GIP_calloc( hfv->unit_size,V_len*nEntry );
	if( dist_4==GIP_NULL )	{
		inter_0 = GIP_alloc( sizeof(double)*nEntry*4 );
	}else	{
		inter_0 = dist_4;
	}
	inter_1=inter_0+nEntry;	intra_0=inter_1+nEntry;	intra_1=intra_0+nEntry;
	dis_sym=GIP_calloc( sizeof(float),nEntry*(nEntry+1)/2 );

	fail_1=0,				fail_2=0;
	for( i = 0; i < nEntry; i++ )	{
		intra_0[i]=0.0;		intra_1[i]=0.0;
		inter_0[i]=0.0;		inter_1[i]=0.0;
	}
	for( i = 0; i <= D_span; i++ )	{
		D_inter[i]=0.0;		D_intra[i]=0.0;
	}
	w_sort=GIP_alloc( sizeof(double)*nEntry*2 );		dis=w_sort+nEntry;
	for( i = 0; i < nEntry; i++ )		w_sort[i]=0.0;
	s_D=D_span/1.0;		//适用于HAMMING DISTANCE
	GU_SetUserCtrl( hUserCtrl,GIP_NULL,0,0x0 );
	G_PRINTF( _T(" title\t  cls  no   intra_0 intra_1 inter_0 inter_1\t  off") );
	C_ptr[0]=0;
	for( i = 0; i < nCls; i++ )	{
		j=C_ptr[i]+1;
		while( j<nEntry && arrEntry[j].cls==i )	j++;
		C_ptr[i+1] = j;
	//	if( W_dim>69 && i!=69 )	 continue;
		for( j=C_ptr[i]; j < C_ptr[i+1]; j++ )	{
			if( hUserCtrl!=GIP_NULL && hUserCtrl->bTerminate==TRUE )
				break;
			no = j;
			hEntry = hLib->arrEntry+no;
			if( hEntry->status!=G_ENTRY_OK )	{
				fail_1++;			fail_2++;
				continue;
			}
		//	G_fv_shift_( nShift,512/4,2,V_len,hLib->fv_dat+no*V_len,hfv->V,0x0 );
			G_fv_shift_( hfv->type,nShift,V_len,hLib->fv_dat+no*V_len,hfv->V,0x0 );
			G_hamming_dist_( nEntry,dis,hLib->fv_dat,mask,hfv,GIP_NULL,0x0 );	
			for( k = 0; k < nEntry; k++ )	{
				pos = _sym_pos_( j,k,nEntry );
				dis_sym[pos] += dis[k];
			}
		}
		GU_SetUserCtrl( hUserCtrl,GIP_NULL,j+1,0x0 );
	}
	for( i = 0; i < nCls; i++ )	{
		for( j=C_ptr[i]; j < C_ptr[i+1]; j++ )	{
			if( arrEntry[j].status!=G_ENTRY_OK )	
				continue;
			if( C_ptr[i+1]-C_ptr[i]>1 )		{
				d_min=DBL_MAX;			d_max=0.0;
				for( k = C_ptr[i]; k < C_ptr[i+1]; k++ )	{		//intra_class distance
					if( k == j || arrEntry[k].status!=G_ENTRY_OK )		continue;
					pos = _sym_pos_( j,k,nEntry );
					d = dis_sym[pos]/2;
					if( d>d_max )		d_max = d;
					if( d<d_min )		d_min = d;
					D_no=G_DOUBLE2INT( d*s_D );
					D_intra[D_no]++;
				}
				//	D_no=G_DOUBLE2INT( d_min*s_D );
				//	D_intra[D_no]++;
				intra_0[j]=d_min;		intra_1[j]=d_max;
			}
//			w_sort[no]=intra_0[j]==0.0 ? 0.0 : intra_1[j]/intra_0[j];
//			G_PRINTF( _T("no=%d,w=%g,path=%s"),no,w_sort[no],hEntry->sFilePath );
			no_e0=-1;
			if( nCls>1 )	{
				d_min=DBL_MAX;			d_max=0.0;
				for( k = 0; k < nEntry; k++ )	{					//inter_class distance
					if( k>=C_ptr[i] && k<C_ptr[i+1] )		continue;
					if( arrEntry[k].status!=G_ENTRY_OK )		continue;
					pos = _sym_pos_( j,k,nEntry );
					d = dis_sym[pos]/2;
					D_no=G_DOUBLE2INT( d*s_D );
					D_inter[D_no]++;
					if( d>d_max )		d_max = d;
					if( d<d_min )		
					{	no_e0=k;	d_min = d;	}
				}
				inter_0[j]=d_min;		inter_1[j]=d_max;
			}
			w_sort[j]=inter_0[j]==0.0 ? 0.0 : intra_0[j]/inter_0[j];
			w_sort[j]=inter_0[j];
			hEntry->x = w_sort[no];
			s=intra_1[j]==0.0 ? -1: inter_0[j]/intra_1[j];
			if( intra_0[j]>=inter_0[j] )		fail_1++;
			if( intra_1[j]>=inter_0[j] )		fail_2++;
			G_PRINTF( _T("  I_%d\t%4d %3d %7.3g %7.3g %7.3g(%d) %7.3g\t%5.3g"),j,i,j-C_ptr[i],
				intra_0[j],intra_1[j],inter_0[j],no_e0,inter_1[j],s );
		}
	}
#ifdef _GIRIST_DEVELOP_
	GIP_LIB_sort_( hLib,w_sort,(int*)(temp),0x0 );
#endif
	err_1 = fail_1*1.0/nEntry;				err_2 = fail_2*1.0/nEntry;
	G_DCRIMI_set( (void*)hLib,hDcrimi,D_span,D_inter,D_intra,1.0-err_1,0x0 );
	hLib->CRR = hDcrimi->CRR;

	GIP_free( C_ptr );
	GIP_free( w_sort );
	GIP_free( mask );
	GIP_free( temp );
	GIP_free( dis_sym );
	if( dist_4==GIP_NULL )	
		GIP_free( inter_0 );

	s = (clock( )-start)*1.0/CLOCKS_PER_SEC;
	hLib->M_spd = s==0.0 ? 0 : (hDcrimi->nz_a+hDcrimi->nz_r)/s;

	G_PRINTF( _T("  total=%d,time=%g,err_1=%g(%d),err_2=%g(%d)"),nEntry,s,err_1,fail_1,err_2,fail_2 );

	return 0x0;
}

/*
[**]
	读小说 苏北之"恋爱"。。。
	http://blog.sina.com.cn/s/blog_4897e12f01000bzh.html
[**]	
	temp[nSample]

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/15/2009	
*/
void GIP_LIB_sort_( GIP_LIB *hLib,double *s,int *temp,int flag )	{
	int i,j,nEntry = hLib->nEntry,*no,t,isAlloc=temp==GIP_NULL;
	GIP_LIB_entry *hEntry;
	double a;
//	FILE *fp;
	
	G_PRINTF( _T("") );
	G_PRINTF( _T("------ LIB sort ------") );
	if( isAlloc )
		no=GIP_alloc( sizeof(int)*nEntry );
	else 
		no = temp;
	for( i=0; i < nEntry; i++ )
		no[i] = i;
	for( i=0; i < nEntry; i++ )	{
		for( j=i+1; j < nEntry; j++ )	{
			if( s[i]<s[j] ){
				a=s[j];		s[j]=s[i];		s[i]=a;	
				t=no[j];	no[j]=no[i];	no[i]=t;	
			}
		}	
	}

//	fp = _tfopen( _T("lib_sort.info"),_T("w") );
	for( i=0; i < nEntry; i++ )	{
		t = no[i];
		hEntry = hLib->arrEntry+t;
		if( hEntry->status!=G_ENTRY_OK )	
			continue;
		G_PRINTF( _T("no=%d,w=%g,path=%s"),no[i],s[i],hEntry->sFilePath );
//		fprintf( fp,_T("%s"),hEntry->sFilePath );
	}
//	fclose( fp );
	if( isAlloc )
		GIP_free( no );
	G_PRINTF( _T("") );
}
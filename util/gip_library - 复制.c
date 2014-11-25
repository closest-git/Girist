#include<windows.h>
#include <memory.h>
#include <FLOAT.h>
#include <stdio.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/gip_util.h"
#include "../util/gip_ss.h"
#include "../util/gip_fv.h"
#include "GIP_library.h"

static char BLAS_N[1]={'N'},BLAS_T[1]={'T'},BLAS_C[1]={'C'},BLAS_L[1]={'L'},BLAS_R[1]={'R'},BLAS_U[1]={'U'};
static int inc_1=1; 
static GIP_FLOATING one_=1.0,fuyi_=-1.0,zero_=0.0;		 
static char drive[_MAX_DRIVE],dir[_MAX_DIR],fname[_MAX_FNAME],ext[_MAX_EXT];

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/19/2009	
*/
G_FILE_TYPE _File_Type( char *ext )		{
	G_FILE_TYPE type;

	if( strcmp(ext,".jpg")==0x0 )	{
		type = G_IMAGE_JPG;
	}else if( strcmp(ext,".bmp")==0x0 )	{
		type = G_IMAGE_BMP;
	}else if( strcmp(ext,".")==0x0 )	{
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
void GIP_LIB_Append( GIP_LIB *hLIB,char *sPath,int no,int cls,int flag )	{
	GIP_LIB_entry *hEntry,*arrEntry=hLIB->arrEntry;

	if( hLIB->nEntry==hLIB->nMaxEntry )	{
		hLIB->nMaxEntry *= 2;
		arrEntry = GIP_alloc( sizeof(GIP_LIB_entry)*hLIB->nMaxEntry );
		GIP_MEMCOPY( arrEntry,hLIB->arrEntry,sizeof(GIP_LIB_entry)*hLIB->nEntry );
		GIP_free( hLIB->arrEntry );
		hLIB->arrEntry = arrEntry;
	}
/*
		if( (src_image=cvLoadImage(sTestPath,iscolor))==0 )    {
			printf("Image was not loaded.\n");
			return -1;
		}
*/	hEntry=hLIB->arrEntry+hLIB->nEntry;
	GIP_MEMCLEAR( hEntry,sizeof(GIP_LIB_entry) );
	ASSERT( strlen(sPath)<MAX_PATH_LEN );
	strcpy( hEntry->sFilePath,sPath );
//	_splitpath( hEntry->sFilePath, drive, dir, fname, ext );
//	hEntry->sTitle=GIP_alloc( sizeof(char)*(strlen(fname)+1) );
//	strcpy( hEntry->sTitle,fname );
	hEntry->no = no;
	hEntry->cls = cls;

	hLIB->nEntry++;
	ASSERT( hLIB->nEntry<hLIB->nMaxEntry );	
}

/*
	http://www.diybl.com/course/3_program/c++/cppjs/2008331/107762.html

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/19/2009	
*/
void _List_File( GIP_LIB *hLIB,int flag )		{
    HANDLE hSearch;
    WIN32_FIND_DATA data;
	G_FILE_TYPE type;
	int no=0;

    hSearch=FindFirstFile("*",&data);
    do{
        if(data.dwFileAttributes==FILE_ATTRIBUTE_DIRECTORY
           &&strcmp(data.cFileName,".")  && strcmp(data.cFileName,"..") )	{
            SetCurrentDirectory(data.cFileName);
            _List_File( hLIB,flag );
            SetCurrentDirectory("..");
        }
        else if(strcmp(data.cFileName,".") && strcmp(data.cFileName,".."))	{
 			_splitpath( data.cFileName, drive, dir, fname, ext );
            type = _File_Type( ext );
			if( type>G_IMAGE_BEGIN && type<G_IMAGE_END ) {
				char sLine[MAX_LINE];
				GetCurrentDirectory( MAX_LINE, sLine);
				sprintf( sLine,"%s\\%s",sLine,data.cFileName );
				GIP_LIB_Append( hLIB,sLine,no,hLIB->nClass,0x0 );
				no++;
			}
        }
	}while(FindNextFile(hSearch,&data));
	if( no>0 )
 		hLIB->nClass++;

    FindClose(hSearch);
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/25/2008	
*/
int GIP_LIB_init( GIP_LIB *hLIB,char *sLibFile,int flag )	{
	char sLine[MAX_LINE],sCurPath[MAX_LINE];
	FILE*  fp;
	int i,nMaxEntry=1000,len,model,method,nLoop,isGetParameter=0;
	G_FILE_TYPE type;
//	GIP_LIB_entry *hEntry;

	memset( hLIB,0x0,sizeof(GIP_LIB) );
	hLIB->arrEntry = GIP_alloc( sizeof(GIP_LIB_entry)*nMaxEntry );
	hLIB->nMaxEntry = nMaxEntry;
	hLIB->hfv=GIP_alloc( sizeof(G_FVector_) );		

	fp=fopen( sLibFile,"r" );
	if( fp==NULL )	{
		printf( "error at open file: %s", sLibFile );
		return -1;
	}   
/*
	while( fReadLine( fp,sLine ) != 0 )	{
		nMaxEntry++;
	}
	fseek( fp,0x0,0x0 );*/
	GetCurrentDirectory( MAX_LINE,sCurPath );
	while( fReadLine( fp,sLine ) != 0 )	{
		len=(int)(strlen( sLine ));
		if( len==0 )			continue;
		if( sLine[0]=='c' )			
		{	printf( "%s\r\n",sLine );			continue;			}
		if( sLine[0]=='!' )		continue;			//注释行		
		if( sLine[0]=='q' && sLine[1]=='u' && sLine[2]=='i' && sLine[3]=='t' )		{
			break;
		}
		if( sLine[0]=='p' )		{	
			if( sscanf( sLine+1,"%d %d %d",&model,&method,&nLoop )!=3 )
			{	break;		}
			ASSERT( nLoop > 0 );
			isGetParameter=1;
			continue;	
		}
		if( isGetParameter==0 )		
			continue;
		_splitpath( sLine, drive, dir, fname, ext );
		type = _File_Type( ext );
		if( type==G_FILE_DIRECTORY )	{
            SetCurrentDirectory(sLine);
			_List_File( hLIB,0x0 );
		}else if( type>G_IMAGE_BEGIN && type<G_IMAGE_END ) {
			GIP_LIB_Append( hLIB,sLine,0,hLIB->nClass,0x0 );
			hLIB->nClass++;
		}
	}
	fclose( fp );	

    SetCurrentDirectory(sCurPath);

	return 0x0;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/25/2008	
*/
void GIP_LIB_clear( GIP_LIB *hLIB )	{
	int i;
	GIP_LIB_entry *hEntry;
	if( hLIB == GIP_NULL )	return;

	if( hLIB->hfv!=GIP_NULL )	{
		G_FVector_clear( hLIB->hfv );
		GIP_free( hLIB->hfv );				hLIB->hfv=GIP_NULL; 
	}

	if( hLIB->arrEntry != GIP_NULL )	{
		for( i = 0; i < hLIB->nEntry; i++ )	{
			hEntry=hLIB->arrEntry+i;
//			GIP_free( hEntry->sTitle );			hEntry->sTitle=GIP_NULL;
		}
		GIP_free( hLIB->arrEntry );				hLIB->arrEntry=GIP_NULL; 
	}
	hLIB->nEntry = -1;

	if( hLIB->W_ptr != GIP_NULL )	{
		GIP_free( hLIB->W_ptr );			hLIB->W_ptr=GIP_NULL;
	}
	if( hLIB->W_no != GIP_NULL )	{
		GIP_free( hLIB->W_no );			hLIB->W_no=GIP_NULL;
	}
	if( hLIB->W_fv != GIP_NULL )	{
		GIP_free( hLIB->W_fv );			hLIB->W_fv=GIP_NULL;
	}
	

}
/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/20/2009	
*/
void GIP_LIB_scan( GIP_LIB *hLib,int flag )	{
	int i,nEntry=hLib->nEntry,cls=-1,nz=0,W_dim;
	int *W_ptr=hLib->W_ptr,*W_no=hLib->W_no;
	GIP_LIB_entry *hEntry;

	W_ptr=GIP_alloc( sizeof(int)*(hLib->nClass+1) );		hLib->W_ptr=W_ptr;
	W_no=GIP_alloc( sizeof(int)*nEntry );					hLib->W_no=W_no;

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

	hLib->W_dim = W_dim;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/26/2008	
*/
void GIP_LIB_reducing( GIP_LIB *hLib,GIP_FLOATING *sample_0, int flag )	{
	int m,nSample = hLib->nSample,N,i,j,u_block;
	int *W_ptr=hLib->W_ptr,*W_no=hLib->W_no;
	double *W_pca,*W_lda,*mean,*u_mean,*u_s,*u;
	GIP_FLOATING *sample,*S_b,*S_w;

	GIP_LIB_scan( hLib,0x0 );

	N = hLib->fv_n_0;
	m = hLib->nClass;
	mean=GIP_alloc( sizeof(GIP_FLOATING)*N );
	W_pca=GIP_alloc( sizeof(GIP_FLOATING)*N*m );
//	hLib->W=GIP_alloc( sizeof(GIP_FLOATING)*N*m );
	GIP_CalcPCA( nSample,N,m,sample,mean,W_pca );		//get PCA matrix
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
//	dgemm( W_pca,W_lda,hLIB->W )
//	dgemm( sample,hLIB->W,hLIB->fv );
	GIP_free( W_pca );				GIP_free( u );

}

/*
	注意：
		1	非常奇怪的bug
		GIP_LIB包含G_FVector_ fv变量时，读写总是有crash，后改为G_FVector_ *hfv			2/13/2009

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/20/2009	
*/
void GIP_LIB_save( GIP_LIB *hLib,char *sLibPath, int flag )	{
	int nEntry,W_dim=hLib->W_dim,W_nz,N,ret=0;
	FILE *fp;
	G_FVector_ *hfv=hLib->hfv;

	W_nz=hLib->W_ptr[W_dim];
	fp=fopen( sLibPath,"w" );
	if( fp==0x0 )
	{	ret=-1;		goto END;		}
//	if( fwrite( &(hLib->W_dim),sizeof(int),1,fp )!= 1 || fwrite( &(hLib->nEntry),sizeof(int),1,fp )!= 1)
	if( fwrite( hLib,1,sizeof(GIP_LIB),fp )!= sizeof(GIP_LIB) )
	{	ret=-2;		goto END;	}
	if( fwrite( hLib->W_ptr,sizeof(int),W_dim+1,fp )!=W_dim+1 ) 
	{	ret=-3;		goto END;	}
	if( fwrite( hLib->W_no,sizeof(int),W_nz,fp )!=W_nz )
	{	ret=-4;		goto END;	}
	if( fwrite( hLib->arrEntry,sizeof(GIP_LIB_entry),hLib->nEntry,fp )!=hLib->nEntry )
	{	ret=-5;		goto END;	}

	N=hfv->V_len;
	hfv->V=GIP_NULL;
	if( fwrite( hfv,sizeof(G_FVector_),1,fp ) != 1 )
	{	ret=-6;		goto END;	}
	if( fwrite( hLib->W_fv,hfv->unit_size,N*W_nz,fp ) != N*W_nz )
	{	ret=-7;		goto END;	}

END:
	if( fp != 0x0 && fclose( fp )!= 0 )
	{	ret=-8;			}
	if( ret==0 )	{
		printf( "SAVE %s:W_dim=%d,W_nz=%d\r\n\r\n",sLibPath,W_dim,hLib->W_ptr[W_dim] );
	}else	{
		printf( "\t!!!Failed to save %s. err=%d\r\n\r\n",sLibPath,ret );
	}
	if( ret==0 && 1 )	{
		GIP_LIB *hLib;
		hLib=GIP_alloc( sizeof(GIP_LIB) );
		GIP_LIB_load( hLib,sLibPath,0x0 );
		GIP_LIB_clear( hLib );				GIP_free( hLib );
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/20/2009	
*/
void GIP_LIB_load( GIP_LIB *hLib,char *sLibPath, int flag )	{
	int N,W_dim,W_nz,ret=0;
	FILE *fp;
	G_FVector_ *hfv=GIP_NULL;

	memset( hLib,0x0,sizeof(GIP_LIB) );
	fp=fopen( sLibPath,"r" );
	if( fp==0x0 )
	{	ret=-1;		goto END;	}
	
	printf( "LOAD %s,fsize=%d\r\n",sLibPath,0x0 );
//	if( fread( &(hLib->W_dim),sizeof(int),1,fp ) !=1 || fread( &(hLib->nEntry),sizeof(int),1,fp ) !=1 )
	if( fread( hLib,1,sizeof(GIP_LIB),fp ) !=sizeof(GIP_LIB) )
	{	ret=-2;		goto END;	}
	W_dim=hLib->W_dim;				ASSERT( W_dim>0 );
	hLib->W_ptr=GIP_alloc( sizeof(int)*(W_dim+1) );
	if( fread( hLib->W_ptr,sizeof(int),W_dim+1,fp )!=W_dim+1 ) 
	{	ret=-3;		goto END;	}
	W_nz=hLib->W_ptr[W_dim];
//	ASSERT( W_nz>=0 && W_nz<=hLib->nEntry );
	hLib->W_no=GIP_alloc( sizeof(int)*W_nz );
	if( fread( hLib->W_no,sizeof(int),W_nz,fp )!=W_nz )  
	{	ret=-4;		goto END;	}

	hLib->arrEntry=GIP_alloc( sizeof(GIP_LIB_entry)*hLib->nEntry );
	if( fread( hLib->arrEntry,sizeof(GIP_LIB_entry),hLib->nEntry,fp )!=hLib->nEntry )  
	{	ret=-5;		goto END;	}

	hfv=GIP_alloc( sizeof(G_FVector_) );		hLib->hfv = hfv;
	GIP_MEMCLEAR( hfv,sizeof(G_FVector_) );
	if( fread( hfv,sizeof(G_FVector_),1,fp ) != 1 )
	{	ret=-6;		goto END;	}
	hfv->V = GIP_NULL; 
	N=hfv->V_len;
	hLib->W_fv=GIP_alloc( hfv->unit_size*N*W_nz );
	if( fread( hLib->W_fv,hfv->unit_size,N*W_nz,fp )!=N*W_nz )  
	{	ret=-7;		goto END;	}

END:
	if( fp != 0x0 && fclose( fp )!= 0 )
	{	ret=-7;			}
	if( ret==0 )	{
		printf( "LOAD %s\r\n",sLibPath );	
	}else	{
		printf( "\t!!!Failed to load %s. err=%d\r\n\r\n",sLibPath,ret );
	}	
}

/*
[**]
	看狄龙之边缘岁月
	“作兄弟的，有今生，无来世”   嘿嘿。。。。
[**]

	Copyright 2008-present, Grusoft.
	v0.1	cys
		2/7/2009	
*/
void GIP_LIB_discriminate( GIP_LIB *hLib,int flag )	{
	int nSample=hLib->nSample,nClass=hLib->nClass,V_len;
	int i,j,k,no_e0;
	G_FVector_ *hfv=hLib->hfv;
	G_FV_TYPE type=hfv->type;
	int *W_ptr=hLib->W_ptr,*W_no=hLib->W_no;
	double *inter_0,*inter_1,*intra_0,*intra_1,d,d_min,d_max,s;
	char *temp;
	
	V_len = hfv->V_len;
	temp = GIP_alloc( hfv->unit_size*V_len );
	inter_0 = GIP_alloc( sizeof(double)*nSample*4 );
	inter_1=inter_0+nSample;	intra_0=inter_1+nSample;	intra_1=intra_0+nSample;

	for( i = 0; i < nSample; i++ )	{
		intra_0[i]=0.0;		intra_1[i]=0.0;
		inter_0[i]=0.0;		inter_1[i]=0.0;
	}
	printf( " title\t  cls  no   intra_0 intra_1 inter_0 inter_1\t  off\r\n" );
	for( i = 0; i < nClass; i++ )	{
		for( j=W_ptr[i]; j < W_ptr[i+1]; j++ )	{
			if( W_ptr[i+1]-W_ptr[i]>1 )		{
				d_min=DBL_MAX;			d_max=0.0;
				for( k = W_ptr[i]; k < W_ptr[i+1]; k++ )	{		//intra_class distance
					if( k == j )		continue;
					d = G_fv_distance( hLib->W_fv,hfv,j,k,temp );
					if( d>d_max )		d_max = d;
					if( d<d_min )		d_min = d;
				}
				intra_0[j]=d_min;		intra_1[j]=d_max;
			}
			no_e0=-1;
			if( nClass>1 )	{
				d_min=DBL_MAX;			d_max=0.0;
				for( k = 0; k < nSample; k++ )	{					//inter_class distance
					if( k>=W_ptr[i] && k<W_ptr[i+1] )		continue;
					d = G_fv_distance( hLib->W_fv,hfv,j,k,temp );
					if( d>d_max )		d_max = d;
					if( d<d_min )		
					{	no_e0=k;	d_min = d;	}
				}
				inter_0[j]=d_min;		inter_1[j]=d_max;
			}
			s=intra_1[j]==0.0 ? -1: inter_0[j]/intra_1[j];
			printf( "  I_%d\t%4d %3d %7.3g %7.3g %7.3g(%d) %7.3g\t%5.3g\r\n",j+1,i,j-W_ptr[i],
				intra_0[j],intra_1[j],inter_0[j],no_e0+1,inter_1[j],s );
		}
	}

	GIP_free( temp );
	GIP_free( inter_0 );
}

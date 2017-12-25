#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <process.h>    
#include <stddef.h>
#include "gip_thread.h"
#include "../pch/gip_def.h"

#define GIP_MAX_THREADS 8
static int nThread=1;
static void *GIP_u0,*GIP_u1;
static int GIP_code;
static HANDLE event_A[GIP_MAX_THREADS],event_B[GIP_MAX_THREADS],threads[GIP_MAX_THREADS];

GU_CTRL_DATA *g_hUserCtrl=GIP_NULL;
/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/4/2009	
*/
void GIP_driver_core( void* hSolver,void **u_0,void **u_1,int *code,int flag  )	{
	int i;
	if( *code==GIP_DRIVER_START )	{
		for( i = 0; i < nThread; i++ )	
			ResumeThread( threads[i] );
	}else	{
		for( i = 0; i < nThread; i++ )	{
			SetEvent( event_B[i] );
		}
	}
	WaitForMultipleObjects( nThread,event_A,TRUE,INFINITE );

//	*u_0 = GIP_u0,	*u_1=GIP_u1;
	*code = GIP_code;
	if( *code==GIP_DRIVER_FINISH )	{

		for( i = 0; i < nThread; i++ )	{
			SetEvent( event_B[i] );
		}
	}

	return;
}

/*
	v0.1	cys
		8/13/2007
*/
int GIP_driver_switch( void *u_0,void *u_1,int code,int noThread )	{
	int ret=GIP_OK,wait;
	ASSERT( nThread==1 && noThread==0 );
	GIP_u0 = u_0,	GIP_u1 = u_1;
	GIP_code = code;

	SetEvent( event_A[noThread] );
	wait = WaitForSingleObject( event_B[noThread],INFINITE );
	ASSERT( wait != WAIT_FAILED );

	return ret;
}

/*
	v0.1	cys
		8/13/2007
*/
int GIP_driver_init( unsigned ( __stdcall *ThreadProc )( void * ),void* pArguments )	{
	int i,ret=GIP_DRIVER_INIT_FAIL;
	unsigned threadID;

	for( i = 0; i < nThread; i++ )	{	
		event_A[i] = CreateEvent( NULL,FALSE,FALSE,NULL );
		if( event_A[i]==NULL )
			goto END;
		event_B[i] = CreateEvent( NULL,FALSE,FALSE,NULL );
		if( event_B[i]==NULL )
			goto END;

		threads[i] = (HANDLE)_beginthreadex(NULL, 0x0, ThreadProc,pArguments, CREATE_SUSPENDED, &threadID ); 
		if( threads[i]==NULL )
			goto END;
	}
	ret=GIP_OK;
END:
	return ret;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		8/13/2007
*/
int GIP_driver_clear( )	{
	int i,ret=GIP_DRIVER_INIT_FAIL,bRet;
//	unsigned threadID;

	for( i = 0; i < nThread; i++ )	{	
		bRet = CloseHandle( event_A[i] );			ASSERT( bRet != 0 );
		bRet = CloseHandle( event_B[i] );			ASSERT( bRet != 0 );
		bRet = CloseHandle( threads[i] );			ASSERT( bRet != 0 );
	}
	ret=GIP_OK;
END:
	return ret;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
void GU_PROGRESS_clear( GU_CTRL_DATA *hUserCtrl )	{
//	ZeroMemory(this,sizeof(_ProgressThreadData));
	GIP_MEMCLEAR( hUserCtrl,sizeof(GU_CTRL_DATA) );
}

/*
	与PROGRESS_CTRL通讯

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
	v0.2	cys
		7/2/2009			增加switch type
*/
void GU_SetUserCtrl( GU_CTRL_DATA *hUserCtrl,LPCTSTR lpszProgressText,UINT_PTR dwProgressbarPos,int flag ) {
	HWND	hWnd;	
	int type;

	if( hUserCtrl == GIP_NULL )
		return;
	hWnd=hUserCtrl->hThreadWnd;
	type = flag;
	switch( type )	{
	case 0:
		if( lpszProgressText!=GIP_NULL )	{
			if( IsWindow(hWnd) && hUserCtrl->bTerminate == FALSE )
				 SendMessage(hWnd,GUM_PROGRESS_TEXTUPDATE,0,(LPARAM)lpszProgressText);
		}
		if( dwProgressbarPos!=-1 )	{
			if( IsWindow(hWnd) && hUserCtrl->bTerminate == FALSE )	{
				SendMessage(hWnd,GUM_PROGRESS_BARUPDATE,dwProgressbarPos,0);
			}
		}
		break;
	case 1:
		SendMessage(hWnd,GUM_PROGRESS_SETRANGE,dwProgressbarPos,0);
		break;
	default:
		break;
	}
}
		
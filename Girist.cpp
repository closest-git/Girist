// Girist.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "Girist.h"
#include <atlconv.h>
#include <gdiplus.h>
using namespace::Gdiplus;

#include "GiristDlg.h"
#include "GiristDlg.h"
#include ".\util\gip_util.h"
#include ".\util\gip_thread.h"
#include ".\Girist\giris_core.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

double GIRIST_control[GIP_USER_CONTROL];

G_DUMP_CODE g_DumpCode=_DUMP_SBAR;

CGiristDlg* g_pFrame;
CString g_sRootPath,g_sLibPath,g_sDefLib,g_sDumpPath;
CString g_sEXE;

static FILE *g_DumpFile=GIP_NULL;

GIP_LIB g_Lib;
G_FVector_ g_fv;
G_DCRIMI_ g_dcrimi;
GIP_LIB_entry *g_hSource=GIP_NULL,*g_hTarget=GIP_NULL;
// CGiristApp

BEGIN_MESSAGE_MAP(CGiristApp, CWinApp)
	ON_COMMAND(ID_HELP, &CWinApp::OnHelp)
END_MESSAGE_MAP()


// CGiristApp construction

CGiristApp::CGiristApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance

}


// The one and only CGiristApp object

CGiristApp theApp;

/*
	基于log文件的重要性，总是在当前目录输出

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/15/2009
*/
void DumpToFileInit( TCHAR *sTitle,int flag )	{
	CString sPath;

	if( g_DumpFile!=GIP_NULL )
		fclose( g_DumpFile );

	sPath.Format( _T(".\\GIRIST_%s.log"),sTitle );
	g_DumpFile=_tfopen( sPath,_T("wt") ); 
	VERIFY( g_DumpFile!=0x0 );
	g_DumpCode=_DUMP_FILE;

	TCHAR  szTime[60];
	SYSTEMTIME sysTime;
	CTime t = CTime::GetCurrentTime();
	t.GetAsSystemTime(sysTime);
	int nLen;
	nLen = GetTimeFormat(LOCALE_USER_DEFAULT, TIME_NOTIMEMARKER,&sysTime, NULL, szTime, sizeof(szTime));
	_ftprintf( g_DumpFile,_T("file open: %s\r\n"),szTime );

	g_sDumpPath = sPath;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/15/2009
*/
void DumpToFileEnd(  )	{
	TCHAR  szTime[60];
	SYSTEMTIME sysTime;
	CTime t = CTime::GetCurrentTime();
	t.GetAsSystemTime(sysTime);
	int nLen;
	nLen = GetTimeFormat(LOCALE_USER_DEFAULT, TIME_NOTIMEMARKER,&sysTime, NULL, szTime, sizeof(szTime));
	_ftprintf( g_DumpFile,_T("file close: %s\r\n"),szTime );

	fclose( g_DumpFile );
	g_DumpFile=GIP_NULL;
}

// CGiristApp initialization

BOOL CGiristApp::InitInstance()
{
	// InitCommonControlsEx() is required on Windows XP if an application
	// manifest specifies use of ComCtl32.dll version 6 or later to enable
	// visual styles.  Otherwise, any window creation will fail.
	CString str,sControlPath;
	int ret;

	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	// Set this to include all the common control classes you want to use
	// in your application.
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);

	CWinApp::InitInstance();

	AfxEnableControlContainer();

	// Standard initialization
	// If you are not using these features and wish to reduce the size
	// of your final executable, you should remove from the following
	// the specific initialization routines you do not need
	// Change the registry key under which our settings are stored
	// TODO: You should modify this string to be something appropriate
	// such as the name of your compa ny or organization
	SetRegistryKey(_T("Local AppWizard-Generated Applications"));

	GIP_init_heapinfo( );

	g_sEXE.LoadString( IDS_APP_NAME );
#ifdef _GIRIST_DEVELOP_
	g_sRootPath=_T("H:\\Girist\\");
//	g_sRootPath=_T("E:\\Girist\\");
//	g_sLibPath=_T("..\\data\\");
//	g_sDefLib=_T("..\\data\\CASIA_Interval.ilb");
	g_sDefLib=_T("..\\data\\UBIRIS_1_3.ilb");
//	g_sDefLib=_T("..\\data\\MMU_3.ilb");
//	g_sDefLib=_T("..\\data\\test_1.ilb");
#else
#endif

	DumpToFileInit( G_STR2TCHAR(g_sEXE),0x0 );
	sControlPath = _T(".\\iris_control.dat");
	if( (ret=G_IRIS_init_control(GIP_NULL,GIRIST_control,0x0 ))!=GIP_OK )	{
		str.LoadString( IDS_FAIL_INIT_CONTROL );//	
		G_PRINTF( _T("%s. \tERR CODE: %d\r\n"),str,ret );
		return FALSE;
	}
//	int nRet = GIP_LIB_init( &g_Lib,"..\\..\\GIPack\\test\\Girist_test.in",0x0 );	ASSERT(nRet==0x0);
//	GIP_LIB_scan( &g_Lib,0x0 );

	TRY
	{
		CGiristDlg* pMainWnd = new CGiristDlg;
		g_pFrame = pMainWnd;
		m_pMainWnd = pMainWnd;
		if( pMainWnd->Create()==FALSE )
			return FALSE;
		if( __argc==1 )	{
			if( !g_sDefLib.IsEmpty( ) )
				g_pFrame->LoadIrisLib( &g_Lib,g_sDefLib,0x0 );
		//		g_pFrame->OnFileNew( );
		}else	{
		//	if( isTest )	{
		//		g_pFrame->OnToolTest( );
		//	}
		}
	}
	CATCH_ALL(e)
	{
		TRACE0("Failed to create main dialog\n");
		return FALSE;
	}
	END_CATCH_ALL
	/*	CGiristDlg dlg;
	m_pMainWnd = &dlg;
	INT_PTR nResponse = dlg.DoModal();
	if (nResponse == IDOK)	{
	}	else if (nResponse == IDCANCEL)	{
	}*/

	// Since the dialog has been closed, return FALSE so that we exit the
	//  application, rather than start the application's message pump.
	return TRUE;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/2/2009		
*/
int CGiristApp::ExitInstance()	{
	if( g_hTarget!=GIP_NULL )	{
		GIP_free( g_hTarget );
		g_hTarget=GIP_NULL;
	}
	GIP_LIB_clear( &g_Lib );
	G_FVector_clear( &g_fv );

	GIP_exit_heapinfo( );

	DumpToFileEnd( );

	return CWinApp::ExitInstance( );	;
}

IProgressDialog* g_pIDlg=GIP_NULL;
int k=0,IP_max;

/*
	教训：
		1 G_CreateIProgress,G_ReleaseIProgress调用IProgressDialog,但很难与工作线程配合使用

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/4/2009		
*/
int G_CreateIProgress( CString sTitle,CWnd* pwndParent,DWORD nMax,DWORD flag )	{
	HRESULT hr;
	DWORD     m_dwDlgFlags=PROGDLG_NORMAL;
	IProgressDialog* g_pIDlg=GIP_NULL;

	hr = CoCreateInstance ( CLSID_ProgressDialog, NULL, CLSCTX_INPROC_SERVER,
                            IID_IProgressDialog, (void**) &g_pIDlg );
	VERIFY( SUCCEEDED(hr) );
    g_pIDlg->SetTitle ( sTitle );
    hr = g_pIDlg->StartProgressDialog ( pwndParent->GetSafeHwnd(),
                                    NULL,
                                    m_dwDlgFlags | PROGDLG_MODAL,
                                    NULL );
	VERIFY( SUCCEEDED(hr) );

	IP_max = nMax;
	g_pIDlg->SetProgress ( 0, nMax );
    g_pIDlg->Timer ( PDTIMER_RESET, NULL );
	::Sleep(1500);

	return 0x0;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/4/2009		
*/
void G_ReleaseIProgress( IProgressDialog* g_pIDlg )	{
    g_pIDlg->StopProgressDialog();
    g_pIDlg->Release();
}

/*
	输出到不同的对像

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/2/2009
*/
extern "C" void G_PRINTF( const TCHAR *format,... )	{
	va_list args;
	int len;
	TCHAR * buffer;

	va_start( args, format );
	len = _vscwprintf( format, args ) // _vscprintf doesn't count
                               + 1; // terminating '\0'
	buffer = (TCHAR*)malloc( len * sizeof(TCHAR) );
	vswprintf( buffer, format, args );
//	TRACE( CString(buffer) );

	if( g_DumpCode==_DUMP_USER_CTRL )	{
		if( g_DumpFile!=GIP_NULL )	{
			_ftprintf( g_DumpFile,_T("%s\n"),buffer );
			fflush( g_DumpFile );
		}
		GU_SetUserCtrl( g_hUserCtrl,buffer,-1,0x0 );
	}else if( g_DumpCode==_DUMP_FILE )	{
		VERIFY( g_DumpFile!=GIP_NULL );
			_ftprintf( g_DumpFile,_T("%s\n"),buffer );
			fflush( g_DumpFile );
	} /*else if( g_pIDlg!=GIP_NULL )	{
		::Sleep(500);
		g_pIDlg->SetProgress ( ++k, IP_max );
	}*/else	{
//		g_pFrame->SetSBInfo( ID_SEPARATORID_INDICATOR_CORE,str );
	}

	free( buffer );
	return ;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/2/2009

extern "C" void G_PRINTF( void *hUserCtrl,const TCHAR *format,... )	{
	va_list args;
	int len;
	TCHAR * buffer;

	va_start( args, format );
	len = _vscwprintf( format, args )+ 1; //  _vscprintf doesn't count terminating '\0'                               
	buffer = (TCHAR*)malloc( len * sizeof(TCHAR) );
	vswprintf( buffer, format, args );
	TRACE( CString(buffer) );
	GU_SetUserCtrl( (GU_CTRL_DATA*)hUserCtrl,buffer,-1,0x0 );

	free( buffer );
}
*/


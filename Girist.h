// Girist.h : main header file for the PROJECT_NAME application
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h"		// main symbols
#include "./pch/gip_def.h"		
#include "./util/gip_fv.h"		
#include "./util/gip_library.h"		

// CGiristApp:
// See Girist.cpp for the implementation of this class
//
typedef enum{
	_DUMP_USER_CTRL=0x100,_DUMP_SBAR,_DUMP_FILE,
}G_DUMP_CODE;
extern G_DUMP_CODE g_DumpCode;

class CGiristApp : public CWinApp
{
public:
	CGiristApp();

// Overrides
	public:
	virtual BOOL InitInstance();
	virtual int ExitInstance(); // return app exit code

// Implementation

	DECLARE_MESSAGE_MAP()
};

extern CGiristApp theApp;
extern CString g_sRootPath,g_sLibPath,g_sDumpPath;
extern CString g_sEXE;

extern double GIRIST_control[GIP_USER_CONTROL];
extern GIP_LIB g_Lib;
extern G_FVector_ g_fv;
extern G_DCRIMI_ g_dcrimi;
extern GIP_LIB_entry *g_hSource,*g_hTarget;

#define G_STR2TCHAR(str)	(TCHAR*)(LPTSTR)(LPCTSTR)(str)

int G_CreateIProgress( CString sTitle,CWnd* pwndParent,DWORD nMax,DWORD flag );
void DumpToFileInit( TCHAR *sTitle,int flag );

/*
	v0.1	cys
		3/21/2006
*/
CString inline GetFileDir( CString sPath )	{
	TCHAR drive[_MAX_DRIVE];
	TCHAR dir[_MAX_DIR];
	TCHAR fname[_MAX_FNAME];
	TCHAR ext[_MAX_EXT];
	CString sDir;

	_wsplitpath_s( sPath.GetBuffer(sPath.GetLength()), drive,_MAX_DRIVE,dir,_MAX_DIR,fname,_MAX_FNAME,ext,_MAX_EXT );
	sDir = CString(drive)+CString(dir);
//	sDir = dir;
	return sDir;
}

/*
	v0.1	cys
		3/21/2006
*/
CString inline GetFileTitle( CString sPath )	{
	TCHAR drive[_MAX_DRIVE];
	TCHAR dir[_MAX_DIR];
	TCHAR fname[_MAX_FNAME];
	TCHAR ext[_MAX_EXT];

	_wsplitpath_s( sPath.GetBuffer(sPath.GetLength()), drive,_MAX_DRIVE,dir,_MAX_DIR,fname,_MAX_FNAME,ext,_MAX_EXT );
	return fname;
}

/*
	v0.1	cys
		3/21/2006
*/
CString inline GetFileExt( CString sPath )	{
	TCHAR drive[_MAX_DRIVE];
	TCHAR dir[_MAX_DIR];
	TCHAR fname[_MAX_FNAME];
	TCHAR ext[_MAX_EXT];

	_wsplitpath_s( sPath.GetBuffer(sPath.GetLength()), drive,_MAX_DRIVE,dir,_MAX_DIR,fname,_MAX_FNAME,ext,_MAX_EXT );
	return ext+1;
}

/*
	v0.1	cys
		5/17/2006
*/
BOOL inline RunEXE( LPTSTR  szCmd  )	{
	STARTUPINFO StartupInfo;
	memset( &StartupInfo,0,sizeof( STARTUPINFO ) );
	StartupInfo.cb = sizeof( STARTUPINFO );
	PROCESS_INFORMATION ProcessInfo;

	BOOL bRet = ::CreateProcess( NULL,szCmd,NULL,NULL,FALSE,
		0,NULL,NULL,&StartupInfo,&ProcessInfo );

	ASSERT( bRet );

	return bRet;
}


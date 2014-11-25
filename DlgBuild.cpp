// DlgBuild.cpp : implementation file
//

#include "stdafx.h"
#include <process.h>    
#include "Girist.h"
#include "DlgBuild.h"
#include "DlgDistri.h"
#include ".\pch\gip_def.h"
#include ".\Girist\giris_core.h"
#include ".\util\gip_thread.h"


// CDlgBuild dialog

IMPLEMENT_DYNAMIC(CDlgBuild, CDialog)

CDlgBuild::CDlgBuild(GIP_LIB *hLib,CWnd* pParent /*=NULL*/)
	: CDialog(CDlgBuild::IDD, pParent)
{
	m_hLib = hLib;
	GU_PROGRESS_clear( &m_Progress );

}

CDlgBuild::~CDlgBuild()
{
}

void CDlgBuild::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_PROGRESS1, m_cProgress);
	DDX_Control(pDX, IDC_LIST2, m_lbInfo);
	DDX_Control(pDX, IDC_BT_SWITCH, m_btSwitch);
	DDX_Control(pDX, IDC_BT_DISTRI, m_btDistri);
	DDX_Control(pDX, IDCANCEL, m_btCancel);
	DDX_Control(pDX, IDOK, m_btOK);
}

#define BUILD_STAGE_NUM 3
static int stage=0;
/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
static unsigned __stdcall  _ThreadProc(LPVOID lpThreadParameter)	{
    GU_CTRL_DATA *hUserCtrl = (GU_CTRL_DATA*) lpThreadParameter;
	INT_PTR nResult=-1;
	GIP_LIB *hLib=(GIP_LIB *)(hUserCtrl->pUserParam);
	int nEntry=hLib->nEntry;
	clock_t start=clock( );
	double t;

	g_hUserCtrl = hUserCtrl;
	hUserCtrl->bAlive = TRUE;
	if(hUserCtrl->bTerminate == TRUE)
		goto TerminateThreadProc;

	stage = 1;				
	G_PRINTF( _T("========== Build started: ==========") );
	G_PRINTF( _T("1>------ IRIS LOCALIZATION & NORMALIZATION ------") );
	start=clock( );	
	nResult = G_IRIS_location( hLib,0x0 );
	t = (clock()-start)*1.0/CLOCKS_PER_SEC;
	G_PRINTF( _T("1>------ TIME %g------"),t );
	if( nResult != 0x0 || hUserCtrl->bTerminate == TRUE )
		goto TerminateThreadProc;

	stage = 2;
	G_PRINTF( _T("2>------ FEATURE EXTRACTION & CODING ------") );
	start=clock( );	
	nResult = G_IRIS_coding( hLib,0x0 );
	t = (clock()-start)*1.0/CLOCKS_PER_SEC;
	G_PRINTF( _T("2>------ TIME %g------"),t );
	if( nResult != 0x0 || hUserCtrl->bTerminate == TRUE )
		goto TerminateThreadProc;

	stage = 3;
	G_PRINTF( _T("3>------ INTRA-CLASS & INTER-CLASS COMPARISON ------") );
	start=clock( );	
	nResult = GIP_LIB_discriminate( hLib,&g_dcrimi,GIP_NULL,0x0 );
	t = (clock()-start)*1.0/CLOCKS_PER_SEC;
	g_dcrimi.time = t;
	G_PRINTF( _T("3>------ TIME %g------"),t );
	if( nResult != 0x0 || hUserCtrl->bTerminate == TRUE )
		goto TerminateThreadProc;	/**/

	G_PRINTF(  _T("========== Build: %d succeeded, %d failed, %d skipped =========="),
		hLib->nSample,nEntry-hLib->nSample,0 );

TerminateThreadProc:
	if(hUserCtrl->bTerminate == FALSE)
		::PostMessage(hUserCtrl->hThreadWnd,GUM_PROGRESS_THREADCOMPLETED,MAKEWPARAM(nResult,0),0);
	else
		::PostMessage(hUserCtrl->hThreadWnd,GUM_CANCEL_PROGRESSTHREAD,MAKEWPARAM(nResult,0),0);

	hUserCtrl->bAlive = FALSE;
	return 0;    
}

/*
	bAlive=TRUE:	工作线程在运行，只可以按m_btSwitch

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
void CDlgBuild::CtrlOnState( BOOL bAlive )	{
	CString	sBt = bAlive==TRUE ? _T("STOP") : _T("Rebuild");
	m_btSwitch.SetWindowTextW( sBt );

//	m_btCancel.EnableWindow( !bAlive );
	m_btOK.EnableWindow( !bAlive );
	m_btDistri.EnableWindow( !bAlive );
	m_lbInfo.EnableWindow( !bAlive );
}

/*
	注意：
		1 存在无法分配线程的可能

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
int CDlgBuild::StartThread( )	{
	unsigned threadID;
	CString str;

	GU_PROGRESS_clear( &m_Progress );
	m_Progress.hThreadWnd = this->GetSafeHwnd( );
	m_Progress.pUserParam = (LPVOID)(m_hLib);
	m_hThread = (HANDLE)_beginthreadex(NULL, 0x0, _ThreadProc,&(m_Progress), CREATE_SUSPENDED, &threadID );
	if( m_hThread==NULL )	{		//存在失败的可能
		str.LoadStringW( IDS_THREAD_CREATE_FAIL );
		str.LoadStringW( IDS_GIRIST_RESTART );
		MessageBox( str,g_sEXE,MB_OK );	
		return -1;
	}
	m_Progress.bAlive = TRUE;
	CtrlOnState( m_Progress.bAlive );
	m_lbInfo.ResetContent( );

	ResumeThread( m_hThread );

	return 0;
}


/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/6/2009
*/
BOOL CDlgBuild::OnInitDialog( ) {
	CDialog::OnInitDialog();

	ASSERT( m_hLib->nEntry>0 );
	m_cProgress.SetRange( 0,m_hLib->nEntry*3 );

	if( StartThread( )!= 0x0 )
		return FALSE;		

	return TRUE;  // return TRUE  unless you set the focus to a control
}	

BEGIN_MESSAGE_MAP(CDlgBuild, CDialog)
	ON_BN_CLICKED(IDC_BT_SWITCH, &CDlgBuild::OnBnClickedBtSwitch)
	ON_MESSAGE( GUM_PROGRESS_BARUPDATE, OnUserBarUpdate) 
	ON_MESSAGE( GUM_PROGRESS_TEXTUPDATE, OnUserTextUpdate) 
	ON_MESSAGE( GUM_PROGRESS_THREADCOMPLETED, OnThreadCompleted) 
	ON_BN_CLICKED(IDC_BT_DISTRI, &CDlgBuild::OnBnClickedBtDistri)
END_MESSAGE_MAP()

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
LRESULT CDlgBuild::OnUserBarUpdate(UINT wParam, LONG lParam)	{
	int pos = wParam;
	m_cProgress.SetPos( pos+m_hLib->nEntry*(stage-1) );
	return 0L; 
} 

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
LRESULT CDlgBuild::OnUserTextUpdate(UINT wParam, LONG lParam)	{
	int cur = m_lbInfo.AddString( (LPCTSTR)lParam );
	m_lbInfo.SetTopIndex( cur );

	return 0L; 
} 

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
LRESULT CDlgBuild::OnThreadCompleted(UINT wParam, LONG lParam)	{
	CString str;
	int delay=500;

	while( m_Progress.bAlive==TRUE )	{
		Sleep(delay);
	}

	CtrlOnState( m_Progress.bAlive );
	return 0L; 
} 
 
/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/6/2009
*/
void CDlgBuild::OnCancel()	{
	CString str;
	int delay=500;

	if( m_Progress.bAlive )	{
		m_Progress.bTerminate=TRUE;
		while( m_Progress.bAlive==TRUE )	{
			Sleep(delay);
		}
		str.LoadStringW( IDS_BUILD_CANCEL );
		m_lbInfo.AddString( str );
		CtrlOnState( m_Progress.bAlive );
	}

	CDialog::OnCancel( );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/6/2009
*/
void CDlgBuild::OnOK()	{
	m_hLib->isUpdate = 1;

	CDialog::OnOK( );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/7/2009
*/
void CDlgBuild::OnBnClickedBtSwitch()
{
	// TODO: Add your control notification handler code here
	CString str;
	int delay=500;

	if( m_Progress.bAlive )	{
		m_Progress.bTerminate=TRUE;
		while( m_Progress.bAlive==TRUE )	{
			Sleep(delay);
		}
		str.LoadStringW( IDS_BUILD_CANCEL );
		m_lbInfo.AddString( str );
		CtrlOnState( m_Progress.bAlive );
	}else	{
		StartThread( );
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/12/2009
*/
void CDlgBuild::OnBnClickedBtDistri()
{
	// TODO: Add your control notification handler code here
	CDlgDistri dlg( m_hLib );

	dlg.DoModal( );
}

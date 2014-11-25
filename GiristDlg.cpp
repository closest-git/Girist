// GiristDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Girist.h"
#include "GiristDlg.h"
#include "DlgProgress.h"
#include "DlgBuild.h"
#include "DlgDistri.h"
#include ".\pch\gip_def.h"
#include ".\Girist\giris_core.h"
#include ".\util\gip_library.h"
#include ".\util\gip_util.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

static UINT auIDStatusBar[] = {
   ID_SEPARATOR,
   ID_INDICATOR_CORE,
   ID_INDICATOR_TIME
};
/*
	[ISSUE-NEEDED]
	1 需由工作线程处理G_IRIS_build，并由CDlgProgress显示信息.参见UPDialog.zip

*/

// CGiristDlg dialog

void CGiristDlg::PostNcDestroy()
{

	delete this;
}

/////////////////////////////////////////////////////////////////////////////
// CDialog::OnCancel
//      Modeless dialogs must not call EndDialog(), as CDialog::OnCancel()
//      does.  All we need to do is destroy the window.

void CGiristDlg::OnCancel()
{
	SaveCurLib( &g_Lib,_PROMPT_SAVE_ );
//	if( OnExit( ) == TRUE )
		DestroyWindow();
}

/////////////////////////////////////////////////////////////////////////////
// CDialog::OnOK
//      Modeless dialogs must not call EndDialog(), as CDialog::OnOK() does
//      We retrieve data from the controls using DDX, then destroy the
//      window.

void CGiristDlg::OnOK()
{
	if (!UpdateData(TRUE))
	{
		TRACE0("UpdateData failed -- modeless dialog terminate\n");
		return;
	}
	SaveCurLib( &g_Lib,_PROMPT_SAVE_ );
//	if( OnExit( ) == TRUE )
		DestroyWindow();
}



CGiristDlg::CGiristDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CGiristDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
	isTargetCheck=FALSE;
}

/*
	v0.1	cys
		3/2/2009		
*/
CGiristDlg::~CGiristDlg()
{
//	Gdiplus::GdiplusShutdown(m_gdiplusToken);
}

/*
	v0.1	cys
		3/2/2009		
*/
BOOL CGiristDlg::Create()
{
//	Gdiplus::GdiplusStartupInput gdiplusStartupInput;
//	Gdiplus::GdiplusStartup(&m_gdiplusToken, &gdiplusStartupInput, NULL);

	m_hIcon = AfxGetApp()->LoadIcon(CGiristDlg::IDD);
	return CDialog::Create(CGiristDlg::IDD);
}

void CGiristDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

	DDX_Control(pDX, IDC_TAB_SOURCE, m_tabSource);
	DDX_Control(pDX, IDC_TAB_TARGET, m_tabTarget);
	DDX_Control(pDX, IDC_TAB_PROJ, m_tabProj);
}

BEGIN_MESSAGE_MAP(CGiristDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()

	ON_COMMAND(ID_IRIS_AUTHEN, OnIrisAuthen )
	ON_COMMAND(ID_IRIS_IDENTIFY, OnIrisIdentify )
	ON_UPDATE_COMMAND_UI(ID_IRIS_IDENTIFY,OnUpdateIrisIdentify)

	ON_COMMAND(ID_LIB_CREATE, OnLibCreate )
	ON_COMMAND(ID_LIB_OPEN, OnLibOpen )
	ON_COMMAND(ID_LIB_CLOSE, OnLibClose )
	ON_COMMAND(ID_FILE_OPEN, OnFileOpen )
	ON_COMMAND(ID_LIB_BUILD, OnLibBuild )
	ON_COMMAND(ID_LIB_SAVE, OnLibSave )
	ON_COMMAND(ID_LIB_DISTRI, OnLibDistri )
	ON_COMMAND(ID_LIB_INFO,OnLibInfo )
//	ON_COMMAND(ID_TOOL_OPTION, OnToolOption )
	ON_COMMAND(ID_APP_ABOUT, OnAppAbout )
	ON_UPDATE_COMMAND_UI(ID_INDICATOR_TIME, OnUpdateTime)
//	ON_CBN_SELCHANGE(ID_CB_MODE, OnSelChangeMode )
//	ON_CBN_SELCHANGE(ID_CB_DEPTH, OnSelChangeDepth )
//	ON_COMMAND_RANGE(ID_CHESS_BEGIN,ID_CHESS_END,OnChessRange)
	ON_UPDATE_COMMAND_UI(ID_IRIS_AUTHEN,OnUpdateIrisAuthen )
	ON_UPDATE_COMMAND_UI(ID_LIB_DISTRI,OnUpdateLibDistri)
	ON_UPDATE_COMMAND_UI(ID_LIB_BUILD,OnUpdateLibBuild)
	ON_UPDATE_COMMAND_UI(ID_LIB_INFO,OnUpdateLibInfo)
	ON_UPDATE_COMMAND_UI(ID_LIB_SAVE,OnUpdateLibSave)
	ON_UPDATE_COMMAND_UI(ID_LIB_CLOSE,OnUpdateLibClose)
//	ON_UPDATE_COMMAND_UI_RANGE(ID_CHESS_BEGIN,ID_CHESS_END,OnUpdateChessRange)

	//}}AFX_MSG_MAP
	ON_NOTIFY(NM_CLICK, IDC_TAB_TARGET, &CGiristDlg::OnNMClickTabTarget)
END_MESSAGE_MAP()

static CString Mode_String[]	= { _T("PLAY"),_T("BOOK"),_T("NEXT") };

// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnStnDblclickBmpLogo();
};

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	ON_STN_DBLCLK(IDC_BMP_LOGO, &CAboutDlg::OnStnDblclickBmpLogo)
END_MESSAGE_MAP()

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

//	HBITMAP   bitmap=(HBITMAP)::LoadImage(AfxGetInstanceHandle(),MAKEINTRESOURCE(IDB_GRUSOFT),IMAGE_BITMAP,0,0,LR_CREATEDIBSECTION);//得到bmp句柄   
//	((CStatic*)GetDlgItem(IDC_BMP_LOGO))->SetBitmap( bitmap );
}

/*
	v0.1	cys
		4/13/2009
*/
void CAboutDlg::OnStnDblclickBmpLogo()
{
	// TODO: Add your control notification handler code here
	ShellExecute(NULL,NULL,_T("http://www.grusoft.com/girist.htm"),NULL,NULL,SW_SHOWNORMAL); 
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/6/2009
*/
void CGiristDlg::OnAppAbout()
{
	CAboutDlg aboutDlg;
	aboutDlg.DoModal();
}
/*
	Add status bar,toolbar and the boardcontrol

	v0.1	cys
		2/28/2009
*/
int CGiristDlg::CreateControls( )	{
	// Create the board control 
	m_tiS.type = TAB_IRIS_SOURCE;
	m_tiS.Create (IDD_TAB_IRIS, this);
	m_tabSource.AddSSLPage (_T("SOURCE"), 0, &m_tiS);
	m_tabSource.ActivateSSLPage( 0 );

	m_tiT.type = TAB_IRIS_TARGET;
	m_tiT.Create (IDD_TAB_IRIS, this);
	m_tabTarget.AddSSLPage (_T("TARGET"), 0, &m_tiT);
	m_tabTarget.ActivateSSLPage( 0 );

	m_tabTree.Create (IDD_TAB_TREE, this);
	m_tabWarn.Create (IDD_TAB_TREE, this);
	m_tabFail.Create (IDD_TAB_TREE, this);
#ifdef _GIRIST_DEVELOP_
	m_tabSort.Create (IDD_TAB_TREE, this);
	m_tabSort.m_tree.m_isSort = 1;
#endif

	m_tabProj.AddSSLPage (_T("Library"), 1, &m_tabTree );
	m_tabProj.AddSSLPage (_T("Warning"), 1, &m_tabWarn );
	m_tabProj.AddSSLPage (_T("Failed"), 1, &m_tabFail );
#ifdef _GIRIST_DEVELOP_
	m_tabProj.AddSSLPage (_T("Sorted"), 1, &m_tabSort );
#endif
	m_tabProj.ActivateSSLPage( 0 );

	if (m_statusBar.Create(this))
	{
		m_statusBar.SetIndicators(auIDStatusBar, sizeof(auIDStatusBar)/sizeof(UINT) );
		OnSetMessageString(AFX_IDS_IDLEMESSAGE);

		// Make a sunken or recessed border around the first pane
		m_statusBar.SetPaneInfo(0, m_statusBar.GetItemID(0),
			SBPS_STRETCH, NULL );
	}

	// Create toolbar at the top of the dialog window
	if ( !m_toolBar.Create(this) || !m_toolBar.LoadToolBar(IDR_TOOLBAR) )	{
		TRACE0("Failed to create toolbar\n");
		return -1;      // fail to create
	}
/*	int index = 0;
	CRect rect;
	while( m_toolBar.GetItemID( index ) != ID_SEL_MODE )
		index++;
	m_toolBar.SetButtonInfo( index,ID_SEL_MODE,TBBS_SEPARATOR,80 );
	m_toolBar.GetItemRect( index,&rect );
	rect.bottom += 160;
	if( !m_cbModes.Create( WS_CHILD | WS_VISIBLE |CBS_AUTOHSCROLL 
		|CBS_DROPDOWNLIST |CBS_HASSTRINGS,rect,&m_toolBar,ID_CB_MODE ) )	{
		TRACE( "Failed to create pattern combo-box\n" );
		return -1;
	}
	for( int i = 0; i < 3; i++ )
		m_cbModes.AddString( Mode_String[i] );
//	m_cbModes.SetCurSel( m_ctrlBoard.mode );
	m_cbModes.EnableWindow( FALSE );
*/
	m_toolBar.SetBarStyle(m_toolBar.GetBarStyle() |
		CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC);

	ResizeControls( );

	return 0;
}

/////////////////////////////////////////////////////////////////////////////
// CRenStarDlg::OnSetMessageString
//      OnSetMessageString updates the status bar text.
//
//      This code is based on CFrameWnd::OnSetMessageString.  We assume
//      a string ID is always passed in wParam.

LRESULT CGiristDlg::OnSetMessageString(WPARAM wParam, LPARAM /*lParam*/)
{
	UINT    nIDMsg = (UINT)wParam;
	CString strMsg;

	if (nIDMsg)
	{
		if (strMsg.LoadString(nIDMsg) != 0)
		{
			CString strStatusText;
			AfxExtractSubString(strStatusText, strMsg, 0, '\n');
			m_statusBar.SetWindowText(strStatusText);
		}
		else
			TRACE1("Warning: no message line prompt for ID %x%04X\n", nIDMsg);
	}

	UINT nIDLast     = m_nIDLastMessage;
	m_nIDLastMessage = nIDMsg;
	m_nIDTracking    = nIDMsg;
	return nIDLast;

}


/*
	CGiristDlg message handlers
	
	v0.1	cys
		2/28/2009
*/
BOOL CGiristDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	if( CreateControls( )==-1 )
		return FALSE;

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here
	m_tabTree.m_tree.OnUpdate( G_ENTRY_ALL );
	m_tabWarn.m_tree.OnUpdate( G_ENTRY_WARNING );
	m_tabFail.m_tree.OnUpdate( G_ENTRY_FAILED );
#ifdef _GIRIST_DEVELOP_
	m_tabSort.m_tree.OnUpdate( G_ENTRY_ALL );
#endif

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CGiristDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
//		CAboutDlg dlgAbout;
//		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CGiristDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CGiristDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

/*
	v0.1	cys
		3/16/2006
*/
void CGiristDlg::ResizeControls( )	{
	// We need to resize the dialog to make room for control bars.
	// First, figure out how big the control bars are.
	CRect rcClientStart,rcTree,rcBoard;
	CRect rcClientNow;
	GetClientRect(rcClientStart);
//	m_TAB.GetWindowRect(rcTree);
	ScreenToClient(rcTree);
	rcBoard = rcClientStart;
	rcBoard.right = rcTree.left;
//	::MoveWindow(m_ctrlBoard.m_hWnd, rcBoard.left, rcBoard.top, rcBoard.right, rcBoard.bottom, TRUE);
	RepositionBars(AFX_IDW_CONTROLBAR_FIRST, AFX_IDW_CONTROLBAR_LAST,0, reposQuery, rcClientNow);

	// Now move all the controls so they are in the same relative
	// position within the remaining client area as they would be
	// with no control bars.
	CPoint ptOffset(rcClientNow.left - rcClientStart.left,rcClientNow.top - rcClientStart.top);

	CRect  rcChild;
	CWnd* pwndChild = GetWindow(GW_CHILD);
	while (pwndChild)
	{
		pwndChild->GetWindowRect(rcChild);
		ScreenToClient(rcChild);
		rcChild.OffsetRect(ptOffset);
		pwndChild->MoveWindow(rcChild, FALSE);
		pwndChild = pwndChild->GetNextWindow();
	}

	// Adjust the dialog window dimensions
	CRect rcWindow;
	GetWindowRect(rcWindow);
	rcWindow.right += rcClientStart.Width() - rcClientNow.Width();
	rcWindow.bottom += rcClientStart.Height() - rcClientNow.Height();
	MoveWindow(rcWindow, FALSE);

	// And position the control bars
	RepositionBars(AFX_IDW_CONTROLBAR_FIRST, AFX_IDW_CONTROLBAR_LAST, 0);	
	// Resize the button control proportionally to dialogs size
}

/*
	v0.1	cys
		3/2/2009
*/
void CGiristDlg::OnUpdateTime(CCmdUI* pCmdUI)
{
	CTime t = CTime::GetCurrentTime();
	TCHAR  szTime[6];
	SYSTEMTIME sysTime;
	t.GetAsSystemTime(sysTime);

	// Localize to user locale, strip seconds and time marker
	int nLen;
	nLen = GetTimeFormat(LOCALE_USER_DEFAULT, TIME_NOSECONDS|TIME_NOTIMEMARKER,
			&sysTime, NULL, szTime, sizeof(szTime));
	SetSBInfo( ID_INDICATOR_TIME,szTime );
	pCmdUI->Enable();
}

/*
	v0.1	cys
		4/9/2009
*/
void CGiristDlg::OnUpdateLibDistri( CCmdUI* pCmdUI )	{
	int mode=(int)(g_Lib.user_control[GIRIS_MODE_]);
	BOOL isOn;
	isOn = g_dcrimi.nz_a>0 || g_dcrimi.nz_r>0 ;

	pCmdUI->Enable( isOn && mode!=GIRIS_MODE_AUTHEN );
}

/*
	v0.1	cys
		4/9/2009
*/
void CGiristDlg::OnUpdateIrisAuthen( CCmdUI* pCmdUI )	{
	int mode=(int)(g_Lib.user_control[GIRIS_MODE_]);

	pCmdUI->SetCheck( mode==GIRIS_MODE_AUTHEN );
}

/*
	v0.1	cys
		4/9/2009
*/
void CGiristDlg::OnUpdateLibBuild( CCmdUI* pCmdUI )	{
	int mode=(int)(g_Lib.user_control[GIRIS_MODE_]);
	pCmdUI->Enable( g_Lib.nEntry>0 && mode!=GIRIS_MODE_AUTHEN );
}

/*
	v0.1	cys
		4/13/2009
*/
void CGiristDlg::OnUpdateLibInfo( CCmdUI* pCmdUI )	{
	pCmdUI->Enable( g_Lib.nEntry>0 );
}

/*
	v0.1	cys
		4/13/2009
*/
void CGiristDlg::OnUpdateLibSave( CCmdUI* pCmdUI )	{
	pCmdUI->Enable( g_Lib.nEntry>0 && g_Lib.isUpdate!=0 );
}
/*
	v0.1	cys
		4/19/2009
*/
void CGiristDlg::OnUpdateLibClose( CCmdUI* pCmdUI )	{
	pCmdUI->Enable( g_Lib.nEntry>0 );
}
/*
	v0.1	cys
		4/13/2009
*/
void CGiristDlg::OnUpdateIrisIdentify( CCmdUI* pCmdUI )	{
	int mode=(int)(g_Lib.user_control[GIRIS_MODE_]);
	pCmdUI->Enable( g_Lib.nEntry>0 && mode!=GIRIS_MODE_AUTHEN );
}

/*
	v0.1	cys
		3/1/2009
*/
void CGiristDlg::OnFileOpen( ) 
{
	// TODO: Add your command handler code here
	CFileDialog dlg( TRUE,NULL,NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
		_T("Iris Image(*.bmp)|*.bmp|")
		_T("All Files(*.*)|*.*|")
#ifdef _GRUS_RENSTAR_EXPERT_
		_T("RenStar Book(*.book)|*.book|")
		_T("RenLib File(*.lib)|*.lib|")
#endif
	);
	if( dlg.DoModal( ) != IDOK )	
		return;

	CString sPath = dlg.GetPathName( ),sExt=GetFileExt( sPath ),sLibPath; 
	if( sExt=="jpg" )	{
	}else if( sExt=="ilb" )	{
	}else	{
		m_tabTree.m_tree.OnUpdate( G_ENTRY_ALL );
	}

	SetWindowText( sPath );

}

/*
	for SHBrowseForFolder
	http://topic.csdn.net/t/20050507/10/3986826.html

	v0.1	cys
		3/5/2009
*/
int CALLBACK   _SHBFF_call_back_(HWND   hwnd,   UINT   msg,   LPARAM   lp,   LPARAM   pData)   { 
	if( msg==BFFM_INITIALIZED )   	{   
		  ::SendMessage(hwnd,BFFM_SETSELECTION,TRUE,(LPARAM)(G_STR2TCHAR(g_sRootPath)) );   
	}  
#ifdef _GIRIST_DEVELOP_
#endif
	return   0;   
}   

/*
	Create New Lib from Directory

	
	ISSUE-NEEDED:
		MB_YESNOCANCEL改为SAVE_NOSVAE_CANCEL

	v0.1	cys
		3/5/2009
*/
void CGiristDlg::OnLibCreate( ) {
	// TODO: Add your command handler code here
	TCHAR   szDir[MAX_PATH],szCurOld[MAX_PATH];   
	BROWSEINFO   bi;   
	ITEMIDLIST   *pidl; 
	int isCreate=0;
	CString strDlgTitle=_T("select the path of iris files:"),sLibPath,str;   

	if( g_Lib.isUpdate > 0 )	{
		if( SaveCurLib( &g_Lib,_PROMPT_SAVE_ )!=0x0 )
			goto EXIT;
	}

	if( m_CreatePath.IsEmpty( ) )	{
		LPMALLOC pMalloc;   
		if(::SHGetMalloc(&pMalloc)!=NOERROR ) 
			goto EXIT;
		GetCurrentDirectory(MAX_PATH,szCurOld);
		bi.hwndOwner   =   this->m_hWnd;   
		bi.pidlRoot   =   NULL;   
		bi.pszDisplayName   =   szDir;   
		bi.lpszTitle   =   strDlgTitle;   
		bi.ulFlags   =   BIF_RETURNONLYFSDIRS;   
		bi.lpfn   =   _SHBFF_call_back_;   
		bi.lParam   =   0;   		bi.iImage   =   0;   

		pidl = SHBrowseForFolder(&bi);  
		if(pidl!= NULL && SHGetPathFromIDList(pidl,szDir) )  { 
			pMalloc->Free(pidl);   
		}else	{
			isCreate=-1;
		}
		pMalloc->Release();  
	}else	{		//	!!testing!!
		wcscpy( szDir,m_CreatePath );
		m_CreatePath = _T("");
	}

	if( isCreate!=0 )
		goto EXIT;
	{		//to use CWaitCursor
	CWaitCursor wait;
	GIP_LIB_clear( &g_Lib );
	G_FVector_clear( &(g_fv) );
	int M_normal=(int)(GIRIST_control[GIRIS_M_NORMAL]),N_normal=(int)(GIRIST_control[GIRIS_N_NORMAL]);
	int nShift=(int)(GIRIST_control[GIRIS_FV_SIFT]);
	G_FVector_init( &(g_fv),M_normal,N_normal,nShift,0x0 );
	int nRet = GIP_LIB_init( &g_Lib,G_STR2TCHAR(szDir),&(g_fv),GIRIST_control,GIP_PATH_DIRECTORY );
	if(nRet!=0x0)	{
		ASSERT( FALSE );
	}
	G_DCRIMI_set( &g_Lib,&g_dcrimi,g_Lib.D_span,g_Lib.D_inter,g_Lib.D_intra,g_Lib.CRR,0x0 );

/*
	int no=1;
	while( 1 )	{		//避免与已有的重复
	//	if( GetFileAttributes(sLibPath)==INVALID_FILE_ATTRIBUTES )           {     
	//		break;                  
    //   }
	};*/
	sLibPath.Format( _T("%sUntitled_%d.ilb"),g_sLibPath,1 );
	_tcscpy( g_Lib.sPath,G_STR2TCHAR(sLibPath) );
	BIT_SET( g_Lib.status,G_LIB_UNTITLE );
	g_Lib.isUpdate = g_Lib.nEntry>0;
	DumpToFileInit( G_STR2TCHAR(sLibPath),0x0 );
//	int ret = GIP_LIB_save( &g_Lib,G_STR2TCHAR(sLibPath), 0x0 );
//	ASSERT( ret==0 );
	if( g_hTarget!=GIP_NULL )	{
		GIP_free( g_hTarget );		g_hTarget=GIP_NULL;
	}	
	Update( _UPDATE_TREE_ );
/*
	m_tabTree.m_tree.OnUpdate( G_ENTRY_ALL );
	m_tabWarn.m_tree.OnUpdate( G_ENTRY_WARNING );
	m_tabFail.m_tree.OnUpdate( G_ENTRY_FAILED );
#ifdef _GIRIST_DEVELOP_
	m_tabSort.m_tree.OnUpdate( G_ENTRY_ALL );
#endif
*/	}

	//	SetWindowText( sPath );
EXIT:
	;
}

/*
	v0.1	cys
		4/13/2009
*/
void CGiristDlg::OnLibClose( ) {
	// TODO: Add your command handler code here
	TCHAR   szDir[MAX_PATH],szCurOld[MAX_PATH];   
	BROWSEINFO   bi;   
	ITEMIDLIST   *pidl; 
	int isCreate=0;
	CString strDlgTitle=_T("select the path of iris files:"),sLibPath,str;   

	CWaitCursor wait;
	if( g_Lib.isUpdate > 0 )	{
		if( SaveCurLib( &g_Lib,_PROMPT_SAVE_ )!=0x0 )
			goto EXIT;
	}
	GIP_LIB_clear( &g_Lib );
	G_FVector_clear( &(g_fv) );
	G_DCRIMI_set( &g_Lib,&g_dcrimi,g_Lib.D_span,g_Lib.D_inter,g_Lib.D_intra,g_Lib.CRR,0x0 );


	sLibPath.Format( _T("%sUntitled_%d.ilb"),g_sLibPath,1 );
	_tcscpy( g_Lib.sPath,G_STR2TCHAR(sLibPath) );
	BIT_SET( g_Lib.status,G_LIB_UNTITLE );
	g_Lib.isUpdate = g_Lib.nEntry>0;

	g_hSource=GIP_NULL;
	if( g_hTarget!=GIP_NULL )	{
		GIP_free( g_hTarget );		g_hTarget=GIP_NULL;
	}	

	Update( _UPDATE_TREE_ );

	//	SetWindowText( sPath );
EXIT:
	;

}

/*
	v0.1	cys
		3/2/2009
*/
void CGiristDlg::OnLibBuild( ) {
	G_DUMP_CODE oldDump = g_DumpCode;

	g_DumpCode = _DUMP_USER_CTRL;
	CDlgBuild dlg( &g_Lib );
	if( dlg.DoModal( )!=IDOK )	{
		goto END;
	}
	m_tabWarn.m_tree.OnUpdate( G_ENTRY_WARNING );
	m_tabFail.m_tree.OnUpdate( G_ENTRY_FAILED );
#ifdef _GIRIST_DEVELOP_
	m_tabSort.m_tree.OnUpdate( G_ENTRY_ALL );
#endif
END:
	g_DumpCode = oldDump;
}

/*
	v0.1	cys
		3/2/2009
*/
void CGiristDlg::Update( int flag ){
	GIP_LIB *hLib = &(g_Lib);
	int mode=(int)(hLib->user_control[GIRIS_MODE_]);

	if( BIT_TEST(flag,_UPDATE_TREE_) )	{
		m_tabTree.m_tree.OnUpdate( G_ENTRY_ALL );
		m_tabWarn.m_tree.OnUpdate( G_ENTRY_WARNING );
		m_tabFail.m_tree.OnUpdate( G_ENTRY_FAILED );
	#ifdef _GIRIST_DEVELOP_
		m_tabSort.m_tree.OnUpdate( G_ENTRY_ALL );
	#endif
	}

	switch( mode )	{
	case GIRIS_MODE_AUTHEN:
		break;
	case GIRIS_MODE_IDENTIFY:
		break;
	default:
//		m_tiS.m_hEntry = g_hCurEntry;
//		m_tiT.m_hEntry = g_hTarget;
		break;
	}
	m_tiS.Update( 0x0 );
	m_tiS.SetWindowText( _T("girist_source") );
	m_tiT.Update( 0x0 );
}

/*
	v0.1	cys
		3/9/2009
*/
void CGiristDlg::OnLibOpen( ) {
//	SaveCurLib( &g_Lib,_PROMPT_SAVE_ );
	if( g_Lib.isUpdate!=0 && SaveCurLib( &g_Lib,_PROMPT_SAVE_ )!=0x0 )	{
		;
	}else	{	

		CFileDialog dlg( TRUE,NULL,NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
				_T("Iris Lib Files(*.ilb)|*.ilb|")
				_T("All Files(*.*)|*.*|")
			);
		if( dlg.DoModal( ) == IDOK )	{		
			CString sPath=dlg.GetPathName( );
			CWaitCursor wait;
			LoadIrisLib( &g_Lib,sPath,0x0 );

			g_hSource=GIP_NULL;
		//	m_tiS.m_hEntry = g_hCurEntry;
		//	m_tiS.Update( 0x0 );
			if( g_hTarget!=GIP_NULL )	{
				GIP_free( g_hTarget );		g_hTarget=GIP_NULL;
			}	
			Update( _UPDATE_TREE_ );

		}
	}

}

/*
	注意:
		1 hLib的clear由调用程序负责

	v0.1	cys
		3/7/2009
*/
void CGiristDlg::LoadIrisLib( GIP_LIB *hLib,CString sLibPath,int flag ) {
	CString str,str_1;
	G_FVector_ *hfv=&g_fv;
//	GIP_LIB_clear( &g_Lib );
	int ret=GIP_LIB_load( hLib,hfv,G_STR2TCHAR( sLibPath ),0x0 );
	if( ret!=GIP_OK )	{
		str_1.LoadStringW( IDS_LIB_LOAD_ERR );
		str.Format( _T("%s \"%s\".\tERR CODE：%d"),str_1,sLibPath,ret );
		MessageBox( str,g_sEXE,MB_OK );	
	}else	{
		str = GetFileTitle( hLib->sPath );
		DumpToFileInit( G_STR2TCHAR(str),0x0 );

//		hfv->nShift = (int)(hLib->user_control[GIRIS_FV_SIFT]);
		G_DCRIMI_set( hLib,&g_dcrimi,hLib->D_span,hLib->D_inter,hLib->D_intra,hLib->CRR,0x0 );

		m_tabTree.m_tree.OnUpdate( G_ENTRY_ALL );
		m_tabWarn.m_tree.OnUpdate( G_ENTRY_WARNING );
		m_tabFail.m_tree.OnUpdate( G_ENTRY_FAILED );
#ifdef _GIRIST_DEVELOP_
		m_tabSort.m_tree.OnUpdate( G_ENTRY_ALL );
#endif
		m_tabProj.ActivateSSLPage( 0 );
	}
}

/*
	v0.1	cys
		3/9/2009
*/
void CGiristDlg::OnLibSave( ) {
	SaveCurLib( &g_Lib,0x0 );	

}

/*
	仅当isUpdate==1时,保存
	v0.1	cys
		3/8/2009
*/
int CGiristDlg::SaveCurLib( GIP_LIB *hLib,int flag )	{
#ifdef _GRUS_RENSTAR_EXPERT_
#endif
	int ret=0,code;
	CString str,sPath;

	if( hLib->isUpdate==0 )	
	{	ret=-1;	goto END;		}

	if( BIT_TEST( flag,_PROMPT_SAVE_ ) )	{
		if( BIT_TEST( hLib->status,G_LIB_UNTITLE ) )	{
			str.LoadString( IDS_SVAE_UNTITLED );
		}	else	{
			str.LoadString( IDS_SVAE_CHANGE );
			str = str+hLib->sPath;
		}
		code = MessageBox( str,g_sEXE,MB_YESNOCANCEL );
	}else
		code = IDYES;

	switch( code )	{
	case IDCANCEL:
		ret=-4;	
		break;
	case IDYES:
		if( BIT_TEST( hLib->status,G_LIB_UNTITLE ) )	{
			CFileDialog dlg( FALSE,NULL,NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
					_T("Iris Lib Files(*.ilb)|*.ilb|")
					_T("All Files(*.*)|*.*|")
				);
			if( dlg.DoModal( ) != IDOK )	
			{	ret=-2;	goto END;		}	
			sPath=dlg.GetPathName( );
		}else
			sPath = hLib->sPath;
		code = GIP_LIB_save( hLib,G_STR2TCHAR(sPath),0x0 );
		if( code==0 )	{
			BIT_RESET( hLib->status,G_LIB_UNTITLE );
			hLib->isUpdate = 0;
			_tcscpy( hLib->sPath,sPath );
			str.LoadString( IDS_SVAE_SUCCEED );
			SetSBInfo( ID_SEPARATOR,str+sPath );
		}else	{		
			ASSERT( FALSE );		//需提示错误原因
			ret=-3;		goto END;		
		}
		break;
	case IDNO:
		break;
	default:
		break;
	}

	m_tabTree.m_tree.OnUpdate( G_ENTRY_ALL );
	m_tabWarn.m_tree.OnUpdate( G_ENTRY_WARNING );
	m_tabFail.m_tree.OnUpdate( G_ENTRY_FAILED );
#ifdef _GIRIST_DEVELOP_
	m_tabSort.m_tree.OnUpdate( G_ENTRY_ALL );
#endif


END:
	return ret;
}

/*
	v0.1	cys
		3/9/2009
*/
void CGiristDlg::OnLibDistri( ) {
	CDlgDistri dlg( &g_Lib );
	if( dlg.DoModal( )!=IDOK )	{
		return;
	}
}

/*
	输出一些信息给用户，重点是hamming distance distribution

	v0.1	cys
		4/13/2009
*/
void CGiristDlg::OnLibInfo( ) {
/*	CString sPath,sTitle,sLine;

	sTitle = GetFileTitle( g_Lib.sPath );
	sPath.Format( _T(".\\%s_info.txt"),sTitle );
	CStdioFile file( sPath,CStdioFile::modeCreate|CStdioFile::modeWrite );
	GIP_LIB *hLib = &g_Lib;
	int nEntry=hLib->nEntry,nSample=hLib->nSample,i,nRet,cls_old=-1;	
	int M_normal,N_normal,mod;
	GIP_LIB_entry *hEntry;
	double s=1.0/CLOCKS_PER_SEC,*w_sort,s_N=1.0/48/512,Mp_0=100,Ms_0=100,*control=hLib->user_control;

	sLine.Format( _T("\r\n\"%s\"\n"),hLib->sPath );
	file.WriteString( sLine );
	sLine.Format( _T("nEntry=%d,nClass=%d\n"),nEntry,hLib->nClass );
	file.WriteString( sLine );
	sLine.Format( _T("nOK=%d,nFail=%d\n"),nSample,nEntry-nSample );
	file.WriteString( sLine );

#ifdef _GIRIST_DEVELOP_
	mod = (int)(control[GIRIS_MODE_]);
	M_normal=(int)(control[GIRIS_M_NORMAL]),		N_normal=(int)(control[GIRIS_N_NORMAL]);
	sLine.Format( _T("mode=%d,M_normal=%d,N_normal=%d\n"),mod,M_normal,N_normal );	
	file.WriteString( sLine );	

	sLine.Format( _T("\r\n no  p_r p_c R_p s_r s_c R_s \n") );
	file.WriteString( sLine );	
	for( i = 0; i < nEntry; i++ )	{
		hEntry = hLib->arrEntry+i;
		if( hEntry->cls!=cls_old )	{
			sLine.Format( _T("%s\n"),hEntry->sFilePath );
			file.WriteString( sLine );	
			cls_old = hEntry->cls;
		}
		sLine.Format( _T("%-4d %3g %3g %3g %3g %3g %3g \n"),i+1,
			hEntry->p_r,hEntry->p_c,hEntry->p_radius,hEntry->s_r,hEntry->s_c,hEntry->s_radius
		 );	
		file.WriteString( sLine );	
	}
#endif
	file.Close( );
*/
	ShellExecute(NULL,NULL,g_sDumpPath,NULL,NULL,SW_SHOWNORMAL); 
}


/*
	one 2 one compare

	v0.1	cys
		4/13/2009
	v0.2	cys
		4/26/2009
*/
void CGiristDlg::OnIrisAuthen( )	{
	GIP_LIB *hLib=&g_Lib;
	CString sPath_1,sPath_2,str,str_1;
	int mode=(int)(g_Lib.user_control[GIRIS_MODE_]),nResult,i;

	if( mode==GIRIS_MODE_AUTHEN && g_Lib.nEntry==2 && isTargetCheck==TRUE )	{
		isTargetCheck=FALSE;
		sPath_1 = hLib->arrEntry[0].sFilePath;
		sPath_2 = hLib->arrEntry[1].sFilePath;
	}else	{
		OnLibClose( );
		do	{
			CFileDialog dlg( TRUE,NULL,NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
				_T("Iris Image Files(*.bmp;*.jpg)|*.bmp;*.jpg|")
				_T("All Files(*.*)|*.*|")
			);
			dlg.m_ofn.lpstrTitle=_T("Please select first iris");   

			if( dlg.DoModal( ) != IDOK )	
				return;
			sPath_1=dlg.GetPathName( );
			GIP_LIB_init( hLib,G_STR2TCHAR(sPath_1),&(g_fv),GIRIST_control,0x0 );
			ASSERT( hLib->nEntry> 0 );
		}while( hLib->nEntry==0 );

		do	{
			CFileDialog dlg( TRUE,NULL,NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
				_T("Iris Image Files(*.bmp;*.jpg)|*.bmp;*.jpg|")
				_T("All Files(*.*)|*.*|")
			);
			dlg.m_ofn.lpstrTitle=_T("Please select second iris");   
			if( dlg.DoModal( ) != IDOK )	
				return;
			sPath_2=dlg.GetPathName( );
			GIP_LIB_Append( hLib,G_STR2TCHAR(sPath_2),0,hLib->nClass,0x0 );
			hLib->nClass++;
			ASSERT( hLib->nEntry> 1 );
		}while( hLib->nEntry==1 );

	//	BIT_SET( hLib->status,G_LIB_UNTITLE );
	//	hLib->isUpdate = 1;
		str=GetFileTitle(sPath_1)+_T("_")+GetFileTitle(sPath_2)+_T("_comp.ilb");
		_tcscpy( hLib->sPath,G_STR2TCHAR(str) );
		hLib->user_control[GIRIS_MODE_]=GIRIS_MODE_AUTHEN;
		str = GetFileTitle( hLib->sPath );
		DumpToFileInit( G_STR2TCHAR(str),0x0 );
		G_FVector_clear( &(g_fv) );
		int M_normal=(int)(GIRIST_control[GIRIS_M_NORMAL]),N_normal=(int)(GIRIST_control[GIRIS_N_NORMAL]);
		int nShift=(int)(GIRIST_control[GIRIS_FV_SIFT]);
		G_FVector_init( &(g_fv),M_normal,N_normal,nShift,0x0 );
	}
	{
	CWaitCursor wait;
	nResult = G_IRIS_location( hLib,0x0 );
	for( i = 0 ; i < 2; i++ )	{
		if( hLib->arrEntry[i].status!=G_ENTRY_FAILED )
			hLib->arrEntry[i].status = G_ENTRY_OK;
		else	{
			str_1.LoadString( IDS_LOCATION_FAIL );
			str.Format( str_1,hLib->arrEntry[i].sFilePath );
			goto G_EXIT;
		}
	}
	nResult = G_IRIS_coding( hLib,0x0 );
	double dist_4[8],*inter_0=dist_4,dist;
	nResult = GIP_LIB_discriminate( hLib,&g_dcrimi,dist_4,0x0 );

//	m_tiS.m_hEntry = &(hLib->arrEntry[0]);
//	m_tiT.m_hEntry = &(hLib->arrEntry[1]);
	g_hSource = &(hLib->arrEntry[0]);
//	g_hTarget = &(hLib->arrEntry[1]);
	if( g_hTarget!=GIP_NULL )	{
			GIP_free( g_hTarget );		g_hTarget=GIP_NULL;
		}
	g_hTarget=(GIP_LIB_entry*)GIP_alloc( sizeof(GIP_LIB_entry) );
	GIP_MEMCOPY( g_hTarget,hLib->arrEntry+1,sizeof(GIP_LIB_entry));
	Update( _UPDATE_TREE_ );
	m_tabProj.ActivateSSLPage( 0 );

	g_dcrimi.sep = hLib->user_control[GIRIS_T_SEP];
	dist = (inter_0[0]+inter_0[1])/2;
	if( dist<g_dcrimi.sep )	{
		str_1.LoadString( IDS_IRIS_AUTHEN_OK );
		str.Format( str_1,
			GetFileTitle(sPath_1),GetFileTitle(sPath_2),dist );
	}else	{
		str_1.LoadString( IDS_IRIS_AUTHEN_FAIL );
		str.Format( str_1,
			GetFileTitle(sPath_1),GetFileTitle(sPath_2),dist );
	}
	str_1.LoadString( IDS_IRIS_AUTHEN_THRESH );
	str += _T("\n")+str_1;
	}
G_EXIT:
	MessageBox( str,g_sEXE,MB_OK );	
}

/*
	one 2 many(LIB) compare

	v0.1	cys
		4/13/2009
*/
void CGiristDlg::OnIrisIdentify( )	{
	GIP_LIB *hLib=&g_Lib;
	CString str,str_1,sPath,sPath_2;
	G_IRIS_ *hIris=GIP_NULL;

	if( hLib->nSample==0 )	{
		str_1.LoadString( IDS_LIB_NOOK );
		str.Format( str_1 );
		MessageBox( str,g_sEXE,MB_OK );	
		return;
	}
	int ret,no=-1;
	double *dist=(double*)GIP_alloc( sizeof(double)*hLib->nEntry ),a;
	if( g_hTarget!=GIP_NULL && isTargetCheck==TRUE )	{
		isTargetCheck=FALSE;
//		GIP_free( g_hTarget );		g_hTarget=GIP_NULL;
	}else	{
		CFileDialog dlg( TRUE,NULL,NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
			_T("Iris Image Files(*.bmp;*.jpg)|*.bmp;*.jpg|")
			_T("All Files(*.*)|*.*|")
		);
		if( dlg.DoModal( ) != IDOK )	
			return;

		if( g_hTarget!=GIP_NULL )	{
			GIP_free( g_hTarget );		g_hTarget=GIP_NULL;
		}
		g_hTarget=(GIP_LIB_entry*)GIP_alloc( sizeof(GIP_LIB_entry) );
		GIP_MEMCLEAR( g_hTarget,sizeof(GIP_LIB_entry) );
		sPath=dlg.GetPathName( );
		_tcscpy( g_hTarget->sFilePath,G_STR2TCHAR(sPath) );
	}
	{	
	CWaitCursor wait;
	hIris = G_IRIS_on_entry( hLib->user_control,g_hTarget,0x0,0x0,&ret );
	if( ret==GIP_OK )	{
		no=G_IRIS_identify( hLib,hIris,dist,0x0 );
	}else	{
	//	str_1.LoadString( IDS_LIB_NOOK );
	//	str.Format( str_1 );
	//	MessageBox( str,g_sEXE,MB_OK );	
	}
	double sep=g_dcrimi.sep==0.0 ? hLib->user_control[GIRIS_T_SEP] : g_dcrimi.sep;
	if( no != -1 )	{
		sPath_2 = hLib->arrEntry[no].sFilePath;
		a = dist[no];
		g_hTarget->x = a;
		if( a<sep )	{
			str.LoadString( IDS_IRIS_IDENTIFY_OK );
			str_1.Format( str,sPath_2,GetFileTitle(sPath),GetFileTitle(sPath_2),a,sep );
		}else	{
			str.LoadString( IDS_IRIS_IDENTIFY_FAIL );
			str_1.Format( str,sPath_2,GetFileTitle(sPath),GetFileTitle(sPath_2),a,sep );
		}	
		g_hSource = &(hLib->arrEntry[no]);
	}else	{
		g_hSource = GIP_NULL;
	}
	str = str_1;
	}
	MessageBox( str,g_sEXE,MB_OK );	

	if( hIris!=GIP_NULL )	{
		G_IRIS_clear( hIris );			GIP_free( hIris );
	}
	GIP_free( dist );

	Update( _UPDATE_TREE_ );

}

/*
	v0.1	cys
		4/15/2009
*/
void CGiristDlg::OnNMClickTabTarget(NMHDR *pNMHDR, LRESULT *pResult)
{
	// TODO: Add your control notification handler code here
#ifdef _GIRIST_DEVELOP_
	isTargetCheck = TRUE;
	int mode=(int)(g_Lib.user_control[GIRIS_MODE_]);

	if( g_hTarget!=GIP_NULL )	{
		if( g_Lib.nEntry>0 && mode!=GIRIS_MODE_AUTHEN )	{
//			BIT_SET( g_hTarget->status,G_ENTRY_SELECT );
			OnIrisIdentify( );
//			BIT_RESET( g_hTarget->status,G_ENTRY_SELECT );
		}
	}
	if( g_Lib.nEntry==2 && mode==GIRIS_MODE_AUTHEN )	{
		OnIrisAuthen( );
	}
#endif
	*pResult = 0;
}

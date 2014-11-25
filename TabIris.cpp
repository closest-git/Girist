// TabIris.cpp : implementation file
//

#include "stdafx.h"
#include "Girist.h"
#include "TabIris.h"
#include "./pch/gip_def.h"
#include "./girist/giris_core.h"
#include "./util/gip_util.h"


// CTabIris dialog

IMPLEMENT_DYNAMIC(CTabIris, CDialog)

CTabIris::CTabIris(CWnd* pParent /*=NULL*/)
	: CTabPageSSL(CTabIris::IDD, pParent)
{
//	m_hEntry=GIP_NULL;
	m_hLib = &g_Lib;
}

CTabIris::~CTabIris()
{
}

void CTabIris::DoDataExchange(CDataExchange* pDX)
{
	CTabPageSSL::DoDataExchange(pDX);

	DDX_Control(pDX, IDC_CUSTOM_NORMAL, m_ctrlNormal);
	DDX_Control(pDX, IDC_CUSTOM_EYE, m_ctrlEye);
	DDX_Control(pDX, IDC_ST_INFO, m_stInfo);
	DDX_Control(pDX, IDC_EDIT1, m_edNotes);
	DDX_Control(pDX, IDC_BUTTON1, m_BT_1);
	DDX_Control(pDX, IDC_BUTTON2, m_BT_2);
}


BEGIN_MESSAGE_MAP(CTabIris, CTabPageSSL)
	ON_STN_CLICKED(IDC_STATIC_NORMAL, &CTabIris::OnStnClickedStaticNormal)
	ON_BN_CLICKED(IDC_BUTTON1, &CTabIris::OnBnClickedButton_1)
	ON_BN_CLICKED(IDC_BUTTON2, &CTabIris::OnBnClickedButton_2)
END_MESSAGE_MAP()


// CTabIris message handlers

void CTabIris::OnStnClickedStaticNormal()
{
	// TODO: Add your control notification handler code here
}

/*
	CGiristDlg message handlers
	
	v0.1	cys
		2/28/2009
*/
BOOL CTabIris::OnInitDialog()
{
	CDialog::OnInitDialog();
	int mode=m_hLib==GIP_NULL ? GIRIS_MODE_DEFAULT : (int)(m_hLib->user_control[GIRIS_MODE_]);
//	m_board.Create(this, CRect(50, 50, 150, 150), 0 );
	m_ctrlEye.sNote = _T("ctrl_Eye");
	m_ctrlNormal.sNote = _T("ctrl_Normal");

	m_ctrlEye.ReCalculateDimensions( );
	m_ctrlNormal.ReCalculateDimensions( );

	CString sInfo=_T("");
	m_stInfo.SetWindowTextW( sInfo );
	if( type==TAB_IRIS_SOURCE )	{
		m_BT_1.EnableWindow( FALSE );
		m_BT_2.EnableWindow( FALSE );
		m_BT_1.ShowWindow( SW_HIDE );
		m_BT_2.ShowWindow( SW_HIDE );
	}else	{
		m_BT_1.SetWindowText( _T("Open File...") );
		m_BT_2.SetWindowText( _T("Target:Source") );
		m_BT_1.EnableWindow( mode!=GIRIS_MODE_AUTHEN && m_hLib!=GIP_NULL && m_hLib->nEntry>0  );
		m_BT_2.EnableWindow( g_hTarget!=GIP_NULL && g_hSource!=GIP_NULL  );
	}

	return TRUE;  // return TRUE  unless you set the focus to a control
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/10/2009
*/
void CTabIris::Clear( int flag ){
	m_stInfo.SetWindowTextW( _T("") );
	
	m_ctrlEye.p_r=-1.0;		m_ctrlEye.p_c=-1.0;		m_ctrlEye.p_radius=-1.0;		
	m_ctrlEye.s_r=-1.0;		m_ctrlEye.s_c=-1.0;		m_ctrlEye.s_radius=-1.0;
	m_ctrlEye.LoadObj( GIP_NULL,0x0 );
	m_ctrlNormal.LoadObj( GIP_NULL,0x0 );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/3/2009
*/
void CTabIris::Update( int flag ){
	G_IRIS_ *hIris=GIP_NULL;
	int ret=-1;
	GIP_LIB_entry *hEntry=GIP_NULL;
	hEntry= type==TAB_IRIS_SOURCE ? g_hSource : type==TAB_IRIS_TARGET ? g_hTarget : GIP_NULL;

	CString sInfo=_T(""),str;
	if( hEntry != GIP_NULL )	{
		CWaitCursor wait;
		hIris = G_IRIS_on_entry( m_hLib->user_control,hEntry,0x0,flag,&ret );
/*		if( m_hEntry->status==G_ENTRY_OK )	{
		}else	if( 1 )	{
//			hIris = G_IRIS_process( m_hLib->user_control,m_hEntry->sFilePath,0x0,flag,&ret );
			m_hEntry->p_r=hIris->p_r;	m_hEntry->p_c=hIris->p_c;	m_hEntry->p_radius=hIris->p_radius;		
			m_hEntry->s_r=hIris->s_r;	m_hEntry->s_c=hIris->s_c;	m_hEntry->s_radius=hIris->s_radius;		
			m_hEntry->p_match=hIris->p_match;						m_hEntry->s_match=hIris->s_match;	
		}*/
		if( ret==GIP_OK )	{
		}else	{
		}
		
	}else
		Clear( 0x0 );		

	if( hIris!=GIP_NULL )	{
		sInfo.Format( 
			_T("\"%s\"\r\n\r\n")
			_T("inner boundary:(%g,%5.3f)\r\nouter boundary:(%g,%5.3f)\r\n"),
			GetFileTitle(hEntry->sFilePath),
			/*hIris->p_r,hIris->p_c,*/hIris->p_radius,hIris->p_match,
			/*hIris->s_r,hIris->s_c,*/hIris->s_radius,hIris->s_match			
		);
//		if( m_hEntry==g_hTarget )	{
//			str.Format( _T("\r\nThe distance to the source iris is %4.3g\r\n"),m_hEntry->x );
//		}
		sInfo+=str;

		m_ctrlEye.p_r=hIris->p_r;	m_ctrlEye.p_c=hIris->p_c;	m_ctrlEye.p_radius=hIris->p_radius;		
		m_ctrlEye.s_r=hIris->s_r;	m_ctrlEye.s_c=hIris->s_c;	m_ctrlEye.s_radius=hIris->s_radius;
		m_ctrlEye.LoadObj( hIris->hObj,0x0 );
		m_ctrlNormal.LoadObj( &(hIris->normal),0x0 );

		G_IRIS_clear( hIris );			GIP_free( hIris );
	}
	m_stInfo.SetWindowTextW( sInfo );

	if( type==TAB_IRIS_TARGET )	{
		int mode=m_hLib==GIP_NULL ? GIRIS_MODE_DEFAULT : (int)(m_hLib->user_control[GIRIS_MODE_]);
		m_BT_1.EnableWindow( mode!=GIRIS_MODE_AUTHEN && m_hLib!=GIP_NULL && m_hLib->nEntry>0 );
		m_BT_2.EnableWindow( g_hTarget!=GIP_NULL && g_hSource!=GIP_NULL  );
	}
}

/*
	BT_1:
		open file

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/26/2009
*/
void CTabIris::OnBnClickedButton_1()
{
	// TODO: Add your control notification handler code here
	if( type==TAB_IRIS_SOURCE )	
		return;
	
	CFileDialog dlg( TRUE,NULL,NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
		_T("Iris Image Files(*.bmp;*.jpg)|*.bmp;*.jpg|")
		_T("All Files(*.*)|*.*|")
	);
	if( dlg.DoModal( ) != IDOK )	
		return;

	CString sPath;
	if( g_hTarget!=GIP_NULL )	{
		GIP_free( g_hTarget );		g_hTarget=GIP_NULL;
	}
	g_hTarget=(GIP_LIB_entry*)GIP_alloc( sizeof(GIP_LIB_entry) );
	GIP_MEMCLEAR( g_hTarget,sizeof(GIP_LIB_entry) );
	sPath=dlg.GetPathName( );
	_tcscpy( g_hTarget->sFilePath,G_STR2TCHAR(sPath) );	
//	m_hEntry = g_hTarget;
	Update( 0x0 );
	
	m_BT_2.EnableWindow( g_hTarget!=GIP_NULL && g_hSource!=GIP_NULL  );
}

/*
	BT_1:
		Distance 2 source and target

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/26/2009
*/
void CTabIris::OnBnClickedButton_2()
{
	// TODO: Add your control notification handler code here
	if( type==TAB_IRIS_SOURCE )	
		return;

	double dist,sep;
	int ret;
	CString str,str_1;
	{
	CWaitCursor wait;
	ret = G_IRIS_distance_2( m_hLib->user_control,m_hLib->hfv,g_hSource,g_hTarget,&dist,0x0 );
	}
	if( ret == GIP_OK )	{
		sep=g_dcrimi.sep==0.0 ? m_hLib->user_control[GIRIS_T_SEP] : g_dcrimi.sep;
		str_1.LoadString( IDS_IRIS_DISTANCE_2 );
		str.Format( str_1,GetFileTitle(g_hSource->sFilePath),GetFileTitle(g_hTarget->sFilePath),dist,sep );
		if( dist<=sep )
			str = _T("SAME EYE\n\n") + str;
		else
			str = _T("DIFFERENT EYE\n\n") + str;
	}else	{
		str_1.LoadString( IDS_IRIS_DISTANCE_2_FAIL );
		str.Format( str_1,ret );
	}
	MessageBox( str,g_sEXE,MB_OK );	
}

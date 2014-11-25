// DlgDistri.cpp : implementation file
//

#include "stdafx.h"
#include <process.h>    
#include "Girist.h"
#include "DlgDistri.h"
#include "util/gip_library.h"
#include "util/gip_util.h"
#include "util/gip_fv.h"
#include "girist/giris_core.h"

/*
	1 If the Hamming distance between two templates is less than the separation point, 
	the templates were generated from the same iris and a match is found. 
	Otherwise if the Hamming distance is greater than the separation point the two templates are 
	considered to have been generated from different irises.

	2 Clearly the separation point will influence the false accept and false reject rates, 
	since a lower separation Hamming distance will decrease FAR while increasing FRR, 
	and vice versa. Therefore, when choosing a separation point it is important to consider 
	both the false accept rate and false reject rate.
*/
#define DISTRI_STAGE_NUM 2
static int _stage_=0;

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
static unsigned __stdcall  _DistriThreadProc(LPVOID lpThreadParameter)	{
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

	G_FV_info( hLib->hfv,0x0 );	
	_stage_ = 1;
	start=clock( );	
	G_PRINTF( _T("1>------ FEATURE EXTRACTION & CODING ------\r\n") );
	nResult = G_IRIS_coding( hLib,0x0 );
	if( nResult != 0x0 || hUserCtrl->bTerminate == TRUE )
		goto TerminateThreadProc;

	_stage_ = 2;
	G_PRINTF( _T("\r\n------ INTRA-CLASS & INTER-CLASS CODING ------\r\n") );
//	nResult = GIP_LIB_scan( hLib,0x0 );
//	if( nResult != 0x0 || hUserCtrl->bTerminate == TRUE )
//		goto TerminateThreadProc;	
	start=clock( );	
	nResult = GIP_LIB_discriminate( hLib,&g_dcrimi,GIP_NULL,0x0 );
	t = (clock()-start)*1.0/CLOCKS_PER_SEC;
	g_dcrimi.time = t;
	G_PRINTF( _T("2>------ TIME %g------\r\n"),t );
	if( nResult != 0x0 || hUserCtrl->bTerminate == TRUE )
		goto TerminateThreadProc;	/**/

TerminateThreadProc:
	if(hUserCtrl->bTerminate == FALSE)
		::PostMessage(hUserCtrl->hThreadWnd,GUM_PROGRESS_THREADCOMPLETED,MAKEWPARAM(nResult,0),0);
	else
		::PostMessage(hUserCtrl->hThreadWnd,GUM_CANCEL_PROGRESSTHREAD,MAKEWPARAM(nResult,0),0);

	hUserCtrl->bAlive = FALSE;
	return 0;    
}

// CDlgDistri dialog

IMPLEMENT_DYNAMIC(CDlgDistri, CDialog)

CDlgDistri::CDlgDistri(GIP_LIB *hLib,CWnd* pParent /*=NULL*/)
	: CDialog(CDlgDistri::IDD, pParent)
	, m_Shift(0)
	, m_WaveLen(0)
{
	m_hLib = hLib;
	GU_PROGRESS_clear( &m_Progress );

}

CDlgDistri::~CDlgDistri()
{
}

void CDlgDistri::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_CUSTOM_DISTRI, m_cDistri);
	DDX_Control(pDX, IDC_CUSTOM_ROC, m_cROC);
	DDX_Control(pDX, IDC_ED_DCRIMI, m_edDrimi);
	DDX_Control(pDX, IDC_PROGRESS1, m_cProgress);
	DDX_Text(pDX, IDC_ED_SHIFT, m_Shift);	
	DDV_MinMaxUInt( pDX, m_Shift,0,15 );	//(0'-15');
	DDX_Text(pDX, IDC_EDIT_SEP, m_Sep);	
	DDV_MinMaxDouble( pDX, m_Sep,0.1,0.99 );	
	DDX_Text(pDX, IDC_ED_WAVELEN, m_WaveLen);	
	DDX_Control(pDX, IDOK, m_btOK);
	DDX_Control(pDX, IDCANCEL, m_btCancel);
	DDX_Control(pDX, IDC_BT_DISTRI, m_btSwitch);
	DDX_Control(pDX, IDC_CB_LEVEL, m_cbLevel);
	DDX_Control(pDX, IDC_ST_SEP, m_stSep);
}

/*
	bAlive=TRUE:	工作线程在运行，只可以按m_btSwitch

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
void CDlgDistri::CtrlOnState( BOOL bAlive )	{
	CString	sBt = bAlive==TRUE ? _T("STOP") : _T("Re Coding");
	m_btSwitch.SetWindowTextW( sBt );

//	m_btCancel.EnableWindow( !bAlive );
	m_btOK.EnableWindow( !bAlive );
	m_cbLevel.EnableWindow( !bAlive );
	m_edDrimi.EnableWindow( !bAlive );
	m_cDistri.EnableWindow( !bAlive );
	m_cROC.EnableWindow( !bAlive );
}

/*
	注意：
		1 存在无法分配线程的可能

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/14/2009
*/
int CDlgDistri::StartThread( )	{
	unsigned threadID;
	CString str;

	GU_PROGRESS_clear( &m_Progress );
	m_Progress.hThreadWnd = this->GetSafeHwnd( );
	m_Progress.pUserParam = (LPVOID)(m_hLib);
	m_hThread = (HANDLE)_beginthreadex(NULL, 0x0, _DistriThreadProc,&(m_Progress), CREATE_SUSPENDED, &threadID );
	if( m_hThread==NULL )	{		//存在失败的可能
		str.LoadStringW( IDS_THREAD_CREATE_FAIL );
		str.LoadStringW( IDS_GIRIST_RESTART );
		MessageBox( str,g_sEXE,MB_OK );	
		return -1;
	}
	m_Progress.bAlive = TRUE;
	CtrlOnState( m_Progress.bAlive );

	G_FVector_ *hfv=m_hLib->hfv;
	hfv->nShift = m_Shift;
	m_hLib->user_control[GIRIS_LG_WAVELEN] = m_WaveLen;

	ResumeThread( m_hThread );

	return 0;
}

/*
	Copyright 2008-present, Grusoft.	
	v0.1	cys
		3/12/2009
*/
CString _Drimi_Info( G_DCRIMI_ *hDcrimi,G_FVector_ *hfv,int flag )	{
	CString info;
	GIP_LIB *hLib=(GIP_LIB *)(hDcrimi->hBase);
	int nEntry=hLib->nEntry;
	double *frr=hDcrimi->f_rr_8+3,*hd=hDcrimi->hd_8+3;		//分别对应1.0e-4,1.0e-3,1.0e-2

	if( hLib->nClass==1 )	{
		info.Format( _T("\t%8s %8s\r\n%8s %8g \t----\r\n%8s %8.3g\t----\r\n%8s %8.3g\t----\r\n")
				_T("-----------------------------------\r\n\r\n"),
//				_T("decidability=%g\r\nFAR=%g\r\nFRR=%g\r\n"),
			_T("intra"),_T("inter"),
			_T("COMP"),hDcrimi->nz_a, 
			_T("MEAN"),hDcrimi->mean_a,
			_T("DEVIA"),hDcrimi->devia_a
		);	
	}else	{
		info.Format( _T("\t\t\t\t|\t\t\t|\t%8s %8s\r\n")
				_T("FAR\t1.0E-4\t1.0E-3\t1.0E-2\t|")_T("CRR=%4.3g%%\t\t|")_T("%8s %8d %8d\r\n")
				_T("FRR\t%-5.3g%%\t%-5.3g%%\t%-5.3g%%\t|")_T("decidability=%5.3g\t\t|")_T("%8s %8.3g %8.3g\r\n") 
				_T("HD\t%-4.3g\t%-4.3g\t%-4.3g\t|")_T("EER=%-6.4g\t\t|")_T("%8s %8.3g %8.3g\r\n") 
				_T("---------------------------------------------------------------------------------------------------------\r\n")
				_T("Location time=%5.3g(s),Matching Speed=%d(irises/sec)\r\n")
#ifdef _GIRIST_DEVELOP_
				_T("%d: f_M=%d(%d),f_N=%d,V_len=%d\r\n")
				_T("%d\r\n")
#endif
				,
		_T("intra"),_T("inter"),
		hDcrimi->CRR*100,_T("COMP"),(int)(hDcrimi->nz_a),(int)(hDcrimi->nz_r), 
		frr[0]*100,frr[1]*100,frr[2]*100,hDcrimi->D_,_T("MEAN"),hDcrimi->mean_a,hDcrimi->mean_r,
		hd[0],hd[1],hd[2],hDcrimi->rEER,_T("DEVIA"),hDcrimi->devia_a,hDcrimi->devia_r,
		hLib->L_spd,(int)(hLib->M_spd)
#ifdef _GIRIST_DEVELOP_
		,hfv->type,hfv->sf_M,hfv->M_step,hfv->sf_N,hfv->V_len
		,(int)(nEntry*(1-hDcrimi->CRR))
#endif
	
		);
	}
	return info;
}

/*
	Copyright 2008-present, Grusoft.	
	v0.1	cys
		3/10/2009
*/
BOOL CDlgDistri::OnInitDialog()
{
	G_FVector_ *hfv=m_hLib->hfv;
	m_Shift = (int)(hfv->nShift);
	m_WaveLen = (int)(m_hLib->user_control[GIRIS_LG_WAVELEN]);
	m_Sep=g_dcrimi.sep;

	CDialog::OnInitDialog();
	CString str;
	int i,level=1;
	for( i = 0; i < 4; i++ )	{
		str.Format( _T("%d"),level );
		m_cbLevel.AddString( str );
		if( level==hfv->M_step )
			m_cbLevel.SetCurSel( i );
		level *=2;
	}
	m_cDistri.ReCalculateDimensions( );
	m_cROC.ReCalculateDimensions( );

	m_cProgress.SetRange( 0,m_hLib->nEntry*DISTRI_STAGE_NUM );
	m_cDistri.Update( 0x0 );
	m_cROC.Update( 0x0 );
	str.Format( _T("Decision Threshold(FAR=%6lf,FRR=%4.3g)"),g_dcrimi.rFAR,g_dcrimi.rFRR);
	m_stSep.SetWindowText( str );
	str = _Drimi_Info( &g_dcrimi,hfv,0x0 );
	m_edDrimi.SetWindowTextW( str );

	return TRUE;  // return TRUE  unless you set the focus to a control
}

BEGIN_MESSAGE_MAP(CDlgDistri, CDialog)
	ON_BN_CLICKED(IDC_BT_DISTRI, &CDlgDistri::OnBnClickedBtDistri)
	ON_MESSAGE( GUM_PROGRESS_BARUPDATE, OnUserBarUpdate) 
	ON_MESSAGE( GUM_PROGRESS_TEXTUPDATE, OnUserTextUpdate) 
	ON_MESSAGE( GUM_PROGRESS_THREADCOMPLETED, OnThreadCompleted) 
	ON_EN_CHANGE(IDC_EDIT_SEP, &CDlgDistri::OnEnChangeEditSep)
END_MESSAGE_MAP()


// CDlgDistri message handlers
/*
	Copyright 2008-present, Grusoft.	
	v0.1	cys
		3/14/2009
*/
void CDlgDistri::OnBnClickedBtDistri()	{
	// TODO: Add your control notification handler code here
	UpdateData( );
//	GIRIST_control[GIRIS_FV_SIFT] = m_Shift;
//	GIRIST_control[GIRIS_LG_WAVELEN] = m_WaveLen;
	CString str;
	int delay=500;

	if( m_Progress.bAlive )	{
		m_Progress.bTerminate=TRUE;
		while( m_Progress.bAlive==TRUE )	{
			Sleep(delay);
		}
		str.LoadStringW( IDS_BUILD_CANCEL );
//		m_lbInfo.AddString( str );
		CtrlOnState( m_Progress.bAlive );
		m_cProgress.SetPos( 0 );
	}else	{
		G_FVector_ *hfv=m_hLib->hfv;
		int no=m_cbLevel.GetCurSel( ),level=-1;
		if( no!=-1 )	{
			level = 0x1<<no;
			if( level>0 && level<hfv->sf_M*hfv->M_step/2 )	{
				G_FVector_set_step( hfv,level,1 );
			}
		}
		StartThread( );
	}

}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
LRESULT CDlgDistri::OnUserBarUpdate(UINT wParam, LONG lParam)	{
	int pos = wParam;
	m_cProgress.SetPos( pos+m_hLib->nEntry*(_stage_-1) );
	return 0L; 
} 

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
LRESULT CDlgDistri::OnUserTextUpdate(UINT wParam, LONG lParam)	{
//	int cur = m_lbInfo.AddString( (LPCTSTR)lParam );
//	m_lbInfo.SetTopIndex( cur );

	return 0L; 
} 

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
LRESULT CDlgDistri::OnThreadCompleted(UINT wParam, LONG lParam)	{
	int delay=500;

	while( m_Progress.bAlive==TRUE  )	{
		Sleep(delay);
	}
	m_cDistri.Update( 0x0 );
	m_cROC.Update( 0x0 );
	G_FVector_ *hfv=m_hLib->hfv;
	CString str = _Drimi_Info( &g_dcrimi,hfv,0x0 );
	m_edDrimi.SetWindowTextW( str );

	m_hLib->isUpdate = TRUE;
	CtrlOnState( m_Progress.bAlive );

	m_Sep=g_dcrimi.sep;
	UpdateData( FALSE );
	return 0L; 
} 

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/16/2009
*/
void CDlgDistri::OnOK()	{
	UpdateData( );
	GIRIST_control[GIRIS_FV_SIFT] = m_Shift;
	GIRIST_control[GIRIS_LG_WAVELEN] = m_WaveLen;

	CWnd *hFocus = GetFocus( );

	if( GetDlgItem(IDC_ED_SHIFT)==hFocus || GetDlgItem(IDC_ED_WAVELEN)==hFocus  )	{
		CWnd* pWndNext = GetNextDlgTabItem(hFocus);
		if (pWndNext) {
			pWndNext->SetFocus();
		}
	}else	{
	//	m_hLib->isUpdate = TRUE;

		CDialog::OnOK( );
	}
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/12/2009
*/
void CDlgDistri::OnCancel()	{
	CString str;
	int delay=500;

	if( m_Progress.bAlive )	{
		m_Progress.bTerminate=TRUE;
		while( m_Progress.bAlive==TRUE )	{
			Sleep(delay);
		}
		str.LoadStringW( IDS_BUILD_CANCEL );
//		m_lbInfo.AddString( str );
		CtrlOnState( m_Progress.bAlive );
		m_cProgress.SetPos( 0 );
	}

	CDialog::OnCancel( );
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/12/2009
*/
void CDlgDistri::OnEnChangeEditSep(){
	double oldSep=m_Sep;
	CString str;

	UpdateData( );
	if( m_Sep<0.1 || m_Sep>0.99 )	{
		m_Sep = oldSep;
		UpdateData( FALSE );
		return;
	}

	if( oldSep!=m_Sep )	{
		g_dcrimi.sep = m_Sep;
		G_DCRIMI_update_sep(  &g_dcrimi,m_hLib->D_span,m_hLib->D_inter,m_hLib->D_intra,0x0 );

		str.Format( _T("Decision Threshold(FAR=%6lf,FRR=%4.3g)"),g_dcrimi.rFAR,g_dcrimi.rFRR);
		m_stSep.SetWindowText( str );
		m_cDistri.Update( 0x0 );
		m_cROC.Update( 0x0 );
		G_FVector_ *hfv=m_hLib->hfv;
		str = _Drimi_Info( &g_dcrimi,hfv,0x0 );
		m_edDrimi.SetWindowTextW( str );

		m_hLib->isUpdate = TRUE;
	}
}

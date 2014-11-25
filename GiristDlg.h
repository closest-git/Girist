// GiristDlg.h : header file
//

#pragma once

#include ".\BoardControl.h"
#include "afxwin.h"
#include "dlgbars.h"
//#include "abtreectrl.h"
#include "afxcmn.h"
#include "TabCtrlSSL.h"
#include "TabIris.h"
#include "TabTree.h"

enum	{
//SaveCurLib的标志位
	_PROMPT_SAVE_=0x100,
};
enum	{
//Update的标志位
	_UPDATE_TREE_=0x100,
};

// CGiristDlg dialog
class CGiristDlg : public CDialog
{
// Construction
public:
	CGiristDlg(CWnd* pParent = NULL);	// standard constructor
	~CGiristDlg( );
	BOOL Create();

// Dialog Data
	enum { IDD = IDD_GIRIST_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	virtual void OnCancel();
	virtual void OnOK();
	virtual void PostNcDestroy();


// Implementation
protected:
	HICON m_hIcon;
	UINT            m_nIDTracking;
	UINT            m_nIDLastMessage;
	ULONG_PTR		m_gdiplusToken;
	BOOL isTargetCheck;

	CDlgStatusBar   m_statusBar;
	CDlgToolBar     m_toolBar;
	int CreateControls( );
	int SaveCurLib( GIP_LIB *hLib,int flag );

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	afx_msg LRESULT OnSetMessageString(WPARAM wParam, LPARAM lParam = 0L);
	afx_msg void OnIrisAuthen();
	afx_msg void OnIrisIdentify();
	afx_msg void OnLibCreate();
	afx_msg void OnLibSave();
	afx_msg void OnLibBuild( );
	afx_msg void OnLibDistri( );
	afx_msg void OnLibInfo( );
	afx_msg void OnLibOpen();
	afx_msg void OnLibClose();
	afx_msg void OnFileOpen();
	afx_msg void OnToolTest();
	afx_msg void OnToolOption( );
	afx_msg void OnAppAbout();
	afx_msg void OnUpdateTime(CCmdUI* pCmdUI);
	afx_msg void OnSelChangeMode();
	afx_msg void OnSelChangeDepth();
	afx_msg void OnChessRange( UINT );
	afx_msg void OnUpdateLibDistri( CCmdUI* );
	afx_msg void OnUpdateLibBuild( CCmdUI* );
	afx_msg void OnUpdateLibInfo( CCmdUI* );
	afx_msg void OnUpdateLibSave( CCmdUI* );
	afx_msg void OnUpdateIrisIdentify( CCmdUI* );
	afx_msg void OnUpdateIrisAuthen( CCmdUI* );
	afx_msg void OnUpdateLibClose( CCmdUI* );
	DECLARE_MESSAGE_MAP()

public:	
	void ResizeControls();		
	void SetSBInfo( int no,CString info )	{ m_statusBar.SetPaneText( m_statusBar.CommandToIndex(no),info );}	
//	void UpdateABTree()	{		m_tabTree.m_treeAB.OnUpdate( );	}
	void OnCommand( int nArg,... );
	void LoadIrisLib( GIP_LIB *hLib,CString sLibPath,int flag );
	void Update( int flag );

public:
	afx_msg void OnGetMinMaxInfo(MINMAXINFO* lpMMI);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
public:
	CString m_CreatePath;
	CTabCtrlSSL m_tabProj,m_tabSource,m_tabTarget;
	CTabIris m_tiS,m_tiT;
	TAB_TREE m_tabTree,m_tabWarn,m_tabSort,m_tabFail;

	CComboBox	m_cbModes;	

friend class CGiristApp;
friend GTreeCtrl;
afx_msg void OnNMClickTabTarget(NMHDR *pNMHDR, LRESULT *pResult);
};

extern CGiristDlg* g_pFrame;

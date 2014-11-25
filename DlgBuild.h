#pragma once
#include "afxwin.h"
#include "./util/gip_thread.h"
#include "afxcmn.h"


// CDlgBuild dialog

class CDlgBuild : public CDialog
{
	DECLARE_DYNAMIC(CDlgBuild)
	HANDLE m_hThread;
	GU_CTRL_DATA m_Progress;
	GIP_LIB *m_hLib;

public:
	CDlgBuild(GIP_LIB *hLib,CWnd* pParent = NULL);   // standard constructor
	virtual ~CDlgBuild();

	virtual void OnCancel();
	virtual void OnOK();

// Dialog Data
	enum { IDD = IDD_BUILD_DLG };

protected:
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	int StartThread( );
//	int StopThread( );
	void CtrlOnState( BOOL state );

	DECLARE_MESSAGE_MAP()
	afx_msg LRESULT OnUserBarUpdate(UINT wParam, LONG lParam); 
	afx_msg LRESULT OnUserTextUpdate(UINT wParam, LONG lParam); 
	afx_msg LRESULT OnThreadCompleted(UINT wParam, LONG lParam); 

public:
	afx_msg void OnBnClickedBtSwitch();
protected:
	CProgressCtrl m_cProgress;
	CListBox m_lbInfo;
private:
	CButton m_btSwitch;
	CButton m_btDistri;
	CButton m_btCancel;
	CButton m_btOK;
public:
	afx_msg void OnBnClickedBtDistri();
};

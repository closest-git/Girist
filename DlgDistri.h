#pragma once
#include "boardcontrol.h"
#include "afxwin.h"
#include "./util/gip_thread.h"
#include "afxcmn.h"


// CDlgDistri dialog

class CDlgDistri : public CDialog
{
	DECLARE_DYNAMIC(CDlgDistri)
	HANDLE m_hThread;
	GIP_LIB *m_hLib;
	GU_CTRL_DATA m_Progress;

	int StartThread( );
	void CtrlOnState( BOOL bAlive )	;
public:
	CDlgDistri(GIP_LIB *hLib,CWnd* pParent = NULL);   // standard constructor
	virtual ~CDlgDistri();

	virtual void OnCancel();
	virtual void OnOK();

// Dialog Data
	enum { IDD = IDD_DISTRI_DLG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	virtual BOOL OnInitDialog();

	DECLARE_MESSAGE_MAP()
	afx_msg LRESULT OnUserBarUpdate(UINT wParam, LONG lParam); 
	afx_msg LRESULT OnUserTextUpdate(UINT wParam, LONG lParam); 
	afx_msg LRESULT OnThreadCompleted(UINT wParam, LONG lParam); 

public:
	CDistriControl m_cDistri;
	CROCControl m_cROC;
protected:
	CEdit m_edDrimi;
	CComboBox m_cbLevel;
public:
	afx_msg void OnBnClickedBtDistri();
	CProgressCtrl m_cProgress;	
protected:
//	CEdit m_edWaveLen;
//	CEdit m_edShift;
	int m_Shift;
	int m_WaveLen;
	double m_Sep;
	CButton m_btOK;
	CButton m_btCancel;
	CButton m_btSwitch;
	CStatic m_stSep;
public:
	afx_msg void OnEnChangeEditSep();
};

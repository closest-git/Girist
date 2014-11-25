#pragma once
#include "afxwin.h"


// CDlgProgress dialog

class CDlgProgress : public CDialog
{
	DECLARE_DYNAMIC(CDlgProgress)

public:
	CDlgProgress(CWnd* pParent = NULL);   // standard constructor
	virtual ~CDlgProgress();
	BOOL Create();

	virtual void OnCancel();
	virtual void OnOK();
	virtual void PostNcDestroy();

// Dialog Data
	enum { IDD = IDD_PROGRESS_DLG };

protected:
	virtual BOOL OnInitDialog();
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	CListBox m_lbInfo;
};

extern CDlgProgress *g_pProgress;

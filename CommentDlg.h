#pragma once


// CCommentDlg dialog

class CCommentDlg : public CDialog
{
	DECLARE_DYNAMIC(CCommentDlg)

public:
	CCommentDlg(CWnd* pParent = NULL);   // standard constructor
	virtual ~CCommentDlg();

// Dialog Data
	enum { IDD = IDD_COMMENT };
	CString m_sComment;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()

};

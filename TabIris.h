#pragma once

#include "TabPageSSL.h"
#include "BoardControl.h"
#include "./util/gip_library.h"
#include "afxwin.h"

typedef enum {
	TAB_IRIS_SOURCE,TAB_IRIS_TARGET
}TAB_IRIS_TYPE;

// CTabIris dialog
class CBoardControl;

class CTabIris : public CTabPageSSL
{
	DECLARE_DYNAMIC(CTabIris)

public:
	CTabIris(CWnd* pParent = NULL);   // standard constructor
	virtual ~CTabIris();

// Dialog Data
	enum { IDD = IDD_TAB_IRIS };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	virtual BOOL OnInitDialog();

	DECLARE_MESSAGE_MAP()
	GIP_LIB *m_hLib;
public:
	afx_msg void OnStnClickedStaticNormal();
	void Update( int flag );
	void Clear( int flag );

	TAB_IRIS_TYPE type;
	CEyeControl m_ctrlEye;
	CBoardControl m_ctrlNormal;
//	GIP_LIB_entry *m_hEntry;
	CStatic m_stInfo;		//¹²6ÐÐ
	CEdit m_edNotes;
	afx_msg void OnBnClickedButton_1();
	afx_msg void OnBnClickedButton_2();
	CButton m_BT_1;
	CButton m_BT_2;
};

// DlgProgress.cpp : implementation file
//

#include "stdafx.h"
#include "Girist.h"
#include "DlgProgress.h"

CDlgProgress *g_pProgress;

// CDlgProgress dialog

IMPLEMENT_DYNAMIC(CDlgProgress, CDialog)

CDlgProgress::CDlgProgress(CWnd* pParent /*=NULL*/)
	: CDialog(CDlgProgress::IDD, pParent)
{

}

CDlgProgress::~CDlgProgress()
{
}

void CDlgProgress::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_LIST2, m_lbInfo);
}

/*
	v0.1	cys
		3/4/2009
*/
void CDlgProgress::PostNcDestroy(){
//	delete this;

	g_pProgress = GIP_NULL;
}

/*
	v0.1	cys
		3/4/2009
*/
void CDlgProgress::OnCancel()	{
	DestroyWindow();
}

/*
	v0.1	cys
		3/4/2009
*/
void CDlgProgress::OnOK()	{
	CDialog::OnOK( );
//		DestroyWindow();
}

BEGIN_MESSAGE_MAP(CDlgProgress, CDialog)
END_MESSAGE_MAP()

/*
	v0.1	cys
		3/4/2009		
*/
BOOL CDlgProgress::Create(){
	BOOL ret = CDialog::Create(CDlgProgress::IDD);
	return ret;
}

/*
	v0.1	cys
		3/4/2009
*/
BOOL CDlgProgress::OnInitDialog()
{
	CDialog::OnInitDialog();

/*	CString str( _T("Hello World") );
	m_lbInfo.AddString(str);
	m_lbInfo.AddString(str);
	m_lbInfo.AddString(str);*/

	g_pProgress = this;
	return TRUE;  // return TRUE  unless you set the focus to a control
}
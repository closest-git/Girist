// CommentDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Girist.h"
#include "CommentDlg.h"


// CCommentDlg dialog

IMPLEMENT_DYNAMIC(CCommentDlg, CDialog)

CCommentDlg::CCommentDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CCommentDlg::IDD, pParent)
{

}

CCommentDlg::~CCommentDlg()
{
}

void CCommentDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);

	DDX_Text( pDX,IDC_ED_COMMENT,m_sComment );

}


BEGIN_MESSAGE_MAP(CCommentDlg, CDialog)
END_MESSAGE_MAP()


// CCommentDlg message handlers

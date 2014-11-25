// TabTree.cpp : implementation file
//

#include "stdafx.h"
#include "Girist.h"
#include "TabTree.h"


// TAB_TREE dialog

IMPLEMENT_DYNAMIC(TAB_TREE, CDialog)

TAB_TREE::TAB_TREE(CWnd* pParent /*=NULL*/)
	: CTabPageSSL(TAB_TREE::IDD, pParent)
{

}

TAB_TREE::~TAB_TREE()
{
}

void TAB_TREE::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_TREE1, m_tree);
}


BEGIN_MESSAGE_MAP(TAB_TREE, CDialog)
END_MESSAGE_MAP()


// TAB_TREE message handlers

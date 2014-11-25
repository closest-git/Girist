#pragma once
#include "TabPageSSL.h"
#include "gtreectrl.h"

class TAB_TREE : public CTabPageSSL
{
	DECLARE_DYNAMIC(TAB_TREE)

public:
	TAB_TREE(CWnd* pParent = NULL);   // standard constructor
	virtual ~TAB_TREE();

// Dialog Data
	enum { IDD = IDD_TAB_TREE };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	GTreeCtrl m_tree;
};

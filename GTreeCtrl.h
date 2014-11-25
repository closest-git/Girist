#pragma once


// GTreeCtrl
class CGiristDlg;

class GTreeCtrl : public CTreeCtrl
{
	DECLARE_DYNAMIC(GTreeCtrl)

	GIP_LIB *m_hLib;
	int m_oldNo;
public:
	GTreeCtrl( );
	virtual ~GTreeCtrl();

	void OnUpdate( G_ENTRY_STATUS status );
	void			SetItemFont(HTREEITEM, LOGFONT&);
	void			SetItemBold(HTREEITEM, BOOL);
	void			SetItemColor(HTREEITEM, COLORREF);
	BOOL			GetItemFont(HTREEITEM, LOGFONT *);
	BOOL			GetItemBold(HTREEITEM);
	COLORREF		GetItemColor(HTREEITEM);

protected:
	struct Color_Font {
		COLORREF color;
		LOGFONT  logfont;
	};
	CMap <void*, void*, Color_Font, Color_Font&> m_mapColorFont;

	int m_isSort;
	//{{AFX_MSG(GTreeCtrl)
	afx_msg void OnLButtonDblClk(UINT , CPoint );
	afx_msg void OnSelchanged(NMHDR* pNMHDR, LRESULT* pResult);
	afx_msg void OnPaint();
	//}}AFX_MSG	
	DECLARE_MESSAGE_MAP()

friend CGiristDlg;
};



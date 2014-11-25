// Filename:     CBoardControl.h
// 
// Created:      2009/2/28
// Author:       CYS
// 
// Description:  CBoardControl is the iris working bord.
// 
// Dependencies: VC 7.1, STL, GDIPLUS
// 
// Credits: 
//
//////////////////////////////////////////////////////////////////////

#pragma once

#include <gdiplus.h>
using namespace::Gdiplus;
#include "./pch/gip_def.h"
#include "./util/gip_fv.h"
#include "./util/gip_library.h"

#define G_BOARD_CLASSNAME    _T("G_BOARD_CTRL")  // Window class name


//显示的类型	cys 6/1/2002	
// 6/28/2002	移至renjudef.h
typedef int SHOW_TYPE;
const SHOW_TYPE SHOW_NONE=0X0;
const SHOW_TYPE SHOW_T2=0X01;
const SHOW_TYPE SHOW_T3=0X02;
const SHOW_TYPE SHOW_T4=0X04;
const SHOW_TYPE SHOW_T43=0X08;

const SHOW_TYPE SHOW_SEGMENT=0X10;
const SHOW_TYPE SHOW_EVAL=0X20;

const SHOW_TYPE SHOW_BOOK=0x01000000;

const SHOW_TYPE SHOW_BLACK=0x80000000;
const SHOW_TYPE SHOW_WHITE=0x40000000;

#define SHOW_CODE( show ) (show&0xFFFF)
#define SHOW_BLACK(sh)	BIT_TEST(sh,SHOW_BLACK) 
#define SHOW_THREAT(sh)	BIT_IS(sh,0x07) 

//颜色
//cys 6/10/2002
const COLORREF crPick = RGB(0,255,0);
const COLORREF crWit = RGB(255,0,0);
const COLORREF crInfo = RGB(0,0,255);

class CTabIris;
class CGiristDlg;
class CBoardControl :public CWnd
{
	DECLARE_DYNAMIC(CBoardControl)

protected:
	const CString m_strClassName;
	bool isUnInitMenu,isSaveToJPG;

	DECLARE_MESSAGE_MAP()	

	void OnInitMenuPopup(CMenu* pPopupMenu, UINT nIndex, BOOL bSysMenu) ;
	void DrawShow( Graphics &graphics );
	BOOL GetPos( CPoint		pt,int& row,int& col );
	void GetPos( int row, int col,PointF& pt );
public:
	afx_msg void OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT , CPoint point);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnPaint();	
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnDestroy();
	afx_msg BOOL OnSetCursor( CWnd* pWnd,UINT nHitTest,UINT message );
	afx_msg LRESULT OnMouseHover(WPARAM wparam, LPARAM lparam);
	afx_msg LRESULT OnUnInitMenuPopup(WPARAM wparam, LPARAM lparam);
	afx_msg void OnBookSave();
	afx_msg void OnBookDelete();
	afx_msg void OnUpdateBookDelete(CCmdUI *pCmdUI);
	afx_msg void OnUpdateBookSave(CCmdUI *pCmdUI);

public:
	CBoardControl(void);
	~CBoardControl(void);
	BOOL Create(CWnd* pParentWnd, const RECT& rect, UINT nID, DWORD dwStyle = WS_VISIBLE|WS_TABSTOP);
	BOOL OnExit();

	int mode;
	/** DrawButton: draws button to the specified Graphics buffer*/
	virtual Bitmap *DrawBoard( );
	void DrawSignalEx( Graphics &graphics, int r,int c, int type,void * pData,ARGB cr=Color::Black );
/** CalculateDimensions: recalculates control region and drawing area */
	void ReCalculateDimensions();	
	void Update( );
	void OnMove( int r,int c );
	bool OpenDocument( CString sPath );
	void LoadObj( GIP_OBJECT *hObj,int flag );

protected:
	REAL b_grid,b_x0,b_y0;
	CString sNote;

	/** m_pCachedBitmap: for better performace, we'll use this for double-buffering */
	CachedBitmap* m_pCachedBitmap; 
	Bitmap* m_bmpObj;
	/** m_bDirtyBuffer: when set to true, indicates that CachedBitmap needs to be updated */
	BOOL m_bDirtyBuffer;
	/** m_ClientRect: stores the client area dimensions */
	Rect m_ClientRect;
	/** m_hRGN: specifies the shape of the control */	
	HRGN m_hRGN;
	/** m_bMouseOver: set to true then the mouse coursor move over the control */
	BOOL m_bMouseOver;
	/** m_bTracking: helps in establishing mouse-over-control event*/
	BOOL m_bTracking;
	SHOW_TYPE show;

friend CTabIris;
friend CGiristDlg;
};

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/3/2009	
*/
class CEyeControl :public CBoardControl
{
	DECLARE_DYNAMIC(CEyeControl)

protected:
	float p_radius,p_r,p_c,s_radius,s_r,s_c;

	DECLARE_MESSAGE_MAP()	

public:
	CEyeControl(void);
	~CEyeControl(void);

	Bitmap *DrawBoard( );

friend CTabIris;
friend CGiristDlg;
};

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/10/2009	
*/
class CDlgDistri;
class CDistriControl :public CBoardControl
{
	DECLARE_DYNAMIC(CDistriControl)
	GIP_LIB *hLib;
	G_DCRIMI_ *hDcimi;
protected:

	DECLARE_MESSAGE_MAP()	

public:
	CDistriControl(void);
	~CDistriControl(void);

	Bitmap *DrawBoard( );
	void Update( int flag );

friend CDlgDistri;
friend CGiristDlg;
};

class CROCControl :public CBoardControl
{
	DECLARE_DYNAMIC(CROCControl)
	GIP_LIB *hLib;
	G_DCRIMI_ *hDcimi;
protected:

	DECLARE_MESSAGE_MAP()	

public:
	CROCControl(void);
	~CROCControl(void);

	Bitmap *DrawBoard( );
	void Update( int flag );

friend CDlgDistri;
friend CGiristDlg;
};
// ButtonControl.cpp : implementation file

#include "StdAfx.h"
#include ".\boardcontrol.h"
#include ".\Girist.h"
#include ".\GiristDlg.h"
#include ".\pch\GIP_def.h"
#include ".\util\GIP_util.h"
#include ".\util\GIP_imag_transform.h"
#include "CommentDlg.h"
#include ".\util\2PassScale.h"

/*
	[ISSUE-NEEDED]
	1 基于GDI+华而不实，转为调用GDI与OpenCV
*/
const int BOARD_SIZE = 15;
const int CROSS_NUM = BOARD_SIZE*BOARD_SIZE;
const int ROW_NUM = BOARD_SIZE*6-2;

#define SIG_NONE	0x00
#define SIG_TOP		0x01
#define SIG_DOWN	0x02
#define SIG_FORBID	0x03
#define SIG_INFO	0x04
#define SIG_BORDER	0x05
#define SIG_IGNORE	0x06		//	6/10/2002 cys
#define SIG_BLACK	0x07
#define SIG_WHITE	0x08		
#define SIG_WEIGHT	0x09	
#define SIG_THREAT	0x10	
#define SIG_BOOK_PLY	0x11

IMPLEMENT_DYNAMIC(CBoardControl, CWnd)
IMPLEMENT_DYNAMIC(CEyeControl, CBoardControl)
IMPLEMENT_DYNAMIC(CDistriControl, CBoardControl)
IMPLEMENT_DYNAMIC(CROCControl, CBoardControl)

#define WM_UNINITMENUPOPUP 0x0125

BEGIN_MESSAGE_MAP(CBoardControl, CWnd)
	ON_WM_INITMENUPOPUP()
	ON_WM_PAINT()
	ON_WM_SIZE()
	ON_WM_ERASEBKGND()
	ON_WM_DESTROY()
	ON_WM_LBUTTONDOWN()
	ON_WM_RBUTTONDOWN()
	ON_WM_MOUSEMOVE()
	ON_WM_KEYUP()
	ON_WM_SETCURSOR()
	ON_MESSAGE(WM_MOUSEHOVER, OnMouseHover)
	ON_MESSAGE(WM_UNINITMENUPOPUP, OnUnInitMenuPopup)
//	ON_COMMAND(ID_BOOK_SAVE, OnBookSave)
//	ON_COMMAND(ID_BOOK_DELETE, OnBookDelete)
//	ON_UPDATE_COMMAND_UI(ID_BOOK_DELETE, OnUpdateBookDelete)
//	ON_UPDATE_COMMAND_UI(ID_BOOK_SAVE, &CBoardControl::OnUpdateBookSave)
END_MESSAGE_MAP()

BEGIN_MESSAGE_MAP(CEyeControl,CBoardControl)
END_MESSAGE_MAP()

BEGIN_MESSAGE_MAP(CDistriControl,CBoardControl)
END_MESSAGE_MAP()

BEGIN_MESSAGE_MAP(CROCControl,CBoardControl)
END_MESSAGE_MAP()

CBoardControl::CBoardControl(void) : m_strClassName(G_BOARD_CLASSNAME), 
	m_bDirtyBuffer(TRUE), m_hRGN(NULL),
	m_bMouseOver(FALSE), m_bTracking(FALSE),show(SHOW_NONE)
{

    WNDCLASS wndcls;
    HINSTANCE hInst = AfxGetInstanceHandle();

    if (!(::GetClassInfo(hInst, m_strClassName, &wndcls)))
    {
        // otherwise we need to register a new class
        wndcls.style            = CS_DBLCLKS | CS_HREDRAW | CS_VREDRAW;
        wndcls.lpfnWndProc      = ::DefWindowProc;
        wndcls.cbClsExtra       = wndcls.cbWndExtra = 0;
        wndcls.hInstance        = hInst;
        wndcls.hIcon            = NULL;
        wndcls.hCursor          = AfxGetApp()->LoadStandardCursor(IDC_ARROW);
        wndcls.hbrBackground    = (HBRUSH) (COLOR_3DFACE + 1);
        wndcls.lpszMenuName     = NULL;
        wndcls.lpszClassName    = m_strClassName;

        if (!AfxRegisterClass(&wndcls))
        {
            AfxThrowResourceException();
 //           return FALSE;
        }
    }

//	mode=MT_PLAY;
	isUnInitMenu = false;
	isSaveToJPG = false;

	m_bmpObj=GIP_NULL;
 //   return TRUE;
}

CBoardControl::~CBoardControl(void)
{		
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/3/2009
*/
CEyeControl::CEyeControl(void) {
	p_radius=-1.0,	p_r=-1.0,	p_c=-1.0;
	s_radius=-1.0,	s_r=-1.0,	s_c=-1.0;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/3/2009
*/
CEyeControl::~CEyeControl(void)
{		
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/10/2009
*/
CDistriControl::CDistriControl(void) {
	hLib = &g_Lib;
	hDcimi = &g_dcrimi;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/10/2009
*/
CDistriControl::~CDistriControl(void)
{		
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/1/2009
*/
CROCControl::~CROCControl(void)
{		
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/1/2009
*/
CROCControl::CROCControl(void) {
	hLib = &g_Lib;
	hDcimi = &g_dcrimi;
}



/*
	v0.1	cys
		4/6/2006
*/
BOOL CBoardControl::Create(CWnd* pParentWnd, const RECT& rect, UINT nID, DWORD dwStyle /*=WS_VISIBLE*/)
{
	return CWnd::Create( m_strClassName, _T(""), dwStyle, rect, pParentWnd, nID);
}

void CBoardControl::OnDestroy()
{
	if( m_bmpObj )	{
		m_bmpObj->DeleteObject( );
		delete m_bmpObj;
		m_bmpObj = NULL;
	}

	if (m_hRGN)
	{
		DeleteObject(m_hRGN);
		m_hRGN = NULL;
	}
}

void CBoardControl::OnSize(UINT, int, int) {
	m_bDirtyBuffer = TRUE;
}

/*
	避免闪烁

	http://www.west263.com/info/html/chengxusheji/C-C--/20080224/12018.html
*/
BOOL CBoardControl::OnEraseBkgnd(CDC*)
{
	// Return true so that the control is not erased
	// (Aside: erase means paint the surface with the window's default brush)
	return TRUE;
}

void CBoardControl::Update( ){
/*	switch( mode )	{
	case MT_PLAY:
		break;
	case MT_BOOK:	
		g_State.Init( &g_Board,true );
		g_Book.Analysis( &g_State );
#ifdef _GRUS_RENSTAR_EXPERT_
		g_pFrame->m_tabBook.OnPly( g_Board.plys[g_Board.curNo] );
#endif		
		break;
	default:
		break;
	}*/
	m_bDirtyBuffer = TRUE;	
	Invalidate(FALSE);		
}

LRESULT CBoardControl::OnMouseHover(WPARAM, LPARAM)
{
	m_bDirtyBuffer = TRUE;
	m_bMouseOver = TRUE;
	m_bDirtyBuffer = TRUE;
	

	Invalidate();	
	
	return 1;
}

void CBoardControl::OnMouseMove( UINT , CPoint point )
{ 
//	CWnd::OnMouseMove( nFlags,point );
//	this->SetFocus( );

	if (!m_bTracking)
	{
		TRACKMOUSEEVENT tme;
		tme.cbSize = sizeof(tme);
		tme.hwndTrack = m_hWnd;
		tme.dwFlags = TME_LEAVE|TME_HOVER;
		tme.dwHoverTime = 1;
		m_bTracking = _TrackMouseEvent(&tme);		
	}	
	if(m_bMouseOver == FALSE)
	{
		m_bMouseOver = TRUE;		
	}
/*
	CString info = _T("out of board"),data;
	int r,c,type;
	if( GetPos( point,r,c )	)	{
		type=g_Board.Type( r,c );
		if( type==ST_EMPTY )	{
			info = GBoard::CROSS_NAME( r,c );
		}
		else
			info.Format( _T("ply-%d:%s"),g_Board.No(r,c)+1,GBoard::CROSS_NAME( r,c ) );
	}
*///	g_pFrame->SetSBInfo( 0,info );
}

/*
	v0.1	cys
		3/18/2006
*/
void CBoardControl::OnMove( int r,int c )	{
	CString info = _T("NO PLY from ENGINE!");

//	g_Board.AddMove( r,c );
	show=SHOW_NONE;
/*	switch( mode )	{
	case MT_PLAY:
		g_pFrame->SetSBInfo( 1,_T("Begin Search...") );
		g_State.Init( &g_Board,true );
		if( true && g_State.numPly>0 )	
			g_hSearch->CheckPlyIsPlaned( &g_State );
		g_hSearch->Init( &g_State,&g_S_Dup,NULL );
		g_hSearch->Search( 0x0 );
		g_pFrame->UpdateABTree( );
		g_pFrame->m_tabTree.m_treeAB.nCut4Save = 0;
		if( g_hSearch->nPly > 0 )	{
			G_POS pos=g_hSearch->plys[0];
			g_Board.AddMove( G_ROW( pos ),G_COL( pos ) );
			g_pFrame->m_tabTree.m_treeAB.nCut4Save = 1;
		}
		g_pFrame->SetSBInfo( 1,g_hSearch->m_sDump );
		break;
	case MT_BOOK:
		if( g_Book.GetCurNode( ) != NULL )	{
			g_pFrame->SetSBInfo( 1,_T("Find results from GBook") );
		}
		break;
	default:
		break;
	}*/
	Update( );

}

/*
	v0.2	cys	
		4/22/2006
*/
void CBoardControl::OnLButtonDown(UINT nFlags, CPoint point )	{	
	if( isUnInitMenu )	{	
		isUnInitMenu=false;
		return CWnd::OnLButtonDown( nFlags,point );	
	}
	m_bDirtyBuffer = TRUE;

	Invalidate();	
}

void CBoardControl::ReCalculateDimensions()
{
	/***********************************************************

		NOTE from the author:

		GDI regions and GDI+ ROUND regions are NOT the same!
		That's why I am going through extra effort to create 
		the region using GDI+ in order for the control to 
		render properly around the edges.

	************************************************************/

	// Calculate CDialog offset
	POINT ptParent = {0, 0};
	::ClientToScreen(m_hWnd, &ptParent);
	::ScreenToClient(GetParent()->m_hWnd, &ptParent);
	// Obtain client rectangle
	GetClientRect(&m_ClientRect);
	int width=m_ClientRect.Width( ),height=m_ClientRect.Height( );
	if( width>0 && height>0 )	{
//		CDC memDC;	
		CPaintDC dc(this);
//		memDC.CreateCompatibleDC(NULL);	 						//建立和屏幕显示兼容的内存显示设备
		m_bmpObj = new CBitmap( );
		m_bmpObj->CreateCompatibleBitmap(&dc,width,height);		//建立一个和屏幕显示兼容的位图
		LoadObj( GIP_NULL,0x0 );
//		memDC.SelectObject(m_bmpObj);
//		memDC.FillSolidRect(0,0,width,height,RGB(255,255,255));
	}

/*	// Get the device context & create Graphics class
	Graphics graphics(m_hWnd);

	// Draw the region
	GraphicsPath path;
	path.AddRectangle(m_ClientRect);
//	path.AddEllipse(m_ClientRect);
	Region rgn(&path);
	SetWindowRgn(rgn.GetHRGN(&graphics), FALSE);

	// Apply CDialog offset
	rgn.Translate(ptParent.x, ptParent.y);
	m_hRGN = rgn.GetHRGN(&graphics);
*/
}

/*
	基于GDI的双缓存	

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/3/2009

*/
void CBoardControl::OnPaint()
{
	// device context for painting
	int m_rW = m_ClientRect.Width( ),m_rH =m_ClientRect.Height( );

	CPaintDC dc(this);
	CDC MemDC;		
	MemDC.CreateCompatibleDC(NULL);	 						//建立和屏幕显示兼容的内存显示设备
//	CBitmap MemBitmap;
//	MemBitmap.CreateCompatibleBitmap(&dc,m_rW,m_rH);		//建立一个和屏幕显示兼容的位图
//	CBitmap *pOldBit=MemDC.SelectObject(&MemBitmap);		//将位图选入到内存显示设备中
//	MemDC.FillSolidRect(0,0,m_rW,m_rH,RGB(255,255,255));
	if( m_bmpObj==NULL )	{
		CWnd::OnPaint( );		//必须调用CWnd::OnPaint( )，否则消息传递出现问题
	}else	{
		BITMAP tagBitMap;
		m_bmpObj->GetBitmap( &tagBitMap );
		int width=tagBitMap.bmWidth,			height=tagBitMap.bmHeight;

		MemDC.SelectObject(m_bmpObj);
		if( width==m_rW && height==m_rH )	{
			dc.BitBlt(0,0,m_rW,m_rH,&MemDC,0,0,SRCCOPY);	
		}else	{
			dc.SetStretchBltMode(HALFTONE);
			dc.SetStretchBltMode(COLORONCOLOR);
			dc.StretchBlt(0,0,m_rW,m_rH,&MemDC,0,0,width,height,SRCCOPY);
		}
	}
//	MemBitmap.DeleteObject();
	MemDC.DeleteDC();

	return;
}

/*
	ISSUE-NEEDED
		1	DASH-DOT PEN

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/5/2009
*/
void CEyeControl::GIP_OBJ_BMP( GIP_OBJECT *hObj,COLORREF *pBMP,int flag  )	{
	if( hObj==GIP_NULL )
		return;
	else	{
		CBoardControl::GIP_OBJ_BMP( hObj,pBMP,flag  );
		GIP_FLOATING *data=hObj->data;
		int *arrPos=(int*)GIP_alloc( sizeof(int)*hObj->M*hObj->N ),nz=-1,i;

		nz = GIP_circle_pos( hObj,p_r,p_c,p_radius,arrPos,0x0 );
		for( i = 0; i < nz; i++ )
			pBMP[arrPos[i]]=RGB( 255,0,0 );
		nz = GIP_circle_pos( hObj,s_r,s_c,s_radius,arrPos,0x0 );
		for( i = 0; i < nz; i++ )
			pBMP[arrPos[i]]=RGB( 255,0,0 );
		GIP_free( arrPos );
	}
}


#ifdef _GRUS_RENSTAR_EXPERT_
/*	
	v0.1	cys
		4/21/2006

	copy form ms-help://MS.VSCC.v80/MS.MSDN.v80/MS.WIN32COM.v10.en/gdicpp/GDIPlus/usingGDIPlus/usingimageencodersanddecoders/retrievingthe.htm
*/
int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)	{
   UINT  num = 0;          // number of image encoders
   UINT  size = 0;         // size of the image encoder array in bytes

   ImageCodecInfo* pImageCodecInfo = NULL;

   GetImageEncodersSize(&num, &size);
   if(size == 0)
      return -1;  // Failure

   pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
   if(pImageCodecInfo == NULL)
      return -1;  // Failure

   GetImageEncoders(num, size, pImageCodecInfo);

   for(UINT j = 0; j < num; ++j)
   {
      if( wcscmp(pImageCodecInfo[j].MimeType, format) == 0 )
      {
         *pClsid = pImageCodecInfo[j].Clsid;
         free(pImageCodecInfo);
         return j;  // Success
      }    
   }

   free(pImageCodecInfo);
   return -1;  // Failure
}
#endif

/*
	v0.1 cys		
		3/27/2006	
*/
void CBoardControl::OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags)	{
	unsigned char kbuf[256];
	VERIFY( GetKeyboardState( kbuf ) != 0 );
	bool isBlack= (GetKeyState(VK_CONTROL)&128) == 128 ;

	m_bDirtyBuffer = TRUE;	
	switch( nChar )	{
#ifdef _GRUS_RENSTAR_EXPERT_
	case 'J':		case 'j':	//save board to jpg file
		{
		CFileDialog	dlg(FALSE, _T("*.jpg"), NULL, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, 
			_T("JPEG文件交换格式(*.jpg)|*.jpg"), this);
		if( dlg.DoModal( ) != IDOK )	
			break;	

		CLSID pngClsid;
		if( GetEncoderClsid(L"image/png", &pngClsid) != -1 )	{
			isSaveToJPG = true;		
			CBitmap	*pBmp=DrawBoard( );	
/*			Graphics *pGraphics = Graphics::FromImage(pBmp);
			pGraphics->DrawCachedBitmap(m_pCachedBitmap, 0, 0);*/
			CString sPath = dlg.GetPathName( );
			pBmp->Save( sPath, &pngClsid, NULL);
			isSaveToJPG = false;
//			delete pGraphics;			
			delete pBmp;
		}
		}		
		break;
	case 'E':		case 'e':			//evaluate
/*		g_State.Init( &g_Board,true );
		g_hSearch->Init( &g_State,&g_S_Dup,NULL );
		g_hSearch->Search( 0x0 );
		g_pFrame->UpdateABTree( );*/
		isBlack = g_State.stCur==ST_BLACK;
		g_hSearch->GetEval( isBlack?gs_eval:NULL,isBlack?NULL:gs_eval );
		show=SHOW_EVAL;
		BIT_SET( show,isBlack ? SHOW_WHITE:SHOW_BLACK );
		break;
	case 'B':		case 'b':
		g_State.Init( &g_Board,true );
		g_State.Shape_Analysis( g_State.stCur,SHAPE_ANLYSIS_ROUTINE );
		g_SEG.Init( G_SEG::S_BLACK|G_SEG::S_THRESH,2 );	/*|G_SEG::S_LEVEL*/
		g_SEG.Grow_CYS( &g_State );
		show=SHOW_SEGMENT;
		BIT_SET( show,isBlack ? SHOW_WHITE:SHOW_BLACK );
		g_pFrame->SetSBInfo( 1,g_SEG.sInfo );
		break;
	case '4':		case '3':		
	case '2':		case '7':
		g_State.Init( &g_Board,true );
		g_Stat.Init( 0x0 );
		if( false )	{
			nOld = g_hSearch->ITER_DEPTH;
			g_hSearch->ITER_DEPTH = 1;
			g_hSearch->Init( &g_State,&g_S_Dup,NULL );
			g_hSearch->Search( 0x0 );
			g_hSearch->ITER_DEPTH = nOld;
			g_pFrame->UpdateABTree( );
		}
		g_State.Shape_Analysis( g_State.stCur,SHAPE_ANLYSIS_ROUTINE );
		//redraw the board		
		show=SHOW_NONE;
		BIT_SET( show,nChar=='7'? SHOW_T43 : nChar=='4'?SHOW_T4 : nChar=='3'?SHOW_T3:SHOW_T2 );
		BIT_SET( show,isBlack ? SHOW_WHITE:SHOW_BLACK );
		break;
#endif
	case VK_F5:
		show=SHOW_NONE;
		break;
	case VK_SHIFT:
//		isOnlineEdit=!isOnlineEdit;
		break;
	case VK_F2:
//		mode=mode==MT_BOOK ? MT_PLAY:MT_BOOK;
//		g_pFrame->m_cbModes.SetCurSel( mode );
//		g_pFrame->m_TAB.ActivateSSLPage( mode==MT_BOOK ? 0:1 );
//		g_pFrame->m_TAB.SetCurSel( mode==MT_BOOK ? 0:1 );
		break;
	case VK_LEFT:
//		if( g_Board.curNo > 0 )
//			g_Board.curNo--;
		break;
	default:
		m_bDirtyBuffer = FALSE;	
		break;
	}
	CWnd::OnKeyUp(nChar, nRepCnt, nFlags);

	if(	m_bDirtyBuffer )		
		Invalidate(FALSE);	
}

BOOL CBoardControl::OnSetCursor(CWnd* pWnd, UINT nHitTest, UINT message ) { 
/*	if( isOnlineEdit )
		::SetCursor( ::LoadCursor(NULL, IDC_HELP) );
	else
		::SetCursor( ::LoadCursor(NULL, IDC_ARROW) );*/

	return CWnd::OnSetCursor(pWnd, nHitTest, message);;
/*	if(nHitTest==HTCLIENT){ 
		SetCursor(m_hNowCursor); 
		return true; 
	} else 
		return CScrollView::OnSetCursor(pWnd, nHitTest, message); */
} 

/*
	v0.1	cys
		4/9/2006
*/
void CBoardControl::OnRButtonDown(UINT nFlags, CPoint point) {
#ifdef _GRUS_RENSTAR_EXPERT_
	// TODO: Add your message handler code here and/or call default
	if( mode==MT_BOOK )	{
		CMenu menu ;
		menu.LoadMenu( IDR_MENU_BOOK );
		ClientToScreen(&point);
		menu.GetSubMenu( 0 )->TrackPopupMenu( TPM_LEFTALIGN | TPM_RIGHTBUTTON, point.x , point.y ,this );
	}

#endif	
	CWnd::OnRButtonDown(nFlags, point);
}

#ifdef _GRUS_RENSTAR_EXPERT_

/*
	v0.1	cys
		4/9/2006
*/
void CBoardControl::OnBookSave(){
	// TODO: Add your command handler code here
	CString sPath =  g_Book.m_sPath;
	if( sPath==_T("") )	{
		CFileDialog dlg( FALSE,NULL,NULL,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
			_T("RenStar Book(*.book)|*.book|All Files(*.*)|*.*|") );
		if( dlg.DoModal( ) != IDOK )	
			return;
		sPath = dlg.GetPathName( ); 
	}
	g_Book.Insert( &g_Board );
	g_Book.Save( sPath );

	Update( );
}

/*
	
*/
void CBoardControl::OnBookDelete()
{
	// TODO: Add your command handler code here
	HB_NODE *hNode = g_Book.GetCurNode( );
	if( mode == MT_BOOK && hNode!= NULL )	{
		g_Board.OnDelete( G_GETPOS( hNode->ply) );
		g_Book.UpdateTree( hNode,HB_NODE::PLY_DELETE,0x0 );

		Update( );
	}
}

/*
	v0.1	cys
		4/9/2006
*/
void CBoardControl::OnUpdateBookDelete(CCmdUI *pCmdUI)
{
	// TODO: Add your command update UI handler code here
	ASSERT( mode==MT_BOOK );
	HB_NODE *hNode = g_Book.GetCurNode( );
	pCmdUI->Enable( hNode != NULL );
}

void CBoardControl::OnUpdateBookSave(CCmdUI *pCmdUI){
	// TODO: Add your command update UI handler code here
	ASSERT( mode==MT_BOOK );
	HB_NODE *hNode = g_Book.GetCurNode( );
	pCmdUI->Enable( (hNode != NULL && g_Book.isModify) || hNode==NULL );
}

#endif
/*
	v0.1	cys
		4/10/2006
*/
BOOL CBoardControl::OnExit( )	{
#ifdef _GRUS_RENSTAR_EXPERT_
	int nCode;
	if( g_Book.isModify )	{
		nCode=AfxMessageBox( _T("RENSTAR_BOOK has been modified,Save it?"),MB_YESNOCANCEL );
		if( nCode==IDCANCEL )
			return FALSE;
		if( nCode==IDYES )	{
			g_Book.Save( g_Book.m_sPath );
		}
	}else if( mode==MT_BOOK && g_Book.round.cur==GBook::NULL_NODE )	{
		nCode=AfxMessageBox( _T("There are new plys,Add it to RENSTAR_BOOK?"),MB_YESNOCANCEL );
		if( nCode==IDCANCEL )
			return FALSE;
		if( nCode==IDYES )	{
			g_Book.Insert( &g_Board );
			g_Book.Save( g_Book.m_sPath );
		}
	}
#endif
	return TRUE;
}

void CBoardControl::OnInitMenuPopup(CMenu* pPopupMenu, UINT nIndex, BOOL bSysMenu) {
  CCmdUI    CmdUI;

  // Skip system menu from processing
  if( FALSE != bSysMenu )
    return;

  CmdUI.m_nIndexMax = pPopupMenu->GetMenuItemCount();
  for( nIndex = 0; nIndex < CmdUI.m_nIndexMax; ++nIndex )
   {
     CmdUI.m_nIndex = nIndex;
     CmdUI.m_nID = pPopupMenu->GetMenuItemID(nIndex);
     CmdUI.m_pMenu = pPopupMenu;

  // There are two options. You can choose the one by removing comments
  // Option 1. All handlers are in dialog
     CmdUI.DoUpdate( this, TRUE );

  // Option 2. There are handlers in dialog and controls
/*
     CmdUI.DoUpdate( this, FALSE );
     // If dialog handler doesn't change state route update
     // request to child controls. The last DoUpdate will
     // disable menu item with no handler
     if( FALSE == CmdUI.m_bEnableChanged )
       CmdUI.DoUpdate( m_pControl_1, FALSE );
     ...
     if( FALSE == CmdUI.m_bEnableChanged )
       CmdUI.DoUpdate( m_pControl_Last, TRUE );
*/
   } 	    
}

LRESULT CBoardControl::OnUnInitMenuPopup(WPARAM, LPARAM){
	isUnInitMenu=true;
	return 1;
}

/*
	v0.1	cys
		4/17/2006
*/
bool CBoardControl::OpenDocument( CString sPath )	{
	CString sExt=GetFileExt( sPath ); 
	if( sExt=="ju" )	{
//		g_Board.Open( sPath );
	}else if( sExt=="lib" )	{
//		g_Book.Init( );
//		g_Book.Load_RenLib( sPath );
	}else if( sExt=="book" )		{
//		g_Book.Init( );
//		g_Book.Open( sPath );
	}else
		return false;

	switch( mode )	{
	case 0:
//		g_hSearch->Init( NULL,NULL,NULL );
		break;
	case 1:
		break;
	default:
		break;
	}	
	Update( );	
	return true;
}

/*
	GIP_OBJECT=>BMP(子类有不同的派生)
	假定灰度图像

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/4/2009
*/
void CBoardControl::GIP_OBJ_BMP( GIP_OBJECT *hObj,COLORREF *pBMP,int flag  )	{
	int x,y,val,pos,M=hObj->M,N=hObj->N;

	for ( x = 0; x < N; x++ )         {
	for ( y = 0; y < M; y++ )            {
		pos = GIP_RC2POS( y,x,N );
		val = G_DOUBLE2INT(hObj->data[pos]);
//		ASSERT( val>=0 && val<=255 );
		if( val<0 )	{
			val=0;
		}else if( val >255 )
			val=255;
		else	{	//0<val<=255		
		}
		pBMP[pos]=RGB(val,val,val);
	}
	}
}

extern "C" 	void GIP_cv_resize( GIP_OBJECT *hObj_0,GIP_OBJECT *hObj,int type,int flag );

/*
	更新m_bmpObj

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/3/2009
*/
void CBoardControl::LoadObj( GIP_OBJECT *hObj_0,int flag )	{
	CBitmap *pBmp=m_bmpObj;
	BITMAP tagBitMap;
	int width=-1,height=-1;
	CDC memDC;
	GIP_OBJECT *hObj=hObj_0;

	ASSERT( pBmp!=0x0 );
	memDC.CreateCompatibleDC(NULL);		
	pBmp->GetBitmap( &tagBitMap );
	width=tagBitMap.bmWidth,			height=tagBitMap.bmHeight;

	if( hObj==GIP_NULL )	{
		memDC.SelectObject(pBmp);
	//	memDC.FillSolidRect(0,0,width,height,RGB(255,255,255)); 
		memDC.FillSolidRect(0,0,width,height,GetSysColor(COLOR_BTNFACE) ); 
			
		m_bDirtyBuffer=TRUE;	
	}else	{
		int x,y;
		COLORREF *pOld,*pNew;
		pOld=new COLORREF[hObj_0->M*hObj_0->N];
		GIP_OBJ_BMP( hObj_0,pOld,0x0 );
		pNew = pOld;
		if( width!=hObj_0->N || height!=hObj_0->M )	{
			C2PassScale<CHammingFilter> ScaleEngine;
	//		C2PassScale<CBlackmanFilter> ScaleEngine;
			pNew = ScaleEngine.AllocAndScale( pNew,hObj_0->N,hObj_0->M,width,height );
			
/*			hObj=GIP_OBJ_create( height,width,GIP_NULL,0x0 );
			pNew=new COLORREF[width*height];
			GIP_cv_resize( hObj_0,hObj,0x0,0x0 );
			GIP_OBJ_BMP( hObj,pNew,0x0 );*/
			delete[] pOld;
		}
		
		memDC.SelectObject(pBmp);
		for ( x = 0; x < width; x++ )         {
		for ( y = 0; y < height; y++ )            {
			memDC.SetPixel( x, y, pNew[y*width+x] );
		}
		}
		delete[] pNew;

		m_bDirtyBuffer=TRUE;
	}
	memDC.DeleteDC( );

	if( hObj!=hObj_0 )		{
		GIP_OBJ_clear(hObj);		GIP_free( hObj );		
	}

	if(	m_bDirtyBuffer )		
		Invalidate(FALSE);	

	return ;
}


/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/10/2009
*/
void CDistriControl::Update(int flag )	{
	int m_rW=-1,m_rH=-1,i,r,c,val,n_grid=15,i_1;
	int D_span=hLib->D_span,h_1,h_2,h_3,h_4,c_0,c_1,peak_1,peak_2,pos_1,pos_2;
	float *D_inter=hLib->D_inter,*D_intra=hLib->D_intra,nz_a,nz_r;
	double D_s,inter_s,intra_s,sep=g_dcrimi.sep;

	VERIFY( hDcimi->hBase==hLib );
	nz_a = hDcimi->nz_a;		nz_r=hDcimi->nz_r;
	if( nz_a<=0.0 && nz_r<=0.0 )
		return;

	CBitmap *pBmp=m_bmpObj;
	BITMAP tagBitMap;
	CDC memDC;

	ASSERT( pBmp!=0x0 );
	memDC.CreateCompatibleDC(NULL);	
	pBmp->GetBitmap( &tagBitMap );
	m_rW=tagBitMap.bmWidth,			m_rH=tagBitMap.bmHeight;
	inter_s=0.0,		intra_s=0.0;
	for( c = n_grid*2; c < m_rW-n_grid; c++ )	{
		if( c==65 )
			i=0;
		c_0=G_DOUBLE2INT( (c-n_grid*2)*1.0/(m_rW-n_grid*3)*D_span );			
		c_1=G_DOUBLE2INT( (c+1-n_grid*2)*1.0/(m_rW-n_grid*3)*D_span );
		h_1=0.0;		h_2=0.0;
		for( i = c_0; i< c_1; i++ )	{
			h_1+=D_inter[i];		h_2+=D_intra[i];
		}
		inter_s=MAX(h_1/nz_r,inter_s);		
		intra_s=MAX(h_2/nz_a,intra_s);
	}
	ASSERT( inter_s>0 || intra_s>0 );
	D_s=(m_rH-n_grid*2)/MAX( inter_s,intra_s );

	memDC.SelectObject(pBmp);
	memDC.FillSolidRect(0,0,m_rW,m_rH,RGB(255,255,255));
	memDC.Rectangle( n_grid*2,n_grid,m_rW-n_grid,m_rH-n_grid );
//绘制distribution
	peak_1=0,	peak_2=0,	pos_1=-1,	pos_2=-1;
	for( c = n_grid*2; c < m_rW-n_grid; c++ )	{
		c_0=G_DOUBLE2INT( (c-n_grid*2)*1.0/(m_rW-n_grid*3)*D_span );			
		c_1=G_DOUBLE2INT( (c+1-n_grid*2)*1.0/(m_rW-n_grid*3)*D_span );
		h_1=0.0;		h_2=0.0;
		for( i = c_0; i< c_1; i++ )	{
			h_1+=D_inter[i];		h_2+=D_intra[i];
		}
		h_1 = nz_r==0 ? 0 : GIP_INTENSITY( h_1/nz_r*D_s );		ASSERT( h_1<=m_rH && h_1>= 0 );
		h_2 = nz_a==0 ? 0 : GIP_INTENSITY( h_2/nz_a*D_s );		ASSERT( h_2<=m_rH && h_2>= 0 );
		if( h_1==0 && h_2==0 )		continue;
		if( h_1>peak_1 )	{
			peak_1=h_1;		pos_1=c;
		}
		if( h_2>peak_2 )	{
			peak_2=h_2;		pos_2=c;
		}
		h_3=min(h_1,h_2);			h_4=max(h_1,h_2 );
		for( i = m_rH-h_3-n_grid; i < m_rH-n_grid; i++ )
			memDC.SetPixel( c, i, RGB( 255, 255, 0 ) );	
		COLORREF cr= h_1>h_2 ? RGB(255,0,0) : RGB(0,255,0 );
		for( i = m_rH-h_4-n_grid; i < m_rH-h_3-n_grid; i++ )	{
			memDC.SetPixel( c, i,cr );
		}
	}
	memDC.TextOutW( 0,0,_T("Hamming Distance Distributions") );
//Axis output
	CFont font_0;
	LOGFONT lf;
	memset(&lf, 0, sizeof(LOGFONT));       
	lf.lfHeight = 14;                      
	_tcscpy(lf.lfFaceName, _T("Arial") );        
	VERIFY(font_0.CreateFontIndirect(&lf));  
	CFont* pOldFont = memDC.SelectObject(&font_0);
	for( i = 1; i < 10; i++ )	{					//Hamming Distance Axis
		c = i*1.0/10*(m_rW-n_grid*3)+n_grid*2;
		sNote.Format( _T("%g"),i*1.0/10 );
		memDC.TextOutW( c-10,m_rH-n_grid,sNote );
		memDC.MoveTo( c,m_rH-n_grid );			memDC.LineTo( c,m_rH-n_grid-5 );
	}
	int nF = (int)(MAX( inter_s,intra_s )*100.0+1);	//Frequency Axis
	nF = MIN( nF,100 );
	i_1 = MIN(nF,10);
	for( i = 0; i < i_1; i++ )	{
		r = (m_rH-n_grid)-i*1.0/i_1*(m_rH-n_grid*2);
		sNote.Format( _T("%d%%"),G_DOUBLE2INT(i*1.0/i_1*nF) );
		memDC.TextOutW( 0,r-5,sNote );
		memDC.MoveTo( n_grid*2,r );			memDC.LineTo( n_grid*2+5,r );
	}
	memDC.SelectObject(pOldFont);
	font_0.DeleteObject();
	CPen pen( PS_DOT,0,RGB(0,0,0) ),*pOldPen;
	pOldPen = memDC.SelectObject( &pen );
	c = sep*(m_rW-n_grid*3)+n_grid*2;
	memDC.MoveTo( c,n_grid );			memDC.LineTo( c,m_rH-n_grid );
	memDC.SelectObject( &pOldPen );

	memDC.DeleteDC( );
//	GIP_save_IplImage_d( hObj->M,hObj->N,GIP_path_(_T("CHT_peak"),g_nzSample),hObj->data,hObj->N,0x0 );
	m_bDirtyBuffer=TRUE;
	if(	m_bDirtyBuffer )		
		Invalidate(FALSE);	

}

/*
	http://wwwiti.cs.uni-magdeburg.de/~sschimke/sose05/15-platasign/Evaluation_of_Biometric_Systems.html

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/1/2009
*/
void CROCControl::Update(int flag )	{
	int m_rW=-1,m_rH=-1,i,r,c,val,n_grid=15;
	int D_span=hLib->D_span,nz_a,nz_r,w_r,w_a,pos_1;
	float *D_inter=hLib->D_inter,*D_intra=hLib->D_intra,*f_ar,*f_rr,frr_1;
	double x_0=1.0e-7,y_0=1.0e-4,x_1,y_1,x,y,a;

	double D_s,inter_s,intra_s;

	CBitmap *pBmp=m_bmpObj;
	BITMAP tagBitMap;
	CDC memDC;

	ASSERT( pBmp!=0x0 );
	memDC.CreateCompatibleDC(NULL);		
	pBmp->GetBitmap( &tagBitMap );
	m_rW=tagBitMap.bmWidth,			m_rH=tagBitMap.bmHeight;
	memDC.SelectObject(pBmp);
	memDC.FillSolidRect(0,0,m_rW,m_rH,RGB(255,255,255));
	memDC.Rectangle( n_grid*2,n_grid,m_rW-n_grid,m_rH-n_grid );

	nz_a = hDcimi->nz_a;		nz_r=hDcimi->nz_r;
	if( nz_a<=0.0 || nz_r<=0.0 )
		return;
	f_ar=(float*)GIP_alloc( sizeof(float)*D_span*2 );		f_rr=f_ar+D_span;
	w_a=0,			w_r=0;
	frr_1=-1.0;		pos_1=-1;
	for( i = 0; i<D_span; i++ )	{
		w_a+=D_intra[i];		w_r+=D_inter[i];
		f_ar[i] = w_r*1.0/nz_r;
		f_rr[i] = (nz_a-w_a)*1.0/nz_a;
		if( f_ar[i]<x_0)	
			continue;
		if( pos_1==-1 )	{
			frr_1 = f_rr[i];	pos_1=i;
		}
//		frr_1=MAX( frr_1,f_rr[i] );
	}
	ASSERT( pos_1!=-1 );
	frr_1 = MIN(frr_1*1.2,0.9 );
	x_1 = log(1.0/x_0);			y_1 = log(frr_1/y_0);
//Axis output
	CFont font_0;
	LOGFONT lf;
	memset(&lf, 0, sizeof(LOGFONT));       
	lf.lfHeight = 12;                      
	_tcscpy(lf.lfFaceName, _T("Arial") );        
	VERIFY(font_0.CreateFontIndirect(&lf));  
	CFont* pOldFont = memDC.SelectObject(&font_0);
	CPen pen( PS_SOLID,0,RGB(100,100,100) ),*pOldPen;
	pOldPen = memDC.SelectObject( &pen );
//FAR axis
	CString sX[8]={
		_T("1E-5%"),_T("1E-4%"),_T("1E-3%"),_T("0.01%"),_T("0.1%"),_T("1%"),_T("10%"),_T("100%")
	};
	for( i = 0; i <= 7; i++ )	{
		c = i*1.0/7*(m_rW-n_grid*3)+n_grid*2;
		memDC.TextOutW( c-10,m_rH-n_grid,sX[i] );
		memDC.MoveTo( c,m_rH-n_grid );			memDC.LineTo( c,n_grid );
	}
//FRR axis
	int nF = 5;	//Frequency Axis
	CString sY[4]={
		_T("0.01%"),_T("0.1%"),_T("1%"),_T("10%")
	};
	a = 1.0e-4;
	for( i = 0; i <= 3; i++ )	{
		if( a>frr_1 )	{
			a = frr_1;
			sY[i].Format( _T("%4.2g%%"),frr_1*100.0 );
		}
		y=log( a/y_0);			
		r=(1-y/y_1)*(m_rH-n_grid*2)+n_grid;
		a*=10.0;
//		r = (m_rH-n_grid)-i*1.0/3*(m_rH-n_grid*2);
		memDC.TextOutW( 0,r-10,sY[i] );
		memDC.MoveTo( n_grid*2,r );			memDC.LineTo( m_rW-n_grid,r );
	}
	if( frr_1>0.15 )	{
		sNote.Format( _T("%4.2g%%"),frr_1*100.0 );
		memDC.TextOutW( 0,n_grid-10,sNote );
	}
	memDC.SelectObject(pOldFont);
	font_0.DeleteObject();
//ROC curve output
	CPen pen_1( PS_SOLID,0,RGB(255,0,0) );
	memDC.SelectObject( &pen_1 );
	r = -1;
	for( i = pos_1; i < D_span; i++ )	{
		ASSERT( f_ar[i]>=1.0e-7);
		if( f_rr[i]<y_0 )
			continue;
		x = log(f_ar[i]/x_0);						y=log( f_rr[i]/y_0);
		x = x/x_1*(m_rW-n_grid*3)+n_grid*2;			y=(1-y/y_1)*(m_rH-n_grid*2)+n_grid;
		if( r==-1 )
		{	memDC.MoveTo( x,y );		r=1;	}
		else
			memDC.LineTo( x,y );
		memDC.SetPixel( x, y, RGB( 255, 0, 0 ) );
	}
	memDC.SelectObject( pOldPen );
	x = log(g_dcrimi.rFAR/x_0);						y=log( g_dcrimi.rFRR/y_0);
	x = x/x_1*(m_rW-n_grid*3)+n_grid*2;				y=(1-y/y_1)*(m_rH-n_grid*2)+n_grid;
	memDC.Ellipse( x-3,y-3,x+3,y+3 );

	GIP_free( f_ar );

	memDC.TextOutW( m_rW/3,0,_T("ROC curves") );
	memDC.TextOutW( m_rW/2,m_rH-n_grid*2,_T("FAR") );
	memDC.TextOutW( 0,m_rH/2,_T("FRR") );
	memDC.DeleteDC( );

	m_bDirtyBuffer=TRUE;
	if(	m_bDirtyBuffer )		
		Invalidate(FALSE);	

}
// ABTreeCtrl.cpp : implementation file
//

#include "stdafx.h"
#include "Girist.h"
#include "GTreeCtrl.h"
#include "giristdlg.h"
#include "./util/gip_util.h"
#include "./girist/giris_core.h"

// GTreeCtrl

IMPLEMENT_DYNAMIC(GTreeCtrl, CTreeCtrl)

GTreeCtrl::GTreeCtrl( )
{
	m_hLib = &g_Lib;
	m_isSort = 0;
	m_oldNo = -1;
}

GTreeCtrl::~GTreeCtrl()
{

}


BEGIN_MESSAGE_MAP(GTreeCtrl, CTreeCtrl)
	ON_WM_PAINT()
	ON_WM_LBUTTONDBLCLK()
	ON_NOTIFY_REFLECT(TVN_SELCHANGED, OnSelchanged)
END_MESSAGE_MAP()

void GTreeCtrl::SetItemFont(HTREEITEM hItem, LOGFONT& logfont)
{
	Color_Font cf;
	if( !m_mapColorFont.Lookup( hItem, cf ) )
		cf.color = (COLORREF)-1;
	cf.logfont = logfont;
	m_mapColorFont[hItem] = cf;
}

//////////////////////////////////////////////////////////////////////
void GTreeCtrl::SetItemBold(HTREEITEM hItem, BOOL bBold)
{
	SetItemState(hItem, bBold ? TVIS_BOLD: 0, TVIS_BOLD);
}

//////////////////////////////////////////////////////////////////////
void GTreeCtrl::SetItemColor(HTREEITEM hItem, COLORREF color)
{
	Color_Font cf;
	if(!m_mapColorFont.Lookup(hItem, cf))
		cf.logfont.lfFaceName[0] = '\0';
	cf.color = color;
	m_mapColorFont[hItem] = cf;
}

//////////////////////////////////////////////////////////////////////
BOOL GTreeCtrl::GetItemFont(HTREEITEM hItem, LOGFONT * plogfont)
{
	Color_Font cf;
	if(!m_mapColorFont.Lookup(hItem, cf))
		return FALSE;
	if(cf.logfont.lfFaceName[0] == '\0') 
		return FALSE;
	*plogfont = cf.logfont;
	return TRUE;

}

//////////////////////////////////////////////////////////////////////
BOOL GTreeCtrl::GetItemBold(HTREEITEM hItem)
{
	return GetItemState(hItem, TVIS_BOLD) & TVIS_BOLD;
}

//////////////////////////////////////////////////////////////////////
COLORREF GTreeCtrl::GetItemColor(HTREEITEM hItem)
{
	// Returns (COLORREF)-1 if color was not set
	Color_Font cf;
	if(!m_mapColorFont.Lookup(hItem, cf))
		return (COLORREF) - 1;
	return cf.color;

}

//////////////////////////////////////////////////////////////////////
void GTreeCtrl::OnPaint() 
{
	CPaintDC dc(this);

	// Create a memory DC compatible with the paint DC
	CDC memDC;
	memDC.CreateCompatibleDC(&dc);

	CRect rcClip, rcClient;
	dc.GetClipBox( &rcClip );
	GetClientRect(&rcClient);

	// Select a compatible bitmap into the memory DC
	CBitmap bitmap;
	bitmap.CreateCompatibleBitmap( &dc, rcClient.Width(), rcClient.Height() );
	memDC.SelectObject( &bitmap );
	
	// Set clip region to be same as that in paint DC
	CRgn rgn;
	rgn.CreateRectRgnIndirect( &rcClip );
	memDC.SelectClipRgn(&rgn);
	rgn.DeleteObject();
	
	// First let the control do its default drawing.
	CWnd::DefWindowProc(WM_PAINT, (WPARAM)memDC.m_hDC, 0);

	HTREEITEM hItem = GetFirstVisibleItem();

	int iItemCount = GetVisibleCount() + 1;
	while(hItem && iItemCount--)
	{		
		CRect rect;

		// Do not meddle with selected items or drop highlighted items
		UINT selflag = TVIS_DROPHILITED | TVIS_SELECTED;
		Color_Font cf;
	
		//if ( !(GetTreeCtrl().GetItemState( hItem, selflag ) & selflag ) 
		//	&& m_mapColorFont.Lookup( hItem, cf ))
		
		if ((GetItemState(hItem, selflag) & selflag) 
			&& ::GetFocus() == m_hWnd)
			;
		else if (m_mapColorFont.Lookup(hItem, cf))
		{
			CFont *pFontDC;
			CFont fontDC;
			LOGFONT logfont;

			if(cf.logfont.lfFaceName[0] != '\0') 
				logfont = cf.logfont;
			else {
				// No font specified, so use window font
				CFont *pFont = GetFont();
				pFont->GetLogFont( &logfont );
			}

			if(GetItemBold(hItem))
				logfont.lfWeight = 700;

			fontDC.CreateFontIndirect(&logfont);
			pFontDC = memDC.SelectObject(&fontDC );

			if(cf.color != (COLORREF) - 1)
				memDC.SetTextColor(cf.color);
			else
				memDC.SetTextColor(GetSysColor(COLOR_WINDOWTEXT));


			CString sItem = GetItemText(hItem);

			GetItemRect(hItem, &rect, TRUE);
			memDC.SetBkColor( GetSysColor(COLOR_WINDOW));
			memDC.TextOut(rect.left + 2, rect.top + 1, sItem);
			
			memDC.SelectObject(pFontDC);
		}
		hItem = GetNextVisibleItem(hItem);
	}


	dc.BitBlt(rcClip.left, rcClip.top, rcClip.Width(), rcClip.Height(), &memDC, 
				rcClip.left, rcClip.top, SRCCOPY);

	memDC.DeleteDC();
}

/*
	no[nEntry+nCls*2+1]

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/11/2009
*/
void _Sort_On_X_( GIP_LIB *hLib,int *no,int flag )	{
	int nEntry=hLib->nEntry,nCls=hLib->nClass,i,j,t,*s_no,*s_ptr,nz,cur;
	GIP_LIB_entry *arrEntry=hLib->arrEntry,*h_1;
	double *w_s=(double*)GIP_calloc( sizeof(double),nCls ),a;
	
	for( i=0; i < nEntry; i++ )
		no[i] = i;

	s_no=no+nEntry;		s_ptr=s_no+nCls;
	for( i=0; i < nCls; i++ )	{
		s_no[i] = i;			s_ptr[i]=0;
	}

	for( i = 0; i < nEntry; i++ )	{
		h_1 = arrEntry+i;
		w_s[h_1->cls] += h_1->x;
		s_ptr[h_1->cls]++; 
	}	
	s_ptr[nCls] = nEntry;
	for( i=nCls-1; i >= 0; i-- )	{
		s_ptr[i] = s_ptr[i+1]-s_ptr[i];
	}
	ASSERT( s_ptr[0]==0 );

	for( i=0; i < nCls; i++ )	{
		for( j=i+1; j < nCls; j++ )	{
			if( w_s[i]<w_s[j] ){
				t=s_no[j];	s_no[j]=s_no[i];	s_no[i]=t;	
				a=w_s[i];	w_s[i]=w_s[j];		w_s[j]=a;
			}
		/*	h_1=arrEntry+no[i];		h_2=arrEntry+no[j];
			if( h_1->x<h_2->x ){
				t=no[j];	no[j]=no[i];	no[i]=t;	
			}*/
		}	
	};
	nz = 0;
	for( i = 0; i < nCls; i++  )	{
		cur = s_no[i];
		for( j=s_ptr[cur]; j < s_ptr[cur+1]; j++ )	{
			no[nz++] = j;
		}
	}
	ASSERT( nz==nEntry );

	GIP_free( w_s );
}

// GTreeCtrl message handlers
/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/18/2006
*/
void GTreeCtrl::OnUpdate( G_ENTRY_STATUS status )	{
	int nNode=m_hLib->nEntry,nCls=m_hLib->nClass,cls_old,*No,nLevel,nItem=0,len;
	CString info,sPath,sTitle;
	int i,pos;
	GIP_LIB_entry *hEntry,*arrEntry=m_hLib->arrEntry;

	DeleteAllItems( );
	if( nNode==0 )
		return;

	sPath = m_hLib->sPath;
	sTitle = GetFileTitle(sPath);
	sTitle = sTitle==_T("") ? _T("New_Library_") : sTitle;
	HTREEITEM hRoot = InsertItem( sTitle ),hItem,hChild;
	VERIFY( SetItemData( hRoot,-1 ) );

	cls_old = -1;
	hItem=GIP_NULL;
	No = (int*)GIP_alloc( sizeof(int)*(nNode+nCls*2+1) );
	if( m_isSort==1 )	{
		_Sort_On_X_( m_hLib,No,0x0 );
	}else	{
		for( i = 0; i < nNode; i++ )
			No[i] = i;
	}
	for( i = 0; i < nNode; i++ )	{
		hEntry = arrEntry+No[i];
		switch( status )	{
		case G_ENTRY_ALL:
			break;
		case G_ENTRY_FAILED:
			if( hEntry->status!=G_ENTRY_FAILED )
				continue;
			break;
		case G_ENTRY_WARNING:
			if( hEntry->status==G_ENTRY_OK )
				continue;
			break;
		default:
			ASSERT( FALSE );
			break;
		}
//		if( status!=G_ENTRY_ALL && hEntry->status!=status )
//			continue;
		nItem++;
		sPath = CString(hEntry->sFilePath);
		if( hEntry->cls!=cls_old )	{
			info = GetFileDir( sPath );
			pos = 0;		nLevel=0;
			while( (pos = info.Find( '\\',pos )) !=-1 )	{
				pos = pos+1;	nLevel++;
			}
			len=info.GetLength( );
			if( nLevel>3 || len>25 )	{		//只需最后3层目录
				pos = 0;		
				while( (pos = info.Find( '\\',pos )) !=-1 )	{
					pos = pos+1;	nLevel--;
					if( len-pos+1<25 && nLevel==3 )
						break;
					if( nLevel==2 )
						break;
				}
				info = info.Right( len-pos+1 );
			}

		//	info.Format( _T("cls_%d"),hEntry->cls+1 );
			hItem = InsertItem( info,hRoot );
			SetItemColor( hItem,RGB(0,0,0) );
			
			VERIFY( SetItemData( hItem,-(hEntry->cls+1) ) );
			cls_old = hEntry->cls;
		}
		info = GetFileTitle(sPath);
		hChild = InsertItem( info,hItem );
		if( hEntry->status==G_ENTRY_FAILED )	{			
			SetItemColor( hChild,RGB(255,0,0) );		
		}else if( hEntry->status==G_ENTRY_WARNING )	{
			SetItemColor( hChild,RGB(155,0,100) );		
		}
		VERIFY( SetItemData( hChild,No[i] ) );
//		SetItemText( hChild,GetFileTitle(sPath) );

	}
	GIP_free( No );
	
	CString sTitle_1;
	sTitle_1.Format( _T("%s (%d items)"),sTitle,nItem );
	SetItemText( hRoot,sTitle_1 );
	Expand( hRoot,TVE_EXPAND );
}

/*
	v0.1	cys
		3/19/2006
*/
void GTreeCtrl::OnSelchanged(NMHDR* pNMHDR, LRESULT* pResult) {
//	NM_TREEVIEW* pNMTreeView = (NM_TREEVIEW*)pNMHDR;
	// TODO: Add your control notification handler code here
	*pResult = 0;
	if( g_pFrame == NULL )	
		return;
	int mode=(int)(g_Lib.user_control[GIRIS_MODE_]);
	if( mode==GIRIS_MODE_AUTHEN )
		return;

	LPNMTREEVIEW pnmtv = (LPNMTREEVIEW) pNMHDR;
	HTREEITEM hItem = pnmtv->itemNew.hItem;
	if( hItem==NULL )
		return;
	int no=(int)(GetItemData( hItem ));
	if( no<0 )	{
		if( g_hTarget==GIP_NULL )
			g_hSource=GIP_NULL;
	}else	{
		ASSERT( no>=0 && no<m_hLib->nEntry );
		g_hSource = m_hLib->arrEntry+no;
	}
	if( no!=m_oldNo )	{		//修正一个奇怪的bug	4/27/2009
		m_oldNo = no;
		g_pFrame->Update( 0x0 );
	}
/*	ALG_TYPE alg=g_hSearch->type;
	CString info;
	int no=GetItemData( hItem ),i,j,child,size,cur,len;
	g_hSearch->curSel=no;
	HASH_TREE& trAB = g_hSearch->trAB;
	HT_NODE* hNode=trAB[no];
	if( no != trAB.curRoot )	{
		info.Format( _T("%s,%d,%d"),hNode->Tree_Dump( alg,true ),no,c_total[no] );
		VERIFY( SetItemText( hItem,info ) );
	}

	child = trAB.child_pos(no ),	size = trAB.child_size(no);
	cur = child;
	if( size > 0 && GetChildItem(hItem)==NULL )	{
		for( i = 0; i < size; i++,cur+=len )	{
			hNode = trAB[cur];
			len=ROUND_PLY_LEN(hNode->atak);
			for( j = 0; j < len; j++ )	{
				hNode = trAB[cur+j];			ASSERT( hNode->parent==no );
				info.Format( _T("%s,%d,%d"),hNode->Tree_Dump( alg,j==0 ),cur+j,c_total[cur+j] );
				hChild = InsertItem( info,hItem );
				VERIFY( SetItemData( hChild,cur+j ) );
			}
		}
	}

	CBoardControl *hBoard = &(g_pFrame->m_ctrlBoard);
	hBoard->Update( );*/
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/18/2009
*/
void GTreeCtrl::OnLButtonDblClk(UINT , CPoint ){
#ifdef 	_GIRIST_DEVELOP_
	if( g_pFrame == NULL || m_hLib==GIP_NULL )	
		return;
	HTREEITEM hItem =GetSelectedItem( );
	if( hItem==NULL )
		return;
	int no=(int)(GetItemData( hItem )),cls,i,nEntry=m_hLib->nEntry;
	if( no>=0 )	{
		ASSERT( no>=0 && no<m_hLib->nEntry );
		g_hSource = m_hLib->arrEntry+no;
		g_hSource->status = G_ENTRY_UNKNOWN;
		g_pFrame->Update( 0x0 );
		return;
	}

	GIP_LIB_entry *hEntry,*arrEntry=m_hLib->arrEntry;
	CString sPath,str;
	cls = -no-1;
	for( i = 0; i < nEntry; i++ )	{
		hEntry = m_hLib->arrEntry+i;		
		if( hEntry->cls!=cls )
			continue;
		sPath = CString(hEntry->sFilePath);
		break;
	}

	sPath = GetFileDir(sPath);
	str.Format( _T("Do you want to create \"%s\".lib and open?"),sPath );
	if( MessageBox( str,g_sEXE,MB_YESNO )!=IDYES )	
		return;

	if( 1 )	{
		g_pFrame->m_CreatePath = sPath;
		g_pFrame->OnLibCreate( );
	}else	{
		GIP_LIB *hNewLib=(GIP_LIB *)GIP_alloc( sizeof(GIP_LIB) );
		int *mark=(int*)GIP_alloc( sizeof(int)*m_hLib->nEntry );
		for( i = 0; i < m_hLib->nEntry; i++ )		mark[i]=0;
		while( hItem!=NULL )	{
			no=(int)(GetItemData( hItem ));
			if( no< 0 )		continue;
			ASSERT( no<m_hLib->nEntry );
			mark[no]=1;
		}
		GIP_LIB_extract( hNewLib,m_hLib,G_STR2TCHAR(sPath),m_hLib->hfv,mark,0x0 );	
		GIP_free( mark );
		GIP_free( hNewLib );

	//	g_pFrame->m_OpenPath = sPath;
	//	g_pFrame->OnLibCreate( );
	}


#endif
	return;

}
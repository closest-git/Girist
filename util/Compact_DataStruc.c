#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include "../PCH/GIP_def.h"
#include "gip_util.h"
#include "Compact_DataStruc.h"

/*
	v0.1	cys
		8/24/2008
*/
void PACT_ol_init( PACT_order_list*hPOL,int compact,int dim,int flag )	{
	int i,nHead;
//	ASSERT( compact < COMPACT_RANK_LIMIT );

	GIP_MEMCLEAR( hPOL,sizeof(PACT_order_list) );
	hPOL->dim=dim;
	hPOL->compact=compact;
	hPOL->nz=0;
	hPOL->base=0;
	nHead=compact;
	hPOL->head=(int*)GIP_alloc(sizeof(int)*nHead);
	hPOL->list=(int*)GIP_alloc(sizeof(int)*dim);
	hPOL->pre=(int*)GIP_alloc(sizeof(int)*dim);
	hPOL->next=(int*)GIP_alloc(sizeof(int)*dim);

	for( i = 0; i < dim; i++ )	{
		hPOL->list[i]=-1;	
		hPOL->pre[i]=-1;	
		hPOL->next[i]=-1;	
	}
	for( i = 0; i < nHead; i++ )	{
		hPOL->head[i]=-1;	
	}
}

/*
	v0.1	cys
		8/24/2008
*/
void PACT_ol_clear( PACT_order_list*hPOL )	{
	hPOL->dim=-1;
	hPOL->compact=-1;
	hPOL->nz++;
	if( hPOL->head!=GIP_NULL )
		GIP_free( hPOL->head );
	if( hPOL->list!=GIP_NULL )
		GIP_free( hPOL->list );
	if( hPOL->pre!=GIP_NULL )
		GIP_free( hPOL->pre );
	if( hPOL->next!=GIP_NULL )
		GIP_free( hPOL->next );
}

/*
	v0.1	cys
		8/24/2008
*/
int PACT_ol_maxpop( PACT_order_list*hPOL )	{
	int i,no=-1,*head=hPOL->head,*next=hPOL->next;

	for( i = hPOL->compact-1; i>=0; i-- )	{
		if( head[i]==-1 )	continue;
		no=head[i];
		head[i]=next[no];
		next[no]=-1;
		break;
	}
	if( no==-1 )		
		ASSERT( hPOL->nz==0 );
	else
		hPOL->nz--;

	return no;
}

/*
	v0.1	cys
		10/23/2008
*/
int PACT_ol_minpop( PACT_order_list*hPOL )	{
	int i,no=-1,*head=hPOL->head,*next=hPOL->next;

	for( i = 0; i < hPOL->compact; i++ )	{
		if( head[i]==-1 )	continue;
		no=head[i];
		head[i]=next[no];
		next[no]=-1;
		break;
	}
	if( no==-1 )		
		ASSERT( hPOL->nz==0 );
	else
		hPOL->nz--;

	return no;
}


/*
	v0.1	cys
		8/24/2008
*/
void PACT_ol_insert( PACT_order_list*hPOL,int no,int val )	{
	int cur,*head=hPOL->head,*next=hPOL->next,*pre=hPOL->pre;
	ASSERT( val>=0 && val<hPOL->compact && hPOL->list[no]==-1 );
	hPOL->list[no]=val;
	if( head[val]==-1 )	{
		head[val]=no;
	}else	{
		cur=head[val];
		next[no]=cur;		pre[cur]=no;
		head[val]=no;
	}
	hPOL->nz++;
	ASSERT( hPOL->nz<=hPOL->dim );
}

/*
	v0.1	cys
		8/24/2008
*/
void PACT_ol_update( PACT_order_list*hPOL,int no,int val )	{
	int cur,*head=hPOL->head,old=hPOL->list[no],*next=hPOL->next,*pre=hPOL->pre;
	ASSERT( val>=0 && val<hPOL->compact );
	ASSERT( old!=val && old!=-1 );
	hPOL->list[no]=val;

	if( pre[old]==-1 )	{		//remove from old bucket
		head[old]=next[old];
	}else	{
		next[pre[old]]=next[old];
	}

	if( head[val]==-1 )	{		//add to new bucket
		head[val]=no;
	}else	{
		cur=head[val];
		next[val]=cur;		pre[cur]=no;
		head[val]=no;
	}	
}

/*
	v0.1	cys
		1/14/2009
*/
int PACT_ol_getvector( PACT_order_list *hPOL,int val,int *vec )	{
	int cur,*head=hPOL->head,*next=hPOL->next,nz=0;

	cur =head[val];
	while( cur != -1 )	{
		vec[nz++] = cur;
		cur = next[cur];
	}
	return nz;
}

/*
	注意
		ptr,ind已分配内存

	v0.1	cys
		1/15/2009
*/
void PACT_ol_2ccs( PACT_order_list *hPOL,int *ptr,int *ind,int flag )	{
	int i,cur,*head=hPOL->head,*next=hPOL->next,nz=0;

	ptr[0]=0;
	for( i = 0; i < hPOL->compact; i++ )	{
		cur =head[i];
		while( cur != -1 )	{
			ind[nz++] = cur;
			cur = next[cur];
		}
		ptr[i+1]=nz;
	}
	ASSERT( nz==hPOL->nz );

	return ;
}
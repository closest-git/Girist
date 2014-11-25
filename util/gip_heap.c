#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include "../PCH/GIP_def.h"
#include "../util/gip_heap.h"
#include "../util/gip_util.h"

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		4/10/2008	
*/
void GIP_HEAP_init( GIP_HEAP *heap,int max_size )	{
	memset( heap,0x0,sizeof(GIP_HEAP) );

	heap->index=(int*)GIP_alloc(sizeof(int)*max_size );
	heap->val=(double*)GIP_alloc(sizeof(double)*max_size );
	heap->maxsize = max_size;
}

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		4/10/2008	
*/
void GIP_HEAP_clear( GIP_HEAP *heap )	{
	if( heap->index != GIP_NULL )
	{	GIP_free(heap->index);		heap->index=GIP_NULL;	} 
	if( heap->val != GIP_NULL )
	{	GIP_free(heap->val);		heap->val=GIP_NULL;	} 
	memset( heap,0x0,sizeof(GIP_HEAP) );
	heap->maxsize = -1;
}

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		4/1/2008	
*/
void GIP_HEAP_up( GIP_HEAP *heap,int base )	{
	double t_v,*hlist=heap->val;
	int parent,t_i,cur=base;

	parent = (cur-1)/2;
	t_v = hlist[cur];		t_i=heap->index[cur];
	while( cur != 0 )	{
		ASSERT( heap->index[parent]!=t_i );
		if( hlist[parent]<=t_v )
			break;
		else{
			hlist[cur] = hlist[parent];			heap->index[cur] = heap->index[parent] ;
			cur = parent;
			parent = (cur-1)/2;
		}
	}

	ASSERT( cur>=0 );
	hlist[cur] = t_v;		heap->index[cur] = t_i;
}

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		4/1/2008	
*/
void GIP_HEAP_down( GIP_HEAP *heap,int base )	{
	double t_v,*hlist=heap->val;
	int child,cur=base,size=heap->size,t_i;

	t_v = hlist[cur];			t_i=heap->index[cur];
	child = cur*2+1;
	while( child < size )	{
		if( child+1<size && hlist[child+1]<=hlist[child] )
			child=child+1;
		ASSERT( heap->index[child]!=t_i );
		if( t_v < hlist[child] )	
			break;
		else	{
			hlist[cur]=hlist[child];		heap->index[cur]=heap->index[child];
			cur = child;
			child = cur*2+1;
		}
	}
	ASSERT( cur<size );
	hlist[cur] = t_v;		heap->index[cur] = t_i;
}


/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		4/1/2008	
*/	
void GIP_HEAP_insert( GIP_HEAP *heap,int no,double val )	{
	int last=heap->size;
	ASSERT( last < heap->maxsize );
	heap->val[last] = val;
	heap->index[last] = no;
	GIP_HEAP_up( heap,last );
	heap->size++;
}

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		4/1/2008	
*/	
int GIP_HEAP_pop( GIP_HEAP *heap,int *no,double *val )		{
	int size=heap->size;
	
	if( size==0 )
		return GIP_HEAP_EMPTY;

	ASSERT( size > 0 );
	*no = heap->index[0];			*val = heap->val[0];

	heap->index[0]=heap->index[size-1];
	heap->val[0]=heap->val[size-1];
	heap->size--;
	if( heap->size > 0 )
		GIP_HEAP_down( heap,0 );

	return GIP_OK;
}

/*
	Copyright 2008-present, Grusoft

	v0.1	cys
		4/24/2008	
*/	
void GIP_HEAP_update( GIP_HEAP *heap,int no,double a )	{
	double t_v,*hlist=heap->val,old;
	int child,cur=0,size=heap->size,t_i,base=-1;

/*	t_v = hlist[cur];			t_i=heap->index[cur];
	while( t_i != no && cur < size )	{
		ASSERT( a >= t_v );
		child = cur*2+1;
		if( heap->index[child]==no )	{
			base=child;		break;
		}else if( child+1<size && heap->index[child+1]==no)	{
			base=child+1;		break;
		}else	{
			if( child+1<size && hlist[child+1]<=hlist[child] )
				child=child+1;			
		}
		cur = child;
		t_i = heap->index[cur];
	}*/
	for( cur = 0; cur < size; cur++ )	{
		if( heap->index[cur]==no )
		{	base=cur;		break;	}
	}
	ASSERT( base != -1 && heap->index[base]==no );
	old = hlist[base];		hlist[base]=a;
	if( a == old)	{
	}else if( a < old )	{
		GIP_HEAP_up( heap,base );
	}else	{
		GIP_HEAP_down( heap,base );
	}
}
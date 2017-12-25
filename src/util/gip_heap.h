#ifndef _GIP_CONTOUR_H_
#define _GIP_CONTOUR_H_

/*
	[50] Heap Sort
	A heap is a list with an ordering that stores a set of data in a complete binary tree.
	
*/

#define GIP_HEAP_EMPTY -1

typedef struct GIP_HEAP_tag GIP_HEAP;
struct GIP_HEAP_tag{
	int maxsize,size;			//行数，列数
	int *index;
	double *val;
};

#ifdef __cplusplus
extern "C"  {
#endif
	void GIP_HEAP_init( GIP_HEAP *heap,int max_size );
	void GIP_HEAP_clear( GIP_HEAP *heap );
	void GIP_HEAP_insert( GIP_HEAP *heap,int no,double val );
	int GIP_HEAP_pop( GIP_HEAP *heap,int *no,double *val  );

	void GIP_HEAP_update( GIP_HEAP *heap,int no,double a );	

#ifdef __cplusplus
}
#endif

#endif
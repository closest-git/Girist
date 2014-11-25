#ifndef _GIP_COMPACT_DATA_STRUC_H_
#define _GIP_COMPACT_DATA_STRUC_H_

/*
	Compact Data Structure
	v0.1	cys
		8/24/2008
	v0.2	cys
		3/18/2009		增加base


	Compact	取值有限适用于稀疏矩阵等
		元素个数<=dim，元素的值<COMPACT_RANK_LIMIT<100，

*/
#ifdef __cplusplus
extern "C"  {
#endif

#define COMPACT_RANK_LIMIT 100

	typedef struct PACT_order_list_tag PACT_order_list;
	struct PACT_order_list_tag{
		int dim,compact,nz,base;
		int *head,*list,*pre,*next;
	};	

	void PACT_ol_init( PACT_order_list*hPOL,int compact,int dim,int flag );
	void PACT_ol_clear( PACT_order_list*hPOL );
	void PACT_ol_insert( PACT_order_list*hPOL,int no,int val );
	void PACT_ol_update( PACT_order_list*hPOL,int no,int val );
	int PACT_ol_maxpop( PACT_order_list*hPOL );
	int PACT_ol_minpop( PACT_order_list*hPOL );

	int PACT_ol_getvector( PACT_order_list *hPOL,int val,int *vec );
	void PACT_ol_2ccs( PACT_order_list *hPOL,int *ptr,int *ind,int flag );

#ifdef __cplusplus
}
#endif

#endif
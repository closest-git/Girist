#ifndef _GIP_MESH_H_
#define _GIP_MESH_H_

typedef enum {		//
	GIP_MESH_TRI=0x10,
}GIP_MESH_TYPE;

typedef struct GIP_MESH_tag GIP_MESH;
struct GIP_MESH_tag{
	int type;
	int M,N,nNode,nEle,n_max;
	int G_nz,G_dim,*G_ptr,*G_ind,*tri_list;
	double *n_val,*n_cord;
	GIP_FLOATING *u;
};


#ifdef __cplusplus
extern "C"  {
#endif
	void GIP_MESH_gen( int M,int N,int ldu,int n_max,GIP_FLOATING *u,GIP_MESH *hMesh,int flag );
	void GIP_MESH_clear( GIP_MESH *hMesh );
	
#ifdef __cplusplus
}
#endif

#endif
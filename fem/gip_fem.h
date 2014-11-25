#ifndef _GIP_FEM_H_
#define _GIP_FEM_H_

/*
	FINITE ELEMENT METHOD

	1 linear,non linear PDE
	2 steady state, transient
	3 尽量设计成BLACK_BOX,以独立于调用

*/

//BOUNDARY CONDITION KIND
#define BOUNDARY_CONDITION_NONE			0x0
#define BOUNDARY_CONDITION_130			130			//special bc from [23] P.130

#ifdef __cplusplus
extern "C"  {
#endif

typedef void GIP_FEM_ELE;

typedef enum {		//
	GIP_FEM_GALERKIN=0x0,
	_FEM_ERR_ENERGY=0x100,_FEM_ERR_RES=0x200,
	GIP_FEM_STEADY=0x1000,GIP_FEM_TRANSIENT=0x2000,GIP_FEM_NONLINEAR=0x4000,

	_FEM_NONLINEAR_NEWTON=0x10000,_FEM_NONLINEAR_FIXPOINT=0x20000,
	_FEM_U_SCALE=0x100000
}GIP_FEM_METHOD;

	typedef struct GIP_FEM_NODE_tag GIP_FEM_NODE;
	struct GIP_FEM_NODE_tag{
		GIP_FLOATING r,c;
	};

	typedef struct GIP_FEM_TRI_tag GIP_FEM_TRI;
	struct GIP_FEM_TRI_tag{
		int nodes[3],g_map[9],flag[3];
		GIP_FLOATING area,N[9],K[9],C[9],f[3];
	};

	typedef struct GIP_FEM_tag GIP_FEM_;
	struct GIP_FEM_tag{
		GIP_MODEL model;
		GIP_FEM_METHOD method;
		int nNode,nEle,K_dim,bc_kind,loop;
		GIP_FEM_NODE *arrNode;			
		GIP_FEM_ELE *arrEle;
		int K_nz,*K_ptr,*K_ind,*I_temp;
		void *hSS;
		double *u,*f,xita,*bc_val,*D_temp,*K_val,dU,dT;
		double u_scale;
		GIP_BOUNDARY_CONDITION *bc_type;
		GIP_OBJECT *hObj;
	};

//	void GIP_FEM_preprocess( GIP_FEM_*hFEM,int M,int N,int ldu,GIP_FEM_METHOD method,GIP_FLOATING *u,int model );
	void GIP_FEM_preprocess( GIP_FEM_*hFEM,GIP_OBJECT *hObj,int model,GIP_FEM_METHOD method );
	void GIP_FEM_dump( GIP_FEM_*hFEM,int flag );
	void GIP_FEM_clear( GIP_FEM_ *hFEM );
	void GIP_FEM_core( GIP_FEM_*hFEM,int,int );
	void GIP_FEM_transient( GIP_FEM_*hFEM,double dT );

#ifdef __cplusplus
}
#endif

#endif
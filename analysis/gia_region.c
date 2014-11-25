/*
ISSUE-NEEDED:
	1 Region Growing:	A new approach
*/

#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/gip_util.h"
#include "../util/gip_ss.h"
#include "GIa_region.h"
#include "GIa_handle.h"

/*
	зЂвт:	
		1 only 4-connected

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/4/2006	
*/	
int _Region_Grow_0( GIa_REGION *hRegion,GIP_OBJECT *hObj,int seed,int *map,int noR,double thrsh,int*temp,int flag )	{
	int i,r,c,nz=0,top,*stack=temp,ldu=hObj->N,pos,pos_1;
	int neibor[]={-1,1,-ldu,ldu};
	GIP_FLOATING *data=hObj->data;
	double d_0=data[seed];

	ASSERT( map[seed]==-1 );

	r=GIP_POS2R( seed,ldu );		c=GIP_POS2C( seed,ldu );
	top=0;
	stack[top]=seed;
	while( top>=0 )	{
		pos = stack[top--];			nz++;		
		for( i = GIP_NEIBOR_LEFT; i <= GIP_NEIBOR_UP; i++ )	{
			pos_1 = pos + neibor[i];
			if( map[pos_1]!=-1 )	continue;
			if( fabs(data[pos_1]-d_0)>thrsh )	continue;
			stack[++top]=pos_1;
			map[pos_1]=noR;
			
		}
	}

	return nz;
}

/*
	get basic property of region

	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/4/2006	
*/	
void _Region_Basic( GIa_REGION *hRegion,GIP_OBJECT *hObj,int *map,int noR,int flag )	{
//	hRegion->area = nz;
	;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/4/2006	
*/	
void _Region_Boundary( GIa_REGION *hRegion,GIP_OBJECT *hObj,int *map,int noR,int flag )	{
	;
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/4/2006	
*/	
void GIa_HANDLE_region( GIa_HANDLE* hIa,GIP_OBJECT *hObj,int nFixSeed,int* seeds,int flag )	{
	int i,seed,nzMap,nR,nz,*map,*temp;
	int isSettingRNum,isSettingSeed;
	GIa_REGION *hRegion;
	double thrsh=hIa->rgn_thrsh;

	isSettingRNum = hIa->nRegion>0;
	isSettingSeed = nFixSeed>0;
	nz=hObj->M*hObj->N;
	map = GIP_alloc( sizeof(int)*nz*2 );		temp=map+nz;
	for( i = 0; i < nz; i++ )	map[i]=-1;
	nR = 0;				nzMap=0;
	while( nzMap < nz )	{
		if( isSettingRNum )
			hRegion = hIa->regions+nR;
		if( isSettingRNum && nR==hIa->nRegion-1 )	{
			seed = -1;
			for( i = 0; i < nz; i++ )	{
				if( map[i]!=-1 )	continue;
				map[i]=nR;			nzMap++;
			}
		}else	{
			if( isSettingSeed )	{
				for( i = 0; i < nFixSeed; i++ )	{
					if( map[seeds[i]]==-1 )	{
						seed = i;
						break;
					}
				}				}else	{
				for( i = 0; i < nz; i++ )	{
					if( map[i]==-1 )	{
						seed = i;
						break;
					}
				}	
			}
			nzMap += _Region_Grow_0( hRegion,hObj,seed,map,nR,thrsh,temp,0x0 );
		}
		_Region_Basic( hRegion,hObj,map,nR,0x0 );
		_Region_Boundary( hRegion,hObj,map,nR,0x0 );
		nR++;
	}

	GIP_free(map);
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/4/2006	
*/	
void GIa_Segment_thresh( GIa_HANDLE *hIa,GIP_OBJECT *hObj,int nTh,double *arrThresh,int flag )	{
	int M=hObj->M,N=hObj->N,MN=M*N,i,j;
	GIP_FLOATING *data=hObj->data;
	int isBorder;
	double off;
	
	for( i = 0; i < MN; i++ )	{
		isBorder = 0;
		for( j = 0; j < nTh; j++ )	{
			off=fabs(data[i]-arrThresh[j]);
			if( off<1.0 )
			{	isBorder=1;		break;	}
		}
		if( isBorder==1 )
			data[i]=0.0;
	}

//	if( 1 )		
//		GIP_save_IplImage_d( hObj->M,hObj->N,".\\trace\\segment_thrsh.jpg",hObj->data,hObj->N,0x0 );
}

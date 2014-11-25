#include <math.h>
#include "../pch/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/GIP_util.h"
#include "../util/GIP_heap.h"
#include "GIP_inpainting.h"

#define GIP_INP_NEIBOR 4
#define NEIBOR_LEFT		0
#define NEIBOR_RIGHT	1
#define NEIBOR_UP		2
#define NEIBOR_DOWN		3

/*
	Selecting the largest possible solution to the quadratic equation
		(max(dxT,0))^2+(max(dyT,0))^2=1.0

	Ref:	
		1 An Image Inpainting Technique Based on the Fast Marching Method.
		2 cvinpaint.cpp	line 233

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/5/2008	
*/
INP_FLOAT _inpaint_solve( GIP_INPAINT *hINP,int pos_1,int pos_2 )	{
	int *flag=hINP->flag;
	INP_FLOAT sol=1.0e6,*T=hINP->T,T_1,T_2,r,s;

	T_1=hINP->T[pos_1];		T_2=hINP->T[pos_2];		
	if( flag[pos_1]==GIP_INP_KNOWN )	{
		if( flag[pos_2]==GIP_INP_KNOWN )	{
			r = sqrt(2.0-(T_1-T_2)*(T_1-T_2) );
			s = (T_1+T_2-r)/2.0;
			if( s>=T_1 && s>=T_2 )
				sol = s;
			else	{
				sol=s+r;
				ASSERT( sol >= T_1 && sol >= T_2 );
			}
		}
		else
			sol = 1+T_1;			
	}else if( flag[pos_2]==GIP_INP_KNOWN ){
		sol = 1+T_2;
	}else{
		sol=1.0+min(T_1,T_2);
	}

	return sol;
}

/*
	inpainting one point

	Ref:	
		1 An Image Inpainting Technique Based on the Fast Marching Method.
		2 cvinpaint.cpp	line 233

	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/7/2008	
*/
void _inpaint_core( GIP_INPAINT *hINP,int pos )	{
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		4/5/2008	
*/
void GIP_inpainting( GIP_INPAINT *hINP )	{
	int nzBand=0,pos,pos_1,i,*flag=hINP->flag;
	int INP_neibor[GIP_INP_NEIBOR];
	GIP_HEAP heap;
	INP_FLOAT T,a1,a2,a3,a4;

	INP_neibor[NEIBOR_LEFT]=-1;				INP_neibor[NEIBOR_RIGHT]=1;
	INP_neibor[NEIBOR_UP]=hINP->ldr;		INP_neibor[NEIBOR_RIGHT]=-hINP->ldr;
	while( nzBand > 0 )	{
		GIP_HEAP_pop( &heap,&pos,&T );
		flag[pos] = GIP_INP_KNOWN;
		for( i = 0; i < GIP_INP_NEIBOR; i++ )	{
			pos_1 = pos + INP_neibor[i];
			if( flag[pos_1]==GIP_INP_KNOWN )	continue;
			if( flag[pos_1]==GIP_INP_INSIDE )		{
				flag[pos_1]=GIP_INP_BAND;
				_inpaint_core( hINP,pos_1 );
			}
			a1 = _inpaint_solve( hINP,pos_1+INP_neibor[NEIBOR_LEFT],pos_1+INP_neibor[NEIBOR_DOWN] );
			a2 = _inpaint_solve( hINP,pos_1+INP_neibor[NEIBOR_RIGHT],pos_1+INP_neibor[NEIBOR_DOWN] );
			a3 = _inpaint_solve( hINP,pos_1+INP_neibor[NEIBOR_LEFT],pos_1+INP_neibor[NEIBOR_UP] );
			a4 = _inpaint_solve( hINP,pos_1+INP_neibor[NEIBOR_RIGHT],pos_1+INP_neibor[NEIBOR_UP] );
		//	T[pos_1] = a;
			GIP_HEAP_insert( &heap,pos_1,T );
		}
	}
}

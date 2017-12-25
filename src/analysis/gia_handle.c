#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/gip_util.h"
#include "../util/gip_ss.h"
#include "GIa_handle.h"

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/4/2006	
*/	
void GIa_HANDLE_init( GIa_HANDLE* hIa,GIP_OBJECT *hObj,int nz_X,int flag)	{
	int M=hObj->M,N=hObj->N,nz=M*N;
	GIa_HISTOGRAM *hHisto=&(hIa->histo);
	double attrib[GIA_ATTRIB_MAX],mean,devia;
	memset( hIa,0x0,sizeof(GIa_HANDLE) );
	
//	GIP_mean_deviation( M,N,N,hObj->data,&mean,&devia );
	hIa->I_buffer = GIP_alloc( sizeof(int)*MAX(nz*3,nz_X) );
	hIa->D_buffer = GIP_alloc( sizeof(double)*MAX(nz*2,nz_X) );
	GIa_Histo_init( hObj,hHisto,0x0,hIa->I_buffer,0x0 );
	GIa_Histo_attrib( hHisto,attrib,0x0 );
	hIa->mean = attrib[GIA_ATTRIB_MEAN];
	hIa->deviation = attrib[GIA_ATTRIB_DEVIA];
	hIa->contrast = attrib[GIA_ATTRIB_CONTRAST];
}

/*
	Copyright 2008-present, Grusoft.
	v0.1	cys
		1/4/2006	
*/	
void GIa_HANDLE_clear( GIa_HANDLE* hIa )	{
	if( hIa->I_buffer!=GIP_NULL )
	{	GIP_free(hIa->I_buffer);	hIa->I_buffer=GIP_NULL;	}
	if( hIa->D_buffer!=GIP_NULL )
	{	GIP_free(hIa->D_buffer);	hIa->D_buffer=GIP_NULL;	}
	if( hIa->regions != GIP_NULL )	{
		GIP_free(hIa->regions);			hIa->regions=GIP_NULL;
	}

	GIa_Histo_clear( &(hIa->histo) );
}


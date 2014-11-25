/*
	level set methods for the PDEs
	1¡¢time-depenent 
*/
#include "grus_def.h"
#include "pde_basic.h"
#include <math.h>

/*
	one step of level set methods

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/19/2008	
*/
void lsm_step( GRUS_FLOATING *u )	{

}

/*
	A basic code for the level set methods

	Copyright 2008-present, Grusoft.
	v0.1	cys
		3/19/2008	
*/
void lsm_basic( )	{
	double t=0,t_final=1.0;
	GRUS_FLOATING *u;

	lsm_init( );
	while( t < t_final )	{
		lsm_step( u );
		if( 0 )
			;
		u = u_old;
	}
}
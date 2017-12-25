#ifndef _GIP_IMAG_PROPERTY_H_
#define _GIP_IMAG_PROPERTY_H_



#ifdef __cplusplus
extern "C"  {
#endif
	void GIP_centre( int M,int N,int ldu,GIP_FLOATING *u,double *temp,double *c_r,double *c_c,int flag );
	void GIP_mean_deviation_abs( int M,int N,int ldu,GIP_FLOATING *u,double *mean,double *devia );
//	void GIP_mean_deviation( int r_0,int r_1,int c_0,int c_1,int ldu,GIP_FLOATING *u,double *mean,double *devia );
	void GIP_mean_deviation( int M,int N,int ldu,GIP_FLOATING *u,double *mean,double *devia );		//ys

#ifdef __cplusplus
}
#endif

#endif
#ifndef _GIP_FV_H_
#define _GIP_FV_H_

/*
	Feature Vector for GIPACK

	v0.1	cys
		2/7/2009
*/

/*
	G_FV_MASEK:	MASEK's 1-D gabor filer
		[ref]: Recognition of Human Iris Patterns for Biometric Identification
	G_FV_CYS:	1-D gabor filer,left half ROI,
		[ref]: Recognition of Human Iris Patterns for Biometric Identification
*/
typedef enum{
	G_FV_DCT_CODE,
	G_FV_ROI=0x100,G_FV_COEFFICIENT,G_FV_MALI,G_FV_MASEK,
	G_FV_CYS,G_FV_CYS_1,G_FV_CYS_HALF,
	G_FV_TYPE_END
}G_FV_TYPE;

typedef enum{
	G_DCRIMI_1
}G_DCRIMI_TYPE;
//feature vector


#ifdef __cplusplus
extern "C"  {
#endif
	typedef struct G_FEATURE_VECTOR_tag G_FVector_;
	struct G_FEATURE_VECTOR_tag{
		void *hBase;
		G_FV_TYPE type;
		int flag,V_len,sf_M,sf_N,M_step,N_step;
		int bpSF,nzSF,unit_size;		//bits per subfeature,number of subfeatures 
		int nShift;
		unsigned char *V,*mask;			//V中不可避免有noise，由mask标记
	};

	typedef struct G_DISCRIMI_tag G_DCRIMI_;
	struct G_DISCRIMI_tag{
		void *hBase;
		G_DCRIMI_TYPE type;
		double rFAR,rFRR,rEER,D_,sep,eer_sep;
		double f_ar_8[8],f_rr_8[8],hd_8[8];	//对应1.0e-7~1.0的8个刻度
		double CRR,time;
		double nz_a,nz_r,mean_a,mean_r,devia_a,devia_r;
	};

	double G_fv_distance( unsigned char *fv,unsigned char *mask,G_FVector_ *hfv,int pos_1,int pos_2,char *temp );
	void G_fv_shift_( G_FV_TYPE type,int nShift,int v_len,G_U8 *V1,G_U8 *V_shift,int flag );

	void G_FVector_init( G_FVector_ *hFV,int M,int N,int nShift,int type );
	void G_FVector_set_step( G_FVector_ *hFV,int M_step,int N_step );
	void G_FVector_clear( G_FVector_ *hFV );
	void G_FV_info( G_FVector_ *hfv,int flag );

	void G_hamming_dist_( int nEntry,double *res,G_U8 *fv,G_U8 *mask,G_FVector_ *hfv,int *I_buffer,int flag );	

	void G_DCRIMI_set( void *hBase,G_DCRIMI_ *hDcrimi,int D_span,float *D_inter,float *D_intra,double crr,int flag );
	void G_DCRIMI_update_sep( G_DCRIMI_ *hDcrimi,int D_span,float *D_inter,float *D_intra,int flag );
#ifdef __cplusplus
}
#endif

#endif
#include <malloc.h>
#include <memory.h>
#include <FLOAT.h>
#include "../PCH/GIP_def.h"
#include "../PCH/GIP_package.h"
#include "../util/gip_util.h"
#include "GIP_pca_lda.h"

/*
	Sample is N-dimensional ==>  featrue space is m-dimensional

	Copyright 2008-present, Grusoft.
	v0.1	cys
		12/24/2008	
*/
void GIP_CalcPCA( int nSample,int N,int m,GIP_FLOATING *sample,GIP_FLOATING *mean,GIP_FLOATING *P )	{
	int i,j,p_s,p_d,step;
	double *db;

	CvMat* pMean = cvCreateMat(1, m, CV_32FC1);
	CvMat* pEigVals = cvCreateMat(1, m, CV_32FC1);
	CvMat* pEigVecs = cvCreateMat(m, m, CV_32FC1);
	CvMat* pData=cvCreateMat(100, 2, CV_64FC1);		//

/*	GIP_VEC2CVMAT( );*/
	db=pData->data.db;
	step=pData->step/sizeof(double);
	for( i = 0; i < nSample; i++)	{
		for( j = 0; j < N; j++ )	{
			p_d=GIP_RC2POS( i,j,step );
			db[p_d] = sample[p_s++];
		}
	}

	cvCalcPCA(pData, pMean, pEigVals, pEigVecs, CV_PCA_DATA_AS_ROW );
	
//	GIP_CVMAT2VEC( );
//	GIP_CVMAT2VEC( );
}

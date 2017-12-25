#ifndef _GIP_DEF_H_
#define _GIP_DEF_H_

#include <tchar.h>

/***/
/****/

#define GIP_APP_NAME		"GRUS IMAGE PACK"
#define GIP_NULL			0x0
#define GIP_INT32_MAX		2147483647    /*MAX_USER_SPACE supported by Win32; maximum (signed) int value*/
#define GIP_MAX_LINE       260

#define PI 3.1415926535897932384626433832795
#define MAX_LINE         300
#define MAX_PATH_LEN     100

#define GIP_SOLVER_CONVERGE 0x10

typedef double GIP_FLOATING;
typedef unsigned char G_U8;		//unsigned 8-bit

typedef struct { 
    double re;
    double im;
}GIP_FLOAT_Z;

typedef enum	{		//basic_attributes
	GIA_ATTRIB_MEAN=0x1,GIA_ATTRIB_DEVIA=0x2,GIA_ATTRIB_CONTRAST=0x3,
	GIA_ATTRIB_MAX=0x8
}GIA_BASIC_ATTRIBUTE;

typedef enum {		//MODELs
	GIP_MODEL_=0x0,
	GIP_SEG_MUMSHAH_2=0x10,GIP_SEG_MUMSHAH_s=0x11,GIP_MUMSHAH_filter=0x12,GIP_MUMSHAH_1=0x13,
	GIP_MODEL_CGM=0x20,
	GIP_SEG_CV=0x21,/*GIP_DECOMP_Hs=0x22,*/
	GIP_MODEL_130=0x130
}GIP_MODEL;

/*
	控制是否包含某个模块
	v0.1	cys
		3/25/2009
*/
typedef enum {		//DISTANCE
	GIP_M_ALL=0xFFFF,
	GIP_M_HISTO=0x010,
	GIP_M_SMOOTH=0x020,
}GIP_MODULE;

typedef enum {		//DISTANCE
	GIP_DISTANCE_EUCLID=0x01,GIP_DISTANCE_HAMMING
}GIP_DISTANCE;

typedef enum {		//标记PATH的多种含义：{	.in,目录，}
	GIP_PATH_INFILE=0x10,GIP_PATH_DIRECTORY=0x20
}GIP_PATH_FLAG;

typedef enum{
	G_FILE_UNKNOWN=0,
	G_IMAGE_BEGIN=0x10,
	G_IMAGE_JPG,G_IMAGE_BMP,
	G_IMAGE_PGM,		//http://netpbm.sourceforge.net/doc/pgm.html
	G_IMAGE_END,
	G_FILE_DIRECTORY=0x100
}G_FILE_TYPE;

/*
	为了便于filter处理，edge map的背景为0，edge为1				1/23/2009
*/
typedef enum {
	GIP_IMAG_ORIGIN=0x0,GIP_IMAG_EDGE,
}GIP_IMAG_TYPE;

typedef enum {		//ALGORITHM TYPE for SIGNED DISTANCE FUNCTION
	GIP_SDF_INTERSECTION,GIP_SDF_PROJECT_1
}GIP_SDF_ALG;

typedef enum {		//BASIC ALGORITHM for image processing: segmentation; 
	GIP_ALG_PYRAMIDS,GIP_ALG_WATERSHED
}GIP_BASIC_ALG;

typedef enum {		//BOUNDARY_CONDITION,影响fem,ss等
	GIP_BC_NONE=0x0,	GIP_BC_FIX=0x1
}GIP_BOUNDARY_CONDITION;

/*
	FINITE DIFFERENCE METHOD
*/
typedef enum {		//
//0xFF00 difference scheme
	GIP_SCHEME_FLAG=0xFF00,GIP_SCHEME_FORWARD=0x100,GIP_SCHEME_BACKWARD=0x200,GIP_SCHEME_CENTRAL=0x400,
}GIP_FD_METHOD;
#define GIP_FD_SCHEME(fdm)	((fdm)&GIP_SCHEME_FLAG)

#define GIP_PI			3.1415926535897932384626433832795			//
#define GIP_ZERO_ROUNDOFF 1.0e-6									//

typedef enum {	
	GIP_NEIBOR_LEFT=0,GIP_NEIBOR_RIGHT,GIP_NEIBOR_DOWN,GIP_NEIBOR_UP
}GIP_NEIBOR;

typedef enum {	
	GIP_ADJACENT_4=0,GIP_ADJACENT_8,GIP_ADJACENT_M
}GIP_ADJACENT_TYPE;

/*
	物体之间的相互关系
*/
enum	{
	GIa_INNER_OBJ=0x100,
};


/*
	利用tag的0-15位		0x00FF
*/
#define GIP_BD_LEFT		0x10
#define GIP_BD_RIGHT	0x20
#define GIP_BD_DOWN		0x40
#define GIP_BD_UP		0x80
/*
	0xF000
*/
#define GIP_INTERFACE	0x1000

/*
typedef struct GIP_IMAGE_tag GIP_IMAGE;
struct GIP_IMAGE_tag{
	int M,N;			//行数，列数
	GIP_FLOATING *val;
};
*/
#undef MAX
#undef MIN
#define MAX(a, b)    ( (a)>(b)?(a):(b) )
#define MIN(a, b)    ( (a)<(b)?(a):(b) )
/*
	//min & max without jumps from OPENCV
	#define  CV_IMIN(a, b)  ((a) ^ (((a)^(b)) & (((a) < (b)) - 1)))
	#define  CV_IMAX(a, b)  ((a) ^ (((a)^(b)) & (((a) > (b)) - 1)))
*/
//位操作
#define BIT_FLAG_1				0x00000001
#define BIT_FLAG_H				0x80
#define BIT_SET( val,flag ) ((val) |= (flag))	
#define BIT_RESET( val,flag ) ((val) &= (~(flag)) ) 
#define BIT_TEST( val,flag ) (((val)&(flag))==(flag))
#define BIT_IS( val,flag ) (((val)&(flag))!=0)

//(r,c) <==> pos	行优先排序		5/30/2008
#define GIP_RC2POS( r,c,ld )	((r)*(ld)+(c))
#define GIP_POS2R( pos,ld )		((pos)/(ld))
#define GIP_POS2C( pos,ld )		((pos)%(ld))
#define GIP_RC_VALID( r,r_0,r_1 )	((r)>=(r_0) && (r)<=(r_1) )
#define G_ENTROPY_( p )		( (p)<=0.0 ? 0.0 : -(p)*log(p) )

//
#define GIP_INTENSITY( a )		(int)( (a)+0.5 )
#define G_RADIEN2DEGREE( a )	(int)( (a)<0.0 ? (6.283185307179+(a))*57.295779513082+0.5 : (a)*57.295779513082+0.5 )
#define G_DOUBLE2INT(a)			(int)( (a)+0.5 )
#define G_CLOCK2SEC(a)			( (a)*1.0/CLOCKS_PER_SEC )
#define GIP_SWAP(a,b,t)			{ (t)=(a);	(a)=(b);	(b)=(t); }

//for independant of system API
#define GIP_MEMCOPY memcpy
#define GIP_MEMCLEAR( _Dst,_Size ) memset( (_Dst),(0x0),(_Size) )
#define GIP_MEM_FREE( mem ) 	if( (mem)!=GIP_NULL )	{	GIP_free( (mem) );		(mem)=GIP_NULL;	}

#define GIP_FWRITE_1( a,fp )	{	if( fwrite( &(a),sizeof(a),1,(fp) )!= 1)		goto _FIO_EXIT_;	}
#define GIP_FREAD_1( a,fp )		{	if( fread( &(a),sizeof(a),1,(fp) )!= 1)			goto _FIO_EXIT_;	}

#define GIP_CONTROL_ITEM 40

/* 前10个control的统一定义*/
#define GIP_ARCHI	1							/* SERIAL,SMP,MMP,OMP */
#define GIP_OS		2							/*  */
#define GIP_DUMP	5							/* DUMP FLAG */
#define GIP_MAT_SIGMA	9						/* SIGMA参数，适用于平移(A-sigma*I)等 */

/*
	BASIC IMAGE LIBRARY
*/
#define GIP_OPENCV_LIB			0x10000
#define GIP_CXIMG_LIB			0x10000

#define GIP_INFO_ITEM 90
#define GIP_STATUS			0	

/********************统一错误代码***************************************/
#define GIP_OK (0)
#define GIP_ERROR_argument_missing (-1)
#define GIP_ERROR_invalid_index (-2)
#define GIP_ERROR_empty_col (-3)
#define GIP_ERROR_invalid_matrix (-4)
#define GIP_ERROR_file_IO (-5)
#define GIP_ERROR_not_square (-7)
#define GIP_ERROR_invalid_value (-8)
#define GIP_ERROR_out_of_memory (-10)
#define GIP_ERROR_out_of_userpace (-11)
#define GIP_ERROR_zero_memory (-12)
#define GIP_ERROR_invalid_lower_tri_part (-13)
#define GIP_ERROR_zero_pivot (-20)
#define GIP_ERROR_0P_unstable (-21)
#define GIP_ERROR_node_level (-22)
#define GIP_ERROR_elim_count (-23)
#define GIP_ERROR_nodal_semap (-24)
#define GIP_ERROR_permut_vector (-25)
#define GIP_ERROR_wait_all (-26)
#define GIP_ERROR_mapping_file (-27)
#define GIP_ERROR_out_of_core (-28)
#define GIP_ERROR_alloc_oc_buffer (-29)

#define GIP_ERR_FMMVAL (-30)		//Strange value from fmm_solve

#define GIP_DRIVER_INIT_FAIL (-100)
#define GIP_LIB_VERSION_ERR (-101)			//GIP_LIB的版本号错误
#define GIP_MEASURE_ERR (-102)				//不能满足标准(阈值)

#define GIP_ERROR_internal_error (-1000) 
#define GIP_USER_BREAK (-2000)

/* 
	global error code for APPLICATIONS
*/
#define GIRIST_ERROR_start (-10000)		//for GIRIST
#define GBST_ERROR_start (-20000)		//for GBOOST
#define GFACE_ERROR_start (-30000)		//for GFACE

/*
	设定用户(APPLICATION)的参数不超过32个

	v0.1	cys
		3/16/2009
*/
#define GIP_USER_CONTROL          32
/* 
	GIRIST control parameter
	具体取值参见giris_core.h
*/
#define GIRIS_MODE_					6			/* 运行模式 */
#define GIRIS_T_SEP					7			/*  HD <= thrshold; decides person is same. */
#define GIRIS_T_GRAD				11			/* thrshold of gradient for finding iris */
#define GIRIS_T_INTN_0				12			/* thrshold of intensity for finding iris  */
#define GIRIS_T_INTN_1				13			/* thrshold of intensity for finding iris  */
#define GIRIS_T_MTCH_S				14			/* thrshold of match ratio of sclera  */
#define GIRIS_T_MTCH_P				15			/* thrshold of match ratio of pupil  */
#define GIRIS_FV_SIFT				16			/*   */
#define GIRIS_LG_WAVELEN			17			/*  wave length for log gabor filter */
#define GIRIS_M_NORMAL				18			/*  rows of normalized area */
#define GIRIS_N_NORMAL				19			/*  columns of normalized area */
#define GIRIS_GRAD_BIAS				20			/*  gradient direction的最大误差 */
#define GIRIS_GRAD_BIAS				20			/*  gradient direction的最大误差 */
#define GIRIS_SPOT_AREA				21			/*  光斑面积的最大比值 */
#define GIRIS_T_R0					22			/*  thrshold of radius */
#define GIRIS_SMOOTH_Q				23			/*  smooth paramenter */

/* 
	GFACE control parameter
	具体取值参见gip_face.h
*/
#define GFACE_MODE_					6			/* 运行模式 */
#define GFACE_SAMPLE_GRID			8			/* 训练例子的采样间隔*/
#define GFACE_CART_NSPLITs			11			/*  NSPLIT of CART tree */
#define GFACE_BOOST_TYPE			12			/*  type of ADABOOST */
#define GFACE_FEATURE_SCALE 		13			/*  rescaling of haar-like features */
#define GFACE_BACK_SAMPLE 			14			/*  number of backgroud sample */
#define GFACE_SAMPLE_ALGORITHM		15
#define GFACE_FEATURE_POOL			16			/*  BASIC,CORE,ALL */
#define GFACE_FALSE_ALAM			17			/*  false alarm rate */
#define GFACE_HIT_RATE				18			

/*
	注意与GIP_IMAGE的区别
	GIP_IMAGE是图像的原始数据
	GIP_DATA记录了预处理后的数据，通常要缩小，插值，归并等操作。
		并含fourier等变换数据
		并记录后处理的结果
		一般以GIP_FLOATING格式存储。
		直接由求解器等引用。求解器等并不负责其分配，释放！
		标准2D存储，无需ldu值！
		暂只支持灰度图像
*/
typedef struct GIP_OBJECT_tag GIP_OBJECT;
struct GIP_OBJECT_tag{
	int M,N;
	GIP_FLOATING *data;
};
enum GIP_OBJECT_FLAG{
	GIP_OBJ_2ZERO=0x100
};

enum	{
	GIP_MATRIX_ROW=0x10,	GIP_MATRIX_COL=0x20
};
/*
	暂借用GIP_OBJECT的结构
	简单即高效
	12/11/2008
*/
typedef GIP_OBJECT GIP_MATRIX;
typedef GIP_MATRIX* GIP_PT_MAT;

#ifdef __cplusplus
extern "C" {
#endif

extern double GIP_CONTROL[GIP_CONTROL_ITEM];
extern double GIP_INFO[GIP_INFO_ITEM];

void GIP_Get_Control_Info( double **control,double **info,int flag );
void GIP_OBJ_clear( GIP_OBJECT *hObj );
void G_PRINTF( const TCHAR *formatstring,... );
/* Sets error status and performs some additonal actions (displaying message box,
   writing message to stderr, terminating application etc.)
   depending on the current error mode */


#ifdef __cplusplus
}
#endif

#ifdef _GIP_MFC_APP_
	
#else
	#define GIP_USE_ASSERT

	#ifdef GIP_USE_ASSERT
		#include <assert.h>
		#define ASSERT(expression) (assert (expression))
	#else
		#define ASSERT(expression)
	#endif
#endif



#define GIP_EXIT            goto GIP_exit
/*
  GIP_ERROR macro unconditionally raises error with passed code and message.
  After raising error, control will be transferred to the exit label.
*/
#define GIP_ERROR( head,Code, Msg )                                 \
{                                                                   \
     GIP_Error( (Code), head, Msg, __FILE__, __LINE__ );				\
     GIP_EXIT;                                                          \
}

#endif

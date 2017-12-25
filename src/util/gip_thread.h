#ifndef _GIP_THREAD_H_
#define _GIP_THREAD_H_

#include <windows.h>

/*
	封装了基于windows的多线程操作

	1 具有与主进程交互协调的功能
*/

/*
	GIP_USER_MESSAGE	与主界面通讯用
*/
#define GIP_USER_MESSAGE		(WM_USER+1234)			//User Defined Message
enum							//These would be used by Thread - Should be inSync with DlgProc Values !!
{
	GUM_DISABLECONTROLS = GIP_USER_MESSAGE,
	GUM_PROGRESS_THREADCOMPLETED,
	GUM_PROGRESS_BARUPDATE,
	GUM_PROGRESS_TEXTUPDATE,
	GUM_PROGRESS_SETRANGE,

	GUM_CANCEL_PROGRESSTHREAD
};

enum{
	GIP_DRIVER_START,GIP_DRIVER_FINISH
};

#ifdef __cplusplus
extern "C" {
#endif
	/*
		与PROGRESS_CTRL通讯

		v0.1	cys
			3/11/2009
	*/
	typedef struct _UserCtrlData GU_CTRL_DATA;
	struct _UserCtrlData	{
		HWND	hThreadWnd;				//The Dialog that Created the Thread !!
		BOOL	bAlive;					//Indicates the Thread State Alive/Dead
		BOOL	bTerminate;				//Thread Monitors this to Know if it has to Terminate itself
		LPVOID	pUserParam;			//Parameter that shoud be sent to the UserProc
//		LP_CUPDIALOG_USERPROC	m_lpUserProc;	//User Progress Procedure - Called by ProgressDialogBox From the ThreadProc
	};
	extern GU_CTRL_DATA *g_hUserCtrl;

	void GU_PROGRESS_clear( GU_CTRL_DATA *hUserCtrl );
	void GU_SetUserCtrl( GU_CTRL_DATA *hUserCtrl,LPCTSTR lpszProgressText,UINT_PTR dwProgressbarPos,int flag );


	int GIP_driver_init( unsigned ( __stdcall *ThreadProc )( void * ),void* pArguments );
	int GIP_driver_clear( );

	void GIP_driver_core( void* hSolver,void **u_0,void **u_1,int *code,int flag  );
	int GIP_driver_switch( void *u_0,void *u_1,int code,int noThread );

#ifdef __cplusplus
  }
#endif



#endif

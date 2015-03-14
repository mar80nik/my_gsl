// my_gsl.cpp : Defines the initialization routines for the DLL.
//

#include "stdafx.h"
#include "my_gslApp.h"
#include "my_gsl.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

struct ParabolaFuncParams: public BaseForMultiFitterFuncParams
{
	enum {ind_c, ind_b, ind_a, ind_max};

	static double func(const double &x, const double *a, const size_t &p) {return 0;};	

	static double df_da(const double &x, const double *a, const size_t &p, double *c) {return 0;};	
	static double df_db(const double &x, const double *a, const size_t &p, double *c) {return 0;};	
	static double df_dc(const double &x, const double *a, const size_t &p, double *c) {return 0;};	

	ParabolaFuncParams( const DoubleArray& y, const DoubleArray& sigma ) : 
		BaseForMultiFitterFuncParams(ind_max, y, sigma)
	{
		pFunction = ParabolaFuncParams::func;
		pDerivatives[ind_a] = df_da; pDerivatives[ind_b] = df_db; pDerivatives[ind_c] = df_dc;
	}
};

class ParabolaFitFunc: public MultiFitterTemplate<ParabolaFuncParams>
{
public:
	double GetTop(double &x);	
};

//
//TODO: If this DLL is dynamically linked against the MFC DLLs,
//		any functions exported from this DLL which call into
//		MFC must have the AFX_MANAGE_STATE macro added at the
//		very beginning of the function.
//
//		For example:
//
//		extern "C" BOOL PASCAL EXPORT ExportedFunction()
//		{
//			AFX_MANAGE_STATE(AfxGetStaticModuleState());
//			// normal function body here
//		}
//
//		It is very important that this macro appear in each
//		function, prior to any calls into MFC.  This means that
//		it must appear as the first statement within the 
//		function, even before any object variable declarations
//		as their constructors may generate calls into the MFC
//		DLL.
//
//		Please see MFC Technical Notes 33 and 58 for additional
//		details.
//

// Cmy_gslApp

BEGIN_MESSAGE_MAP(Cmy_gslApp, CWinApp)
END_MESSAGE_MAP()


// Cmy_gslApp construction

Cmy_gslApp::Cmy_gslApp()
{
	// TODO: add construction code here,
	// Place all significant initialization in InitInstance
}


// The one and only Cmy_gslApp object

Cmy_gslApp theApp;


// Cmy_gslApp initialization

BOOL Cmy_gslApp::InitInstance()
{
	CWinApp::InitInstance();

	return TRUE;
}

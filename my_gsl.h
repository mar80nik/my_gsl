// my_gsl.h : main header file for the my_gsl DLL
//

#pragma once

#ifndef __AFXWIN_H__
	#error "include 'stdafx.h' before including this file for PCH"
#endif


// Cmy_gslApp
// See my_gsl.cpp for the implementation of this class
//

class Cmy_gslApp : public CWinApp
{
public:
	Cmy_gslApp();

// Overrides
public:
	virtual BOOL InitInstance();

	DECLARE_MESSAGE_MAP()
};

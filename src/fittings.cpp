#include "stdafx.h"
#include "fittings.h"

//////////////////////////////////////////////////////////////////////////
double ParabolaFuncParams::func( const double &x, const double *a, const size_t &p )
{
	return a[ind_a]*x*x + a[ind_b]*x + a[ind_c];
}
double ParabolaFuncParams::df_da( const double &x, const double *a, const size_t &p, double *c) { return x*x; }
double ParabolaFuncParams::df_db( const double &x, const double *a, const size_t &p, double *c) { return x; }
double ParabolaFuncParams::df_dc( const double &x, const double *a, const size_t &p, double *c) { return 1; }
double ParabolaFitFunc::GetTop(double &x)
{
	x = (-a[ParabolaFuncParams::ind_b]/(2*a[ParabolaFuncParams::ind_a])); 
	return GetXrelY(x);
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//f = C + (A/(1 + exp(-2*k*(x - B))));
double KneeFuncParams::func( const double &x, const double *a, const size_t &p )
{
	return a[ind_C] + (a[ind_A]/( 1 + exp(-2*a[ind_k]*(x - a[ind_B])) ));	
}
double KneeFuncParams::df_dA(const double &x, const double *a, const size_t &p, double *c)
{
	double t0 = *c;
	return 1/(1 + t0);	
}
double KneeFuncParams::df_dB(const double &x, const double *a, const size_t &p, double *c)
{
	double t0 = *c;
	return -2*a[ind_k]*a[ind_A]*t0/((1 + t0)*(1 + t0));
}
double KneeFuncParams::df_dC(const double &x, const double *a, const size_t &p, double *c)
{
	return 1;
}
double KneeFuncParams::df_dk(const double &x, const double *a, const size_t &p, double *c)
{
	double t0 = *c;
	return 2*(x - a[ind_B])*a[ind_A]*t0/((1 + t0)*(1 + t0));
}
double KneeFitFunc::GetInflection( double &x, const double &level )
{
	x = a[KneeFuncParams::ind_B] + log(level/(1-level))/(2*a[KneeFuncParams::ind_k]);
	return GetXrelY(x);
}

double * KneeFuncParams::PrepareDerivBuf( const double &x, const double *a, const size_t &p )
{
	buf = exp(-2*a[ind_k]*(x - a[ind_B]));
	return &buf;
}
/////////////////////////////////////////////////////////////////
////////////   f = A + C*exp(-(x - x0)^2/b)		     ////////////
/////////////////////////////////////////////////////////////////
double GaussFuncParams::func( const double &x, const double *a, const size_t &p )
{
	return a[ind_A] + a[ind_C]*exp(-(x - a[ind_x0])*(x - a[ind_x0])/a[ind_b]);
}
double GaussFuncParams::df_dA( const double &x, const double *a, const size_t &p, double *c )  { return 1; }
double GaussFuncParams::df_dC( const double &x, const double *a, const size_t &p, double *c )  { return c[1];}
double GaussFuncParams::df_dx0( const double &x, const double *a, const size_t &p, double *c ) { return 2*a[ind_C]*c[0]*c[1]; }
double GaussFuncParams::df_db( const double &x, const double *a, const size_t &p, double *c )  { return a[ind_C]*c[0]*c[0]*c[1]; }
double * GaussFuncParams::PrepareDerivBuf( const double &x, const double *a, const size_t &p )
{
	buf[0] = (x - a[ind_x0])/a[ind_b]; buf[1] = exp(-a[ind_b]*buf[0]*buf[0]);
	return &buf[0];
}
double GaussFitFunc::GetWidth()
{
	return 2*sqrt(a[GaussFuncParams::ind_b]);
}
/////////////////////////////////////////////////////////////////
////////////   f = A*sin^2(W*x + x0) + C			     ////////////
/////////////////////////////////////////////////////////////////
double Sin2FuncParams::func( const double &x, const double *a, const size_t &p )
{
	double Sin = sin(a[ind_W]*x + a[ind_x0]);
	return a[ind_A]*Sin*Sin + a[ind_C];
}

double Sin2FuncParams::df_dA( const double &x, const double *a, const size_t &p, double *c )
{
	return c[0]*c[0];
}

double Sin2FuncParams::df_dW( const double &x, const double *a, const size_t &p, double *c )
{
	return a[ind_A]*2*c[0]*c[1]*x;
}

double Sin2FuncParams::df_dx0( const double &x, const double *a, const size_t &p, double *c )
{
	return a[ind_A]*2*c[0]*c[1];
}

double Sin2FuncParams::df_dC( const double &x, const double *a, const size_t &p, double *c )
{
	return 1;
}

double * Sin2FuncParams::PrepareDerivBuf( const double &x, const double *a, const size_t &p )
{
	buf[0] = sin(a[ind_W]*x + a[ind_x0]); buf[1] = cos(a[ind_W]*x + a[ind_x0]);
	return &buf[0];	
}

double Sin2FitFunc::GetShift()
{
	return (M_PI - a[Sin2FuncParams::ind_x0])/(2*a[Sin2FuncParams::ind_W]);
}
#pragma  once

#include "my_gsl.h"

/////////////////////////////////////////////////////////////////
////////////    f = a*x*x + b*x + c    //////////////////////////
/////////////////////////////////////////////////////////////////
struct ParabolaFuncParams: public BaseForMultiFitterFuncParams
{
	enum {ind_c, ind_b, ind_a, ind_max};

	static double func(const double &x, const double *a, const size_t &p);	

	static double df_da(const double &x, const double *a, const size_t &p, double *c);	
	static double df_db(const double &x, const double *a, const size_t &p, double *c);	
	static double df_dc(const double &x, const double *a, const size_t &p, double *c);	

	ParabolaFuncParams( const DoubleArray& x, const DoubleArray& y, const DoubleArray& sigma ) : 
		BaseForMultiFitterFuncParams(ind_max, x, y, sigma)
	{
		pFunction = ParabolaFuncParams::func;
		pDerivatives[ind_a] = df_da; pDerivatives[ind_b] = df_db; pDerivatives[ind_c] = df_dc;
	}
};

class ParabolaFitFunc: public MultiFitterTemplate<ParabolaFuncParams>
{
public:
	virtual void GetReport(ControledLogMessage &log)
	{
		MultiFitterTemplate<ParabolaFuncParams>::GetReport(log);
		double x, y; y = GetTop(x); 
		log.T.Format("xmin = %g ymin = %g", x, y); log << log.T;
	}
	double GetTop(double &x);	
};
/////////////////////////////////////////////////////////////////
////////////   f = C + (A/(1 + exp(-2*k*(x - B))))   ////////////
/////////////////////////////////////////////////////////////////
struct KneeFuncParams: public BaseForMultiFitterFuncParams
{
	enum {ind_A, ind_B, ind_C, ind_k, ind_max};

	double buf;

	static double func(const double &x, const double *a, const size_t &p);	

	static double df_dA(const double &x, const double *a, const size_t &p, double *c);	
	static double df_dB(const double &x, const double *a, const size_t &p, double *c);	
	static double df_dC(const double &x, const double *a, const size_t &p, double *c);	
	static double df_dk(const double &x, const double *a, const size_t &p, double *c);	

	KneeFuncParams( const DoubleArray& x, const DoubleArray& y, const DoubleArray& sigma ) : 
		BaseForMultiFitterFuncParams(ind_max, x, y, sigma)
	{
		pFunction = KneeFuncParams::func;
		pDerivatives[ind_A] = df_dA; pDerivatives[ind_B] = df_dB; 
		pDerivatives[ind_C] = df_dC; pDerivatives[ind_k] = df_dk;
	}
	virtual double * PrepareDerivBuf(const double &x, const double *a, const size_t &p);	
};//////////////////////////////////////////////////////////////////////////


class KneeFitFunc: public MultiFitterTemplate<KneeFuncParams>
{
public:
	virtual void GetReport(ControledLogMessage &log, double level)
	{
		MultiFitterTemplate<KneeFuncParams>::GetReport(log);
		double x, y; y = GetInflection(x, level);
		log.T.Format("xmin = %g ymin = %g", x, y); log << log.T;
	}
	double GetInflection(double &x, const double &level);	
};
/////////////////////////////////////////////////////////////////
////////////   f = A + C*exp(-(x - x0)^2/b)		     ////////////
/////////////////////////////////////////////////////////////////
struct GaussFuncParams: public BaseForMultiFitterFuncParams
{
	enum {ind_A, ind_C, ind_x0, ind_b, ind_max};

	double buf[2];

	static double func(const double &x, const double *a, const size_t &p);	

	static double df_dA(const double &x, const double *a, const size_t &p, double *c);	
	static double df_dC(const double &x, const double *a, const size_t &p, double *c);	
	static double df_dx0(const double &x, const double *a, const size_t &p, double *c);	
	static double df_db(const double &x, const double *a, const size_t &p, double *c);	

	GaussFuncParams( const DoubleArray& x, const DoubleArray& y, const DoubleArray& sigma ) : 
		BaseForMultiFitterFuncParams(ind_max, x, y, sigma)
	{
		pFunction = GaussFuncParams::func;
		pDerivatives[ind_A] = df_dA; pDerivatives[ind_x0] = df_dx0; 
		pDerivatives[ind_C] = df_dC; pDerivatives[ind_b] = df_db;
	}
	virtual double * PrepareDerivBuf(const double &x, const double *a, const size_t &p);	
};//////////////////////////////////////////////////////////////////////////


class GaussFitFunc: public MultiFitterTemplate<GaussFuncParams>
{
public:
	double GetWidth();	
};
/////////////////////////////////////////////////////////////////
////////////   f = A*sin^2(W*x + x0) + C			     ////////////
/////////////////////////////////////////////////////////////////
struct Sin2FuncParams: public BaseForMultiFitterFuncParams
{
	enum {ind_A, ind_W, ind_x0, ind_C, ind_max};

	double buf[2];

	static double func(const double &x, const double *a, const size_t &p);	

	static double df_dA(const double &x, const double *a, const size_t &p, double *c);	
	static double df_dW(const double &x, const double *a, const size_t &p, double *c);	
	static double df_dx0(const double &x, const double *a, const size_t &p, double *c);	
	static double df_dC(const double &x, const double *a, const size_t &p, double *c);	

	virtual double * PrepareDerivBuf(const double &x, const double *a, const size_t &p);	

	Sin2FuncParams(  const DoubleArray& x, const DoubleArray& y, const DoubleArray& sigma ) : 
		BaseForMultiFitterFuncParams(ind_max, x, y, sigma)
	{
		pFunction = Sin2FuncParams::func;
		pDerivatives[ind_A] = df_dA;	pDerivatives[ind_W] = df_dW; 
		pDerivatives[ind_x0] = df_dx0;	pDerivatives[ind_C] = df_dC; 
	}
};
//////////////////////////////////////////////////////////////////////////


class Sin2FitFunc: public MultiFitterTemplate<Sin2FuncParams>
{
public:	
	double GetShift();
};
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "type_array.h"
#include "MyTime.h"
#include "MessageInspector.h"
#include "gsl\gsl_math.h"
#include "gsl\gsl_multimin.h"
#include "gsl\gsl_errno.h"
#include "gsl\gsl_math.h"
#include "gsl\gsl_roots.h"
#include "gsl\gsl_complex.h"
#include "gsl\gsl_complex_math.h"
#include "gsl\gsl_vector.h"
#include "gsl\gsl_blas.h"
#include "gsl\gsl_multifit_nlin.h"
#include "gsl\gsl_fft_real.h"
#include "gsl\gsl_fft_halfcomplex.h"

#define DEGREE (M_PI/180.)

void Convert_gsl_vector_to_DoubleArray(const gsl_vector* vector, DoubleArray& arr);
gsl_vector* CreateGSLReplica(const DoubleArray& arr);

struct SolverErrors 
{
	double abs, rel; 
	SolverErrors() {CleanUp();} 
	SolverErrors(double _abs, double _rel = 0) {CleanUp(); abs = _abs; rel = _rel; } 
	void CleanUp() {abs = rel = 0.;}
};
//////////////////////////////////////////////////////////////////////////
struct ComplexImGSL {double Im; ComplexImGSL(double i=1.) {Im=i;}};
struct ComplexReGSL {double Re; ComplexReGSL(double r=0) {Re=r;}};

#define cJ ComplexImGSL()

class ComplexGSL
{
public:
	gsl_complex z;

	ComplexGSL(double Re=0, double Im=0) { z=gsl_complex_rect(Re,Im); }
	ComplexGSL(gsl_complex c) { z=c; }
	ComplexGSL(ComplexReGSL r) { z=gsl_complex_rect(r.Re,0); }
	ComplexGSL(ComplexImGSL i) { z=gsl_complex_rect(0,i.Im); }
	~ComplexGSL() {};
	ComplexGSL operator-(const ComplexGSL& c)	{ return ComplexGSL(gsl_complex_sub(z,c.z)); }
	ComplexGSL operator-(const double r)		{ return ComplexGSL(gsl_complex_sub_real(z,r)); }
	ComplexGSL operator+(const ComplexGSL& c)	{ return ComplexGSL(gsl_complex_add(z,c.z)); }
	ComplexGSL operator+(const double r)		{ return ComplexGSL(gsl_complex_add_real(z,r)); }
	ComplexGSL operator/(const ComplexGSL& c)	{ return ComplexGSL(gsl_complex_div(z,c.z)); }
	ComplexGSL operator*(const ComplexGSL& c)	{ return ComplexGSL(gsl_complex_mul(z,c.z)); }
	void operator*=(const ComplexGSL& c)		{ z = gsl_complex_mul(z,c.z); }
	ComplexGSL operator*(const double r)		{ return ComplexGSL(gsl_complex_mul_real(z,r)); }
	ComplexGSL operator*(const ComplexImGSL& i)	{ return ComplexGSL(gsl_complex_mul_imag(z,i.Im)); }
	void operator*=(const double r)				{ z=gsl_complex_mul_real(z,r); return;}
	void operator*=(const ComplexImGSL i)		{ z=gsl_complex_mul_imag(z,i.Im); return;}
	double abs2()							{ return gsl_complex_abs2(z); }

};
ComplexGSL sqrt(ComplexGSL& c);
ComplexGSL pow2(ComplexGSL& c);
ComplexGSL exp(ComplexGSL& c);

//////////////////////////////////////////////////////////////////////////
//	This class is not allowed for use.
//	It is necessary to define [F.function] and [F.params] to be passed in this function.

struct SolverData 
{
	SolverErrors err;
	struct Counters 
	{
		size_t func_call, iter; 
		Counters() {CleanUp();} 
		void CleanUp() {func_call = iter = 0;}
	} cntr;
	ms dt; int status; size_t max_iter;

	SolverData(size_t _max_iter = 100)
	{
		CleanUp(); max_iter = _max_iter;
	}
	virtual ~SolverData() 
	{ 
		CleanUp(); 
	}
	virtual void CleanUp() {status = GSL_FAILURE; err.CleanUp(); cntr.CleanUp();}
};

enum SolverRegime {SINGLE_ROOT, MULTI_ROOT};

class BaseForFuncParams
{
public:
	virtual void PrepareBuffers() {}
	virtual void DestroyBuffers() {}
	BaseForFuncParams() {CleanUp();}
	virtual ~BaseForFuncParams() {DestroyBuffers();}
	virtual void CleanUp() {}
};

struct BoundaryConditions
{
	double min, max;
	BoundaryConditions() {min = max = 0; }
	BoundaryConditions(double _min, double _max) {min = _min; max = _max; }
};

template <class FuncParams>
class Solver1dTemplate: public SolverData
{
//************************************************//	
protected:
	static double func(double x, void * data)
	{
		Solver1dTemplate<FuncParams> *solver = (Solver1dTemplate<FuncParams>*)data; FuncParams *params = solver->params;
		solver->cntr.func_call++;		
		return (params->*(params->funcCB))(x);
	};
	typedef CArray<BoundaryConditions> BoundaryConditionsArray;
//************************************************//
protected:
	gsl_root_fsolver *s;
	const gsl_root_fsolver_type *fsolver_type;
	gsl_function F; size_t iter;
	SolverRegime rgm; int subrgns_max;	
	BoundaryConditionsArray SubRgns;
	FuncParams* params;
public:	
	DoubleArray Roots;	
public:
	Solver1dTemplate(const SolverRegime _rgm, const int _max_iter=100):
		SolverData(_max_iter), rgm(_rgm)
	{ 
		fsolver_type=gsl_root_fsolver_brent; s = NULL; 		
		F.function = func; F.params = this; params = NULL;
		subrgns_max = 50;
	}
	virtual ~Solver1dTemplate()
	{
		CleanUp();
	}
	virtual int Run(FuncParams *_params, const BoundaryConditions &_X, const SolverErrors &Err)
	{
		MyTimer Timer1; Timer1.Start(); CleanUp(); 
		params = _params; ASSERT(params); err = Err; size_t iter;
		params->PrepareBuffers(); FindSubRgns(_X, SubRgns); 
		s = gsl_root_fsolver_alloc (fsolver_type);

		for(int i = 0; i < SubRgns.GetSize(); i++)
		{
			BoundaryConditions& X = SubRgns[i]; iter = 0;
			gsl_root_fsolver_set (s, &F, X.min, X.max);
			do
			{
				iter++;
				status = gsl_root_fsolver_iterate (s);
				X.min = gsl_root_fsolver_x_lower (s); X.max = gsl_root_fsolver_x_upper (s);
				status = gsl_root_test_interval (X.min, X.max, err.abs, err.rel);
			}
			while (status == GSL_CONTINUE && iter < max_iter);	

			if( status == GSL_SUCCESS) 
			{		
				Roots << gsl_root_fsolver_root(s);
			}
			cntr.iter += iter;
		}
		dt=Timer1.StopStart(); 
		return status;
	}	
	virtual void CleanUp()
	{
		if (s != NULL) { gsl_root_fsolver_free(s); s = NULL; }
		Roots.RemoveAll(); SubRgns.RemoveAll(); 
		if (params != NULL)
		{
			params->CleanUp(); params = NULL;	
		}		
	}
public:
	virtual int FindSubRgns(const BoundaryConditions &X) {return FindSubRgns(X, SubRgns);}
	virtual int FindSubRgns (const BoundaryConditions &X, BoundaryConditionsArray& SubRgns)
	{
		if (rgm == SINGLE_ROOT)
		{
			SubRgns.Add(X);
		}
		else
		{
			double x, dx = (X.max - X.min)/(subrgns_max - 1); 
			BoundaryConditions y(func(X.min, F.params), 0);

			for(int i = 1; i < subrgns_max; i++)
			{
				x = X.max - i*dx; y.max = func(x, F.params);
				if ((y.min < 0 && y.max > 0) || (y.min > 0 && y.max < 0))
				{
					SubRgns.Add(BoundaryConditions(x, x + dx));
				}
				y.min = y.max;
			}
		}
		return GSL_SUCCESS;
	}
};

//////////////////////////////////////////////////////////////////////////
#define MAX_DELTA 1000
//////////////////////////////////////////////////////////////////////////
template <class FuncParams>
class MultiDimMinimizerTemplate: public SolverData
{
//************************************************//
public:
	static double func(const gsl_vector * x, void * data)
	{
		MultiDimMinimizerTemplate<FuncParams>* solver = (MultiDimMinimizerTemplate<FuncParams>*)data;
		solver->cntr.func_call++;
		return solver->params->func(x);
	};
//************************************************//
private:
	gsl_vector *X, *dX;
	gsl_multimin_fminimizer *s;	
	const gsl_multimin_fminimizer_type *fminimizer_type;	
	gsl_multimin_function F; 
	FuncParams* params;
public:
	DoubleArray Roots; double minimum_value;
//************************************************//
	MultiDimMinimizerTemplate(const int _max_iter=100): SolverData(_max_iter)
	{
		fminimizer_type = gsl_multimin_fminimizer_nmsimplex2; 
		s = NULL; X = dX = NULL; 
		F.f = func; F.params = this; params = NULL;
	}
	~MultiDimMinimizerTemplate()
	{
		CleanUp();
	}
	int Run(FuncParams* _params, DoubleArray &initX, DoubleArray &initdX, const SolverErrors &Err)
	{
		MyTimer Timer1; Timer1.Start(); double size; CleanUp(); 
		X = CreateGSLReplica(initX); dX = CreateGSLReplica(initdX); F.n = initX.GetSize(); 
		params = _params; ASSERT(params); params->PrepareBuffers(); err = Err;
		s = gsl_multimin_fminimizer_alloc (fminimizer_type, F.n);
		gsl_multimin_fminimizer_set (s, &F, X, dX); 
		do
		{
			cntr.iter++;
			status = gsl_multimin_fminimizer_iterate(s);
			if (status) break;
			size = gsl_multimin_fminimizer_size (s);
			status = gsl_multimin_test_size (size, err.abs);
		}
		while (status == GSL_CONTINUE && cntr.iter < max_iter);

		Convert_gsl_vector_to_DoubleArray(s->x, Roots); minimum_value = s->fval;
		params->DestroyBuffers();
		dt = Timer1.StopStart();
		return status;
	}
	virtual void CleanUp()
	{
		SolverData::CleanUp();
		if (s != NULL)		{gsl_multimin_fminimizer_free (s); s = NULL;}
		if (X != NULL)		{gsl_vector_free(X); X = NULL;}
		if (dX != NULL)	{gsl_vector_free(dX); dX = NULL;}
		Roots.RemoveAll(); 	
		if (params != NULL)
		{
			params->CleanUp(); params = NULL;	
		}
	}
	FuncParams * Attach_params(FuncParams *_params)
	{
		FuncParams *ret = params; params = _params; return ret;
	}
};
//////////////////////////////////////////////////////////////////////////
template <class FuncParams>
class Simple2DMinimizerTemplate: public SolverData
{
	//************************************************//
public:
	static double func(const gsl_vector * x, void * data)
	{
		Simple2DMinimizerTemplate<FuncParams>* solver = (Simple2DMinimizerTemplate<FuncParams>*)data;
		solver->cntr.func_call++;
		return solver->params->func(x);
	};
	//************************************************//
private:
	gsl_vector *film;
protected:
	FuncParams* params;
	DoubleArray range_min, range_max, dd; CArray<size_t> iter;
	size_t dim_size;
public:
	DoubleArray Roots; double minimum_value;
	//************************************************//
	Simple2DMinimizerTemplate(const int _max_iter=100): SolverData(_max_iter)
	{
		params = NULL; film = NULL;
	}
	~Simple2DMinimizerTemplate()
	{
		CleanUp();
	}
	int Run(FuncParams* _params, const DoubleArray& _range_min, const DoubleArray& _range_max, const DoubleArray& _dd)
	{
		MyTimer Timer1; Timer1.Start(); CleanUp(); 
		params = _params; ASSERT(params); params->PrepareBuffers(); 
		range_min = _range_min; range_max = _range_max; dd = _dd;
		dim_size = dd.GetSize(); ASSERT(dim_size == 2);
		film = gsl_vector_alloc (dim_size); ASSERT(film);	

		for (size_t i = 0; i< dim_size; i++) 
		{
			iter.Add(1 + (int)((range_max[i] - range_min[i])/dd[i]));
		}
		for (size_t i = 0; i < iter[0]; i++)
		{
			gsl_vector_set (film, 0, range_max[0] - dd[0]*i); 
			for (size_t j = 0; j < iter[1]; j++)
			{
				gsl_vector_set (film, 1, range_max[1] - dd[1]*j);

				double diff = func(film, this);
				if (diff == MAX_DELTA)
				{
					break;
				}
				else 
				{
					cntr.iter++;
					if (diff < minimum_value)
					{
						minimum_value = diff;
						Convert_gsl_vector_to_DoubleArray(film, Roots);					
					}
				}
			}
		}

		params->DestroyBuffers();
		dt = Timer1.StopStart();
		status = GSL_SUCCESS;
		return status;
	}	
	virtual void CleanUp()
	{
		SolverData::CleanUp();
		Roots.RemoveAll(); range_max.RemoveAll(); range_min.RemoveAll(); dd.RemoveAll(); iter.RemoveAll(); 
		minimum_value = MAX_DELTA;
		if (params != NULL) {params->CleanUp(); params = NULL;}
		if (film != NULL) {gsl_vector_free(film); film = NULL;}			
	}
	FuncParams * Attach_params(FuncParams *_params)
	{
		FuncParams *ret = params; params = _params; return ret;
	}
};
//////////////////////////////////////////////////////////////////////////
typedef double (*pFunc)(const double &x, const double *a, const size_t &p);
typedef double (*pDerivFunc)(const double &x, const double *a, const size_t &p, double *c);
struct BaseForMultiFitterFuncParams:public BaseForFuncParams
{
public:
	size_t n, p; const double *y; double *sigma;
	pDerivFunc *pDerivatives; pFunc pFunction;

	BaseForMultiFitterFuncParams(const size_t _p, const DoubleArray &_y, const DoubleArray &_sigma);
	virtual ~BaseForMultiFitterFuncParams();
	virtual double * PrepareDerivBuf(const double &x, const double *a, const size_t &p) { return NULL; };
	size_t GetPointsNum() {return n;}
	int f(const gsl_vector * a, gsl_vector * f);
	int df(const gsl_vector * a, gsl_matrix * J);
	int FillSigma(const DoubleArray &sigma);
};

struct BaseForFitFunc
{
	double leftmostX, rightmostX, dx;
	DoubleArray a, da; pFunc pFunction;

	BaseForFitFunc() { leftmostX = rightmostX = dx = 0; pFunction = NULL;}
	double GetXabsY(const double &x);
	double GetXrelY(double &x);
	HRESULT MakeGraph(DoubleArray &x, DoubleArray &y);
};

template <class FuncParams>
class MultiFitterTemplate: public SolverData, public BaseForFitFunc
{
private:
	const gsl_multifit_fdfsolver_type *multifit_fdfsolver_type; 
	gsl_multifit_fdfsolver *s; 
	size_t n, p;	
	gsl_multifit_function_fdf F;
	gsl_vector* initX;
	FuncParams*	params;
protected:
	static int f(const gsl_vector * x, void *data, gsl_vector * f)
	{
		MultiFitterTemplate<FuncParams>* solver = (MultiFitterTemplate<FuncParams>*)data;
		solver->cntr.func_call++;
		return solver->params->f(x, f);
	}
	static int df(const gsl_vector * x, void *data, gsl_matrix * J)
	{
		MultiFitterTemplate<FuncParams>* solver = (MultiFitterTemplate<FuncParams>*)data;
		solver->cntr.func_call++;
		return solver->params->df(x, J);
	}
	static int fdf(const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J)	
	{
		MultiFitterTemplate<FuncParams>* solver = (MultiFitterTemplate<FuncParams>*)data;
		solver->params->f(x, f);
		solver->params->df(x, J);
		return GSL_SUCCESS;
	}
public:	
	MultiFitterTemplate(const int _max_iter = 100): SolverData(_max_iter)
	{
		multifit_fdfsolver_type=gsl_multifit_fdfsolver_lmsder; s = NULL; initX = NULL;
		F.f = f; F.df = df; F.fdf = fdf; F.params = this; 
		params = NULL;		
	}
	virtual ~MultiFitterTemplate() {CleanUp();}
	virtual void CleanUp()
	{
		SolverData::CleanUp();
		if (s != NULL) { gsl_multifit_fdfsolver_free (s); s = NULL; }
		if (initX != NULL) { gsl_vector_free(initX); initX = NULL; }
		a.RemoveAll(); da.RemoveAll(); params = NULL;
	}
	int CalculateFrom(const DoubleArray& x, const DoubleArray& y, const DoubleArray& sigma, DoubleArray& init_a)
	{
		FuncParams params(y, sigma); 
		leftmostX = x[0]; rightmostX = x[params.n - 1]; dx = (rightmostX - leftmostX)/(params.n - 1);
		pFunction = params.pFunction;
		Run(&params, init_a, SolverErrors(1e-6));
		return status;
	}
	void GetReport(ControledLogMessage &log)
	{		
		log.T.Format("status = %s", gsl_strerror (status)); log << log.T;
		log.T.Format("----------------------------------"); log << log.T;
		if (status == GSL_SUCCESS)
		{
			for(int i = 0; i < a.GetSize(); i++)
			{
				log.T.Format("a%d = %g +/- %g%%", i, a[i], 100*da[i]/a[i]); log << log.T;
			}
		}
		else
		{
			log.SetPriority(lmprHIGH);
			for(int i = 0; i < a.GetSize(); i++)
			{
				log.T.Format("a%d = %g", i, a[i]); log << log.T;
			}
		}
		log.T.Format("time = %g ms", dt.val()); log << log.T;
		log.T.Format("func_cals = %d", cntr.func_call); log << log.T;
		log.T.Format("iter_num = %d", cntr.iter); log << log.T;
	}
protected:
	int Run(FuncParams* _params, DoubleArray& init_a, const SolverErrors &Err)
	{
		MyTimer Timer1; Timer1.Start(); CleanUp(); 
		if ((params = _params) != NULL)
		{
			if (init_a.GetSize() == FuncParams::ind_max)
			{
				params->PrepareBuffers(); 
				F.p = p = init_a.GetSize(); F.n = n = params->GetPointsNum(); 
				initX = CreateGSLReplica(init_a); err = Err;
				s = gsl_multifit_fdfsolver_alloc (multifit_fdfsolver_type, n, p);
				gsl_multifit_fdfsolver_set (s, &F, initX);
				do
				{
					cntr.iter++;
					status = gsl_multifit_fdfsolver_iterate (s);
					status = gsl_multifit_test_delta (s->dx, s->x, err.abs, err.rel);
				}
				while (status == GSL_CONTINUE && cntr.iter < max_iter);
				Convert_gsl_vector_to_DoubleArray(s->x, a);
				if (status == GSL_SUCCESS)
				{
					gsl_matrix *covar = gsl_matrix_alloc(p, p);
					if (covar != NULL)
					{
						gsl_multifit_covar(s->J, 0.0, covar);
						double c = GSL_MAX_DBL(1, gsl_blas_dnrm2(s->f) / sqrt((double)(n - p))); 
						Convert_gsl_vector_to_DoubleArray(&((gsl_matrix_diagonal (covar)).vector), da);
						for (int i = 0; i < da.GetSize(); i++)
						{
							da[i] = fabs(c*sqrt(da[i]));
						}
						gsl_matrix_free(covar);
					}
				}	
				params->DestroyBuffers(); 
			}
		}
		else
		{
				status = GSL_EINVAL;
		}
		dt=Timer1.StopStart(); params = NULL;
		return status;
	}
	HRESULT Fill_FitFunc(BaseForFitFunc* FitFunc)
	{
		if (FitFunc != NULL && params != NULL)
		{
			FitFunc->a = a; 
			FitFunc->leftmostX = params->leftmostX;
			FitFunc->rightmostX = params->rightmostX; 
			FitFunc->dx = params->dx; FitFunc->pFunction = params->pFunction;
			*((SolverData*)FitFunc) = *((SolverData*)this); 
			return S_OK;
		}
		return E_FAIL;		
	}
};
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
class FFTRealTransform: public SolverData
{
public:
	enum Direction {FORWARD, BACKWARD};
	struct Params: public SolverData
	{
		DoubleArray* y; Direction dir;
		Params() { y=NULL; dir=FORWARD;}
		Params(DoubleArray& _y) { y=&_y; }
	} params;
protected:
	gsl_fft_real_wavetable * real;
	gsl_fft_halfcomplex_wavetable * hc;
	gsl_fft_real_workspace * work;
	double *y; int n;

	int Init(Params& _params, Direction _dir);
	int Main();
	void CleanUp();
public:
	int status;
public:
	FFTRealTransform(int _max_iter=100) { real=NULL; hc=NULL; work=NULL; }
	virtual ~FFTRealTransform() {CleanUp();}
	int Run(Params& in,Direction dir);
};
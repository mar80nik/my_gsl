#pragma once

#include "MyTime.h"
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

template<class type>
class TypeArray: public CArray<type>
{
public:
	TypeArray()	{}
	TypeArray(const TypeArray& ref)	{ RemoveAll(); Copy(ref); }
	TypeArray& operator << (const type &d) {(*this).Add(d); return (*this);}
	TypeArray& operator =(const TypeArray& arr)
	{
		RemoveAll(); Copy(arr);
		return *this;
	}
	operator type*() {return CArray<type>::GetData();}
	BOOL operator==(const TypeArray &ref)
	{
		if (GetSize() != ref.GetSize()) return FALSE;
		for (int i = 0; i < GetSize(); i++)
		{
			if (operator[](i) != ref[i]) return FALSE;
		}
		return TRUE;
	}
	type get(int n) const
	{
		if (n < GetSize())
		{
			return operator[](n);
		}
		else
		{
			return 0;
		}
	}
	BOOL HasValues() const {return (GetSize() != 0);}
};

class DoubleArray: public TypeArray<double>
{
public:
	DoubleArray() {};
	DoubleArray(const DoubleArray& ref): TypeArray<double>(ref) {};

	gsl_vector* CreateGSLReplica();
	void operator= (const gsl_vector& vector);
	virtual void Serialize(CArchive& ar);
};

struct SolverErrors 
{
	double abs, rel; 
	SolverErrors() {CleanUp();} 
	SolverErrors(double _abs, double _rel = 0) {CleanUp(); abs = _abs; rel = _rel; } 
	void CleanUp() {abs = rel = 0.;}
};

struct PerfomanceInfoMk1 
{
	SolverErrors err;
	struct Counters {size_t func_call, iter; Counters() {CleanUp();} void CleanUp() {func_call = iter = 0;}} cntr;
	ms dt; int status; size_t max_iter;
	PerfomanceInfoMk1() {status = GSL_FAILURE; max_iter = 0; CleanUp();}
	virtual ~PerfomanceInfoMk1() { CleanUp(); }
	virtual void CleanUp() {status = GSL_FAILURE; max_iter = 0; err.CleanUp(); cntr.CleanUp();}
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
	ComplexGSL operator-(ComplexGSL& c)		{ return ComplexGSL(gsl_complex_sub(z,c.z)); }
	ComplexGSL operator-(double r)			{ return ComplexGSL(gsl_complex_sub_real(z,r)); }
	ComplexGSL operator+(ComplexGSL& c)		{ return ComplexGSL(gsl_complex_add(z,c.z)); }
	ComplexGSL operator+(double r)			{ return ComplexGSL(gsl_complex_add_real(z,r)); }
	ComplexGSL operator/(ComplexGSL& c)		{ return ComplexGSL(gsl_complex_div(z,c.z)); }
	ComplexGSL operator*(ComplexGSL& c)		{ return ComplexGSL(gsl_complex_mul(z,c.z)); }
	ComplexGSL operator*(double r)			{ return ComplexGSL(gsl_complex_mul_real(z,r)); }
	ComplexGSL operator*(ComplexImGSL& i)	{ return ComplexGSL(gsl_complex_mul_imag(z,i.Im)); }
	void operator*=(double r)				{ z=gsl_complex_mul_real(z,r); return;}
	void operator*=(ComplexImGSL i)			{ z=gsl_complex_mul_imag(z,i.Im); return;}
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
	double GetWidth() const { return (max - min);}
};
typedef CArray<BoundaryConditions> BoundaryConditionsArray;

template <class FuncParams>
class Solver1dTemplate: public SolverData
{
//************************************************//	
protected:
	static double func(double x, void * data)
	{
		Solver1dTemplate<FuncParams> *solver = (Solver1dTemplate<FuncParams>*)data; 
		FuncParams *params = solver->params;
		solver->cntr.func_call++;		
		return (params->*(params->funcCB))(x);
	};
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
			double y_l, y_r = func(X.max, F.params);

			for(int i = 1; i < subrgns_max; i++)
			{
				x = X.max - i*dx; y_l = func(x, F.params);
				if ((y_l < 0 && y_r > 0) || (y_l > 0 && y_r < 0))
				{
					SubRgns.Add(BoundaryConditions(x, x + dx));
				}
				y_r = y_l;
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
		X = initX.CreateGSLReplica(); dX = initdX.CreateGSLReplica(); F.n = initX.GetSize(); 
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

		Roots = *(s->x); minimum_value = s->fval;
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
						Roots = *film;					
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
	size_t n, p; const double *y; double *sigma; double leftmostX, rightmostX, dx;
	pDerivFunc *pDerivatives; pFunc pFunction;

	BaseForMultiFitterFuncParams(const size_t _p, const DoubleArray &_x, const DoubleArray &_y, const DoubleArray &_sigma);
	virtual ~BaseForMultiFitterFuncParams();
	virtual double * PrepareDerivBuf(const double &x, const double *a, const size_t &p) { return NULL; };
	size_t GetPointsNum() {return n;}
	int f(const gsl_vector * x, gsl_vector * f);
	int df(const gsl_vector * x, gsl_matrix * J);
	int FillSigma(const DoubleArray &sigma);
};


struct BaseForFitFunc: public SolverData
{
	double leftmostX, rightmostX, dx;
	DoubleArray a, da; pFunc pFunction;

	BaseForFitFunc(): SolverData() { leftmostX = rightmostX = dx = 0; pFunction = NULL;}
	double GetXabsY(const double &x);
	double GetXrelY(double &x);
};

template <class FuncParams>
class MultiFitterTemplate: public SolverData
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
	DoubleArray a, da;

public:
	MultiFitterTemplate(const int _max_iter = 100): SolverData(_max_iter)
	{
		multifit_fdfsolver_type=gsl_multifit_fdfsolver_lmsder; s = NULL; initX = NULL;
		F.f = f; F.df = df; F.fdf = fdf; 
		F.params = this; params = NULL;		
	}
	virtual ~MultiFitterTemplate() {CleanUp();}
	virtual void CleanUp()
	{
		SolverData::CleanUp();
		if (s != NULL) { gsl_multifit_fdfsolver_free (s); s = NULL; }
		if (initX != NULL) { gsl_vector_free(initX); initX = NULL; }
		a.RemoveAll(); da.RemoveAll(); params = NULL;
	}
	int Run(FuncParams* _params, DoubleArray& init_a, const SolverErrors &Err)
	{
		MyTimer Timer1; Timer1.Start(); CleanUp(); 
		params = _params; ASSERT(params); params->PrepareBuffers(); 
		F.p = p = init_a.GetSize(); F.n = n = params->GetPointsNum(); 
		initX = init_a.CreateGSLReplica(); err = Err;
		s = gsl_multifit_fdfsolver_alloc (multifit_fdfsolver_type, n, p);
		gsl_multifit_fdfsolver_set (s, &F, initX);
		do
		{
			cntr.iter++;
			status = gsl_multifit_fdfsolver_iterate (s);
			status = gsl_multifit_test_delta (s->dx, s->x, err.abs, err.rel);
		}
		while (status == GSL_CONTINUE && cntr.iter < max_iter);
		a = *(s->x);
		if (status == GSL_SUCCESS)
		{
			gsl_matrix *covar = gsl_matrix_alloc(p, p);
			if (covar != NULL)
			{
				gsl_multifit_covar(s->J, 0.0, covar);
				double c = GSL_MAX_DBL(1, gsl_blas_dnrm2(s->f) / sqrt((double)(n - p))); 
				da = (gsl_matrix_diagonal (covar)).vector;
				for (int i = 0; i < da.GetSize(); i++)
				{
					da[i] = fabs(c*sqrt(da[i]));
				}
				gsl_matrix_free(covar);
			}
		}	
		dt=Timer1.StopStart(); params->DestroyBuffers();
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

class FFTRealTransform: public PerfomanceInfoMk1
{
public:
	enum Direction {FORWARD, BACKWARD};
	struct Params: public PerfomanceInfoMk1
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
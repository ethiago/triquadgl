//------------------------------------------------
// Curvature.h
//------------------------------------------------
//
// Curvature types definition
// Version 0.1 - 19/04/2004
//
// Thomas Lewiner: thomas.lewiner@polytechnique.org
// Joao Gomes    : jgomes@mat.puc-rio.br
// Marcos Craizer: craizer@mat.puc-rio.br
// Helio Lopes   : lopes@mat.puc-rio.br
// Math Dept, PUC-Rio
//
//________________________________________________


#ifndef _CURVATURE_H_
#define _CURVATURE_H_

#ifdef WIN32
#pragma warning(disable:4786)
#endif // WIN32

#define DOUBLE_REAL
//#define FLOAT_REAL
//#define INTERVAL_REAL
//#define GMP_REAL


#ifdef DOUBLE_REAL
#undef FLOAT_REAL
#undef INTERVAL_REAL
#undef GMP_REAL

#include <math.h>

//________________________________________________
// real type : double
typedef double real ;
//------------------------------------------------
//________________________________________________

#endif // DOUBLE_REAL


//________________________________________________
//________________________________________________



#ifdef FLOAT_REAL
#undef DOUBLE_REAL
#undef INTERVAL_REAL
#undef GMP_REAL

#include <math.h>

//________________________________________________
// real type : float
typedef float real ;
//------------------------------------------------
//________________________________________________

#endif // FLOAT_REAL


//________________________________________________
//________________________________________________



#ifdef INTERVAL_REAL
#undef DOUBLE_REAL
#undef FLOAT_REAL
#undef GMP_REAL

#include "RInterval.h"

//________________________________________________
// real type : double interval
typedef RInterval real ;
//------------------------------------------------
//________________________________________________

#endif // INTERVAL_REAL


//________________________________________________
//________________________________________________


#ifdef GMP_REAL
#undef FLOAT_REAL
#undef INTERVAL_REAL
#undef DOUBLE_REAL

#include </home/tomlew/gmp-4.1.3/mpfr++.h>

//________________________________________________
// real type : mpf_class
typedef mpfr_class real ;
//------------------------------------------------
inline const real fabs ( const real &x ) { return abs(x) ; }
inline const real pow  ( const real &x, const real   &y ) { real p(x) ; return p.pow(y) ; }
inline const real pow  ( const real &x, const unsigned long &y ) { real p(x) ; return p.pow(y) ; }
inline const real pow  ( const real &x, const long   &y ) { real p(x) ; return p.pow(y) ; }
inline const real pow  ( const real &x, const int    &y ) { real p(x) ; return p.pow((long)y) ; }
inline const real pow  ( const real &x, const double &y ) { real p(x) ; return p.pow((long)y) ; }
inline const real atan2( const real &x, const real &y ) { return atan( x/y ) ; }
inline const int  rint ( const real &x ) { return (int)x ; }
//________________________________________________

#endif // GMP_REAL


#include <float.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#ifndef INTERVAL_REAL

#ifndef WIN32
inline int isnan( real x ) { return _isnan(x) ; }
static const real NaN  = (real) (-sqrt(-1.0));
inline const real nan() { return NaN ; }
inline int rint ( real x ) { return (int)x ; }

#pragma warning(disable:4786)

#else  // WIN32
#define nan() nan("n-char-sequence")
#endif // WIN32

#else // INTERVAL_REAL

#include "DIeee.h"
#define isnan IsNan
//static const real NaN  = real(nan(""));
//inline const real nan() { return NaN ; }
inline int rint ( real x ) { return (int)((float)x) ; }

#endif // INTERVAL_REAL

#endif // _CURVATURE_H_

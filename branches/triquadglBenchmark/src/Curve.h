//------------------------------------------------
// Curve.h
//------------------------------------------------
//
// Curve representation and basic operations
// Version 0.1 - 19/04/2004
//
// Thomas Lewiner: thomas.lewiner@polytechnique.org
// Joao Gomes    : jgomes@mat.puc-rio.br
// Marcos Craizer: craizer@mat.puc-rio.br
// Helio Lopes   : lopes@mat.puc-rio.br
// Math Dept, PUC-Rio
//
//________________________________________________


#ifndef _CURVE_H_
#define _CURVE_H_

#ifndef WIN32
#pragma interface
#endif // WIN32

#include <vector>
#include <string>
using std::vector ;
using std::string ;

#include "Curvature.h"


//________________________________________________
//
class Curve
//------------------------------------------------
{
// Members
private:
  int          _n ;            // number of samples
  vector<real> _x , _y  ;      // sample points
  vector<real> _tx, _ty ;      // tangents   of sample points
  vector<real> _nx, _ny ;      // normals    of sample points
  vector<real> _k       ;      // curvatures of sample points
  vector<real> _R2      ;      // R2 of sample points

  real         _xc, _yc ;      // central point of the curve
  real         _s , _ds ;      // curve extension
  real         _mk      ;      // mean value of k

//--------------------------------------------------------------------------
// Constructors
public:
  Curve ()
    :  _n(0), _x(0), _y(0), _tx(0), _ty(0), _nx(0), _ny(0), _k(0), _R2(0), _xc(0), _yc(0), _s(0), _ds(0), _mk(0) {}
  Curve (const Curve &c)
    :  _n(c._n), _x(c._x), _y(c._y), _tx(c._tx), _ty(c._ty), _nx(c._nx), _ny(c._ny), _k(c._k), _R2(c._R2), _xc(c._xc), _yc(c._yc), _s(c._s), _ds(c._ds), _mk(c._mk) {}
  Curve (const int n) ;

  Curve& operator= (const Curve&) ;
  ~Curve() {}
  void clear()
  {
    _n=0; _x.clear(); _y.clear(); _tx.clear(); _ty.clear(); _nx.clear(); _ny.clear(); _k.clear(); _R2.clear();
  }


//--------------------------------------------------------------------------
// Constant accessors
public :
  const int    &n      () const { return _n       ; }

  const real &xc     () const { return _xc      ; }
  const real &yc     () const { return _yc      ; }
  const real &ds     () const { return _ds      ; }
  const real &s      () const { return _s       ; }
  const real &mk     () const { return _mk      ; }

  const vector<real> &x  () const { return _x  ; }
  const vector<real> &y  () const { return _y  ; }
  const vector<real> &tx () const { return _tx ; }
  const vector<real> &ty () const { return _ty ; }
  const vector<real> &nx () const { return _nx ; }
  const vector<real> &ny () const { return _ny ; }
  const vector<real> &k  () const { return _k  ; }
  const vector<real> &R2 () const { return _R2 ; }

  const real &x  (const int &p) const { return _x [p] ; }
  const real &y  (const int &p) const { return _y [p] ; }
  const real &tx (const int &p) const { return _tx[p] ; }
  const real &ty (const int &p) const { return _ty[p] ; }
  const real &nx (const int &p) const { return _nx[p] ; }
  const real &ny (const int &p) const { return _ny[p] ; }
  const real &k  (const int &p) const { return _k [p] ; }
  const real &R2 (const int &p) const { return _R2[p] ; }


//--------------------------------------------------------------------------
// Accessors
public :
  vector<real> &x  () { return _x  ; }
  vector<real> &y  () { return _y  ; }
  vector<real> &tx () { return _tx ; }
  vector<real> &ty () { return _ty ; }
  vector<real> &nx () { return _nx ; }
  vector<real> &ny () { return _ny ; }
  vector<real> &k  () { return _k  ; }

  real &x  (const int &p) { return _x [p] ; }
  real &y  (const int &p) { return _y [p] ; }
  real &tx (const int &p) { return _tx[p] ; }
  real &ty (const int &p) { return _ty[p] ; }
  real &nx (const int &p) { return _nx[p] ; }
  real &ny (const int &p) { return _ny[p] ; }
  real &k  (const int &p) { return _k [p] ; }
  real &R2 (const int &p) { return _R2[p] ; }

  bool set_n(const int &np) ;
  void reset(const int &beg = 0, const int &end = -1) ;

//--------------------------------------------------------------------------
// Calculii
public :
  real norm_t   (const int &p) { return tx(p) * tx(p) + ty(p) * ty(p) ; }
  real dot_t_n  (const int &p) { return tx(p) * nx(p) + ty(p) * ny(p) ; }
  real compute_s(const int &i, const bool &improved = false ) ;
  static real R2_Estimator (const int &q, const real &x1, const real &x2, const real *const Dx, const real *const Ds,
                              const real &alpha = 1.0, const real &sigma = 0.0, const int &pw = 0 ) ;

  void compute_nosso_curv       (const int &q,
                                 const real &alpha = 1.0, const real &sigma = 0.0, const int &pw = 0,
                                 const int &coupled = 0, const bool &rot = false, const bool &improved = false) ;
  void compute_pouget_curv      (const int &q) ; // Pouget based on parabola fitting
  void compute_coeurjollyI_curv (const int &q) ; // Coeurjolly based on angle
  void compute_coeurjollyII_curv(const int &q) ; // Coeurjolly based on area
  void compute_worringI_curv    (const int &q,   // Worring orientation based (chain code)
                                 const real &sigma = 1.0/(16*16)) ;
  void compute_worringII_curv   (const int &q,   // Worring orientation based (resampling)
                                 const real &sigma = 1.0/(16*16)) ;
  void compute_worringIII_curv  (const int &q,   // Worring based on line fitting
                                 const real &sigma = 1.0/(16*16)) ;
  void compute_worringIV_curv   (const int &q,   // Worring based on path
                                 const real &sigma = 1.0/(16*16)) ;
  void compute_utcke_curv       (const int &q) ; // Utcke based on circle fitting
  void compute_belyaev_curv     (const int &q) ; // Belyaev via derivatives
  void compute_gumhold_curv     (const int &q,   // Gumhold based on external angles
                                 const real &alpha = 1.0) ;
  void compute_estrozi_curv     (const int &q) ; // Estrozi based on FFT

  int  refine() ;

  void buildT();

protected :
  int  resample(const int &p, const int &q, real *&xr, real *&yr) ;
};
//________________________________________________



#endif // _CURVE_H_

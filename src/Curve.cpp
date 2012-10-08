//------------------------------------------------
// Curve.cpp
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


#ifndef WIN32
#pragma implementation
#endif // WIN32

#include <time.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_randist.h>
#include "Curvature.h"
#include "Curve.h"


//___________________________________________________________________________________
//
Curve::Curve(const int n)
  :  _n(0), _x(0), _y(0), _tx(0), _ty(0), _nx(0), _ny(0), _k(0), _xc(0), _yc(0), _s(0), _ds(0), _mk(0)
//-----------------------------------------------------------------------------------
{
  set_n(n) ;
}
//___________________________________________________________________________________


//___________________________________________________________________________________
//
Curve& Curve::operator=(const Curve& c)
//-----------------------------------------------------------------------------------
{
  if( &c == this ) return *this ;
  _n       = c._n       ;
  _x       = c._x       ;
  _y       = c._y       ;
  _tx      = c._tx      ;
  _ty      = c._ty      ;
  _nx      = c._nx      ;
  _ny      = c._ny      ;
  _k       = c._k       ;
  _R2      = c._R2      ;
  _xc      = c._xc      ;
  _yc      = c._yc      ;
  _s       = c._s       ;
  _ds      = c._ds      ;
  _mk      = c._mk      ;

  return *this ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
//
bool Curve::set_n(const int &np)
//-----------------------------------------------------------------------------------
{
  if( np == _n ) return false ;
  _n = np ;
  _x .resize(_n,0) ;
  _y .resize(_n,0) ;
  _tx.resize(_n,0) ;
  _ty.resize(_n,0) ;
  _nx.resize(_n,0) ;
  _ny.resize(_n,0) ;
  _k .resize(_n,0) ;
  _R2.resize(_n,0) ;
  return true ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
//
void Curve::reset(const int &beg /*= 0*/, const int &end /*= -1*/ )
//-----------------------------------------------------------------------------------
{
  int _end = (end==-1) ? _n : end ;
  for( int p = beg ; p < _end ; ++p )
    _tx[p] = _ty[p] = _nx[p] = _ny[p] = _k[p] = _R2[p] = 0.0 ;
  _mk = 0 ;
}
//___________________________________________________________________________________




//___________________________________________________________________________________
//
inline real sign_of(const real &a) { return (a < 0.0) ? -1.0 : 1.0 ; }
//-----------------------------------------------------------------------------------


//___________________________________________________________________________________
//
#define NSUBPOINTS 10
#define SUBPOINTS_WITH_CURV 0
real Curve::compute_s(const int &i, const bool &improved /*= false*/)
//-----------------------------------------------------------------------------------
{
  real ex  = _x[i] - _x[i-1];
  real ey  = _y[i] - _y[i-1];
  real l   = hypot( ex,ey ) ;
  if( improved )
  {
    if( isnan( _tx[ i ] ) || isnan( _ty[ i ] ) || isnan( _nx[ i ] ) || isnan( _ny[ i ] ) ||
        isnan( _tx[i-1] ) || isnan( _ty[i-1] ) || isnan( _nx[i-1] ) || isnan( _ny[i-1] ) )
      return l ;

    ex /= l ;
    ey /= l ;

    real x1 = (_nx[i-1] * ex + _ny[i-1] * ey) / (_tx[i-1] * ex + _ty[i-1] * ey) ;
    real y1 = (_nx[ i ] * ex + _ny[ i ] * ey) / (_tx[ i ] * ex + _ty[ i ] * ey) ;

#ifdef SUBPOINTS_WITH_CURV
    real nrm ;
    nrm = hypot( (real)1.0,x1 ) ;
    real x2 = _k[i-1] * nrm * nrm * nrm ;
    nrm = hypot( (real)1.0,y1 ) ;
    real y2 = _k[ i ] * nrm * nrm * nrm ;
#endif // SUBPOINTS_WITH_CURV

    real p  = 0 ;
    real q  = 0 ;
    real t  = 0 ;
    real dt = 1.0/NSUBPOINTS ;
    real L  = 0 ;
    for( int j = 0 ; j < NSUBPOINTS ; ++j, t += dt )
    {
      p = q ;

#ifdef SUBPOINTS_WITH_CURV
      q = x1*t
        + 0.5*x2*t*t
        + (-(real)6*x1 - (real)4*y1 + (real)0.5*y2 - (real)1.5*x2)*t*t*t
        + ((real)8*x1 + (real)7*y1 - y2 + (real)1.5*x2)*t*t*t*t
        + (-(real)3*y1 + (real)0.5*y2 - (real)0.5*x2 - (real)3*x1)*t*t*t*t*t ;
#else  // SUBPOINTS_WITH_CURV
      q = x1*t + (-(real)2*x1 - y1)*t*t + (x1 + y1)*t*t*t ;
#endif // SUBPOINTS_WITH_CURV


      L += hypot( dt,p-q ) ;
    }
    return l * L ;
  }
  else
  {
    return l ;
  }
}
//___________________________________________________________________________________



//___________________________________________________________________________________
//
inline real weight (const real &d, const real &alpha /*= 1.0*/, const real &sigma /*= 0.0*/, const int &pw /*= 0*/ )
{
  return alpha * exp( -d*d * sigma ) / pow( d, pw ) ;
}
//-----------------------------------------------------------------------------------
real Curve::R2_Estimator (const int &q, const real &x1, const real &x2, const real *const Dx, const real *const Ds,
                            const real &alpha /*= 1.0*/, const real &sigma /*= 0.0*/, const int &pw /*= 0*/ )
//-----------------------------------------------------------------------------------
{
  real dxm = 0.0, dxem = 0.0, s1 = 0.0, s2 = 0.0, s3= 0.0, t1, t2 ;
  int i;

  real sw    = 0.0 ;
  for( i = 1; i <= 2*q; ++i )
  {
    real w = weight( Ds[i], alpha, sigma, pw ) ;
    t1 = Dx[i] * w ;
    t2 = ( Ds[i]*x1 + 0.5*Ds[i]*Ds[i]*x2 ) * w ;
    dxm  += t1 ;
    dxem += t2 ;
    s1   += t1 * t2 ;
    s2   += t1 * t1 ;
    s3   += t2 * t2 ;
    sw   += w  ;
  }
  dxm  /= sw ;
  dxem /= sw ;
  s1   -= dxm  * dxem ;
  s2   -= dxm  * dxm  ;
  s3   -= dxem * dxem ;

  return s1*s1/(s2*s3);
}
//___________________________________________________________________________________



//___________________________________________________________________________________
// functional: u = cos(alpha), v = sin(alpha)
bool minfun( const real m, const real n, const real u, const real v )
{
  return m*u + n*v < u*v ;
}
//-----------------------------------------------------------------------------------
// solves the functional: u = cos(alpha), v = sin(alpha)
void bissec( const real m, const real n, real &u, real &v )
{
  real un ; // negative approximation
  real up ; // positive approximation

  if( (double)(m*n) > 0 )
  {
    // negative approximation -> -mn/(m^2+n^2)
    un = -  n   / hypot( n ,m) ;
    if( (double)m > 0 )
      up = -(n+1.0) / hypot(n+1.0,m) ; // positive approximation ->  m*(sqrt(m^2+(n+1)^2)-(n+1))/(m^2+(n+1)^2)
    else
      up = -(n-1.0) / hypot(n-1.0,m) ; // positive approximation -> -m*(sqrt(m^2+(n-1)^2)-(n-1))/(m^2+(n-1)^2)
  }
  else
  {
    // positive approximation -> -mn/(m^2+n^2)
    up = -  n   / hypot( n ,m) ;
    if( (double)m > 0 )
      un = -(n-1.0) / hypot(n-1.0,m) ; // negative approximation -> -m*(sqrt(m^2+(n-1)^2)-(n-1))/(m^2+(n-1)^2)
    else
      un = -(n+1.0) / hypot(n+1.0,m) ; // negative approximation ->  m*(sqrt(m^2+(n+1)^2)-(n+1))/(m^2+(n+1)^2)
  }
/*
  v = sign_of(m) * sqrt( 1.0 - up*up ) ;
  if( !minfun(m,n,up,v) ) printf( "init+ u:%f, v:%f => %f\n", up,v, (m*up+n*v) - up*v ) ;
  v = sign_of(m) * sqrt( 1.0 - un*un ) ;
  if(  minfun(m,n,un,v) ) printf( "init- u:%f, v:%f => %f\n", un,v, (m*un+n*v) - un*v ) ;
*/
  while( true )
  {
    u  = (up + un) / 2.0 ;
    // applies the signal rules: u*n < 0, v*m > 0
    v = sign_of(m) * sqrt( 1.0 - u*u ) ;

//    if( u*n > 0 ) printf( "sign error: u:%f, n:%f\n", u, n ) ;
//    printf( "u:%f, v:%f => %f\n", u,v, (m*u+n*v) - u*v ) ;

    if( fabs((double)up - (double)un) < 1e-15 ) break ;

    if( minfun( m,n, u,v ) ) up = u ;
    else                     un = u ;
  }

//  printf( "%f ", (m*u+n*v) - u*v ) ;
}
//-----------------------------------------------------------------------------------
void Curve::compute_nosso_curv(const int &q,
                               const real &alpha /*= 1.0*/, const real &sigma /*= 0.0*/, const int &pw /*= 0*/,
                               const int &coupled /*= 0*/, const bool &rot /*= false*/, const bool &improved /*= false*/)
//-----------------------------------------------------------------------------------
{
  if( q < 1 ) return ;

  reset(  0   , q) ;
  reset(_n-q+1,_n) ;
  real *Ds = new real [2*q+1];
  real *Dx = new real [2*q+1];
  real *Dy = new real [2*q+1];

  _mk = 0 ;
  clock_t begin = clock() ;
  for( int p = q ; p < _n-q ; ++p )
  {
    int i,j ;

    Ds[0] = 0.0;
    Dx[0] = _x[p-q] - _x[p];
    Dy[0] = _y[p-q] - _y[p];

    for( i=p-q+1, j=1; i<=p+q; ++i, ++j )
    {
      Dx[j] = _x[i] - _x[p];
      Dy[j] = _y[i] - _y[p];
      Ds[j] = Ds[j-1] + compute_s( i, improved ) ;
    }

    real rx,ry, Rx,Ry ;

    if( rot )
    {
      rx = _tx[p] ;
      ry = _ty[p] ;
      real nrm = hypot( rx,ry ) ;

      if( nrm != 0.0 )
      { rx /= nrm ; ry /= nrm ; }
      else
      { rx = 1.0 ; ry = 0.0 ; }

      for(i=0;i<=(2*q);++i)
      {
        Rx = Dx[i] * rx + Dy[i] * ry ;
        Ry = Dy[i] * rx - Dx[i] * ry ;
        Dx[i] = Rx ;
        Dy[i] = Ry ;
      }

      if( coupled == 2 )
      {
        Ds[0] = 0.0 ;
        for(i=1;i<=(2*q);++i )
        {
          Ds[i] = Ds[i-1] + fabs( Dx[i] - Dx[i-1] ) ;
        }
      }
    }

    real a = 0.0;
    real b = 0.0;
    real c = 0.0;
    real e = 0.0;
    real f = 0.0;
    real g = 0.0;
    real h = 0.0;
    real m = Ds[q] ;
    for(i=0;i<=(2*q);++i)
    {
      Ds[i] -= m ;
      if( i == q ) continue ;
      real w = weight( Ds[i], alpha, sigma, pw ) ;
      w *= w ;
      a += Ds[i] * Ds[i] * w ;
      b += Ds[i] * Ds[i] * Ds[i] * w ;
      c += Ds[i] * Ds[i] * Ds[i] * Ds[i] * w ;
      e += Ds[i] * Dx[i] * w ;
      f += Ds[i] * Dy[i] * w ;
      g += Ds[i] * Ds[i] * Dx[i] * w ;
      h += Ds[i] * Ds[i] * Dy[i] * w ;
    }
    c/=4.0; b/=2.0; g/=2.0; h/=2.0;

    real quot = a*c-b*b ;

    if( coupled == 1 )
    {
      _tx[p] = (c*e-b*g)/quot;
      _ty[p] = (c*f-b*h)/quot;
      _nx[p] = (a*g-b*e)/quot;
      _ny[p] = (a*h-b*f)/quot;
      _k [p] = (e*h-g*f)/quot;

      quot   = hypot(_tx[p],_ty[p]) ;
      _k [p]/= quot*quot*quot ;

      quot = sign_of(_k[p]) * hypot( _nx[p],_ny[p] ) ;
      _nx[p] /= quot ;
      _ny[p] /= quot ;
    }
    else if( coupled == 0 )
    {
      real x1,y1, x2,y2, r2 ;

      x1 = (c*e-b*g)/quot;
      y1 = (c*f-b*h)/quot;
      if (fabs(x1) < fabs(y1))
      {
        y1 = sign_of(y1)*sqrt(1.0-x1*x1);
        x2 = (a*g-b*e)/quot;
        y2 = -(x1*x2)/y1;
        r2 = R2_Estimator (q,x1,x2,Dx,Ds,alpha,sigma,pw);
      }
      else
      {
        x1 = sign_of(x1)*sqrt(1.0-y1*y1);
        y2 = (a*h-b*f)/quot;
        x2 = -(y1*y2)/x1;
        r2 = R2_Estimator (q,y1,y2,Dy,Ds,alpha,sigma,pw);
      }

      _tx[p] = x1 ;
      _ty[p] = y1 ;
      _k [p] = x1*y2-x2*y1;
      _nx[p] = -y1 ;
      _ny[p] =  x1 ;
      _R2[p] = r2 ;
    }
    else if( coupled == 2 )
    {
/*
      _tx[p] = e / hypot( e,f ) ;
      _ty[p] = f / hypot( e,f ) ;
      _nx[p] = -_ty[p] ;
      _ny[p] =  _tx[p] ;
      _k [p] = ( h * _tx[p] - g * _ty[p] ) / c ;
*/
      // computes the functional parameters for cos(alpha)
      real nrm = hypot( g,h ) ;
      real m = (c/nrm) * ( ( f*(h/nrm) + e*(g/nrm) ) / nrm ) ;
      real n = (c/nrm) * ( ( f*(g/nrm) - e*(h/nrm) ) / nrm ) ;

      // tests the validity of the parameters
      real test = pow(fabs(m),2.0/3.0) + pow(fabs(n),2.0/3.0) ;
      if( test < 1.0 ) printf( "%f < 1!\n", (double)test ) ;

      // solves the functional: u = cos(alpha)
      real u,v ;
      bissec( m,n, u,v ) ;

//      printf( "u:%g, v:%g, c:%g, e:%g, f:%g, g:%g, h:%g, nrm:%g, m:%g,  n:%g\n", u,v,c,e,f,g,h,nrm,m,n ) ;

      // converts the result
      _tx[p] = ( h/nrm) * u + (g/nrm) * v ;
      _ty[p] = (-g/nrm) * u + (h/nrm) * v ;
      _nx[p] = -_ty[p] ;
      _ny[p] =  _tx[p] ;
      _k [p] = (nrm/c) * u ; // _tx[p] * (a*h-b*f)/quot - _ty[p] * (a*g-b*e)/quot ; // ( f * _tx[p] - e * _ty[p] ) / ( g * _tx[p] + h * _ty[p] ) ; // ( h * _tx[p] - g * _ty[p] ) / c ;
    }

    _mk += fabs(_k[p]) ;

    if( rot )
    {
      Rx = _tx[p] * rx - _ty[p] * ry ;
      Ry = _ty[p] * rx + _tx[p] * ry ;
      _tx[p] = Rx ;
      _ty[p] = Ry ;

      Rx = _nx[p] * rx - _ny[p] * ry ;
      Ry = _ny[p] * rx + _nx[p] * ry ;
      _nx[p] = Rx ;
      _ny[p] = Ry ;
    }
  }
  _mk /= _n - 2*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_nosso_curv (%scoupled, %srotated, %simproved): %d points (q=%d,alpha=%g,sigma=%g,pw=%d), %.0f msec\n", coupled?"":"not ", rot?"":"not ", improved?"":"not ", _n, q, (double)alpha, (double)sigma, pw, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;

  delete [] Ds ;
  delete [] Dx ;
  delete [] Dy ;
}
//___________________________________________________________________________________




//___________________________________________________________________________________
// Pouget based on parabola fitting
void Curve::compute_pouget_curv(const int &q)
//-----------------------------------------------------------------------------------
{
  if( q < 1 ) return ;
  reset() ;

  _mk = 0 ;
  clock_t begin = clock() ;
  for( int p = q ; p < _n-q ; ++p )
  {
    // estimate the variation of X and Y
    int i ;
    real mx = 0, my = 0,  vx = 0, vy = 0 ;
    for( i=p-q; i<=p+q; ++i )
    {
      mx += _x[i] ;
      my += _y[i] ;
      vx += _x[i] * _x[i] ;
      vy += _y[i] * _y[i] ;
    }
    mx /= (2*q+1) ;
    my /= (2*q+1) ;
    vx = vx/(2.0*q) - mx*mx * (2.0*q+1)/(2.0*q) ;
    vy = vy/(2.0*q) - my*my * (2.0*q+1)/(2.0*q) ;

    // estimate the derivatives
    real Dx, Dy ;
    real a = 0.0;
    real b = 0.0;
    real c = 0.0;
    real f = 0.0;
    real h = 0.0;
    if( vx > vy )
    {
      for( i=p-q; i<=p+q; ++i )
      {
        Dx = _x[i] - _x[p];
        Dy = _y[i] - _y[p];
        a += Dx * Dx ;
        b += Dx * Dx * Dx ;
        c += Dx * Dx * Dx * Dx ;
        f += Dx * Dy ;
        h += Dx * Dx * Dy ;
      }
      c/=4.0; b/=2.0; h/=2.0;
      real quot = a*c-b*b ;

      _ty[p] = (c*f-b*h)/quot;
      _tx[p] = sqrt( 1.0 - _ty[p] * _ty[p] ) ;
      _ny[p] = (a*h-b*f)/quot;
//      _nx[p] = sqrt( 1.0 - _ny[p] * _ny[p] ) ;
      quot = sqrt( 1.0 + _ty[p] * _ty[p] ) ;
      _k [p] = _ny[p] / ( quot * quot * quot ) ;
      _mk += fabs(_k[p]) ;
//      printf( "y(x): y = %g, x = %g\tdy = %g, d2y = %g, k = %g\n", _y[p], _x[p], _ty[p], _ny[p], _k[p] ) ;
    }
    else
    {
      for( i=p-q; i<=p+q; ++i )
      {
        Dx = _x[i] - _x[p];
        Dy = _y[i] - _y[p];
        a += Dy * Dy ;
        b += Dy * Dy * Dy ;
        c += Dy * Dy * Dy * Dy ;
        f += Dy * Dx ;
        h += Dy * Dy * Dx ;
      }
      c/=4.0; b/=2.0; h/=2.0;
      real quot = a*c-b*b ;

      _tx[p] = (c*f-b*h)/quot;
      _ty[p] = sqrt( 1.0 - _tx[p] * _tx[p] ) ;
      _nx[p] = (a*h-b*f)/quot;
//      _ny[p] = sqrt( 1.0 - _nx[p] * _nx[p] ) ;
      quot = sqrt( 1.0 + _tx[p] * _tx[p] ) ;
      _k [p] = -_nx[p] / ( quot * quot * quot ) ;
      _mk += fabs(_k[p]) ;
//      printf( "x(y): x = %g, y = %g\tdx = %g, d2x = %g, k = %g\n", _x[p], _y[p], _tx[p], _nx[p], _k[p] ) ;
    }

    if( isnan( _tx[p] ) || isnan( _ty[p] ) || isnan( _k[p] ) )
    {
      _tx[p] = _ty[p] = _nx[p] = _ny[p] = 0 ;
    }
    else
    {
      _nx[p] = -_ty[p] ;
      _ny[p] =  _tx[p] ;
    }
  }
  _mk /= _n - 2*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_pouget_curv: %d points (q=%d), %.0f msec\n", _n, q, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;
}
//___________________________________________________________________________________




//___________________________________________________________________________________
// Coeurjolly based on angle
void Curve::compute_coeurjollyI_curv(const int &q)
//-----------------------------------------------------------------------------------
{
  real vx, vy, vn,  wx, wy, wn ;
  reset() ;

  _mk = 0 ;
  clock_t begin = clock() ;
  for( int p = q ; p < _n-q ; ++p )
  {
    vx = _x[p+q] - _x[p] ;
    vy = _y[p+q] - _y[p] ;
    wx = _x[p-q] - _x[p] ;
    wy = _y[p-q] - _y[p] ;
    vn = hypot( vx,vy ) ;
    wn = hypot( wx,wy ) ;
    _k[p] = acos((vx*wx+vy*wy)/(vn*wn)) / ( vn + wn ) ;
    _mk += fabs(_k[p]) ;
  }
  _mk /= _n - 2*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_coeurjollyI_curv: %d points (q=%d), %.0f msec\n", _n, q, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
// Coeurjolly based on area
void Curve::compute_coeurjollyII_curv(const int &q)
//-----------------------------------------------------------------------------------
{
  real vx, vy, vn,  wx, wy, wn,  ux, uy, un ;
  reset() ;

  _mk = 0 ;
  clock_t begin = clock() ;
  for( int p = q ; p < _n-q ; ++p )
  {
    vx = _x[p+q] - _x[ p ] ;
    vy = _y[p+q] - _y[ p ] ;
    wx = _x[p-q] - _x[ p ] ;
    wy = _y[p-q] - _y[ p ] ;
    ux = _x[p-q] - _x[p+q] ;
    uy = _y[p-q] - _y[p+q] ;
    vn = hypot( vx,vy ) ;
    wn = hypot( wx,wy ) ;
    un = hypot( ux,uy ) ;

    _k[p] = sqrt( ( (wn+un)*(wn+un) - vn*vn ) * ( vn*vn - (wn-un)*(wn-un) ) ) / ( vn*wn*un ) ;
    _mk += fabs(_k[p]) ;
  }
  _mk /= _n - 2*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_coeurjollyII_curv: %d points (q=%d), %.0f msec\n", _n, q, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
// Worring orientation based (chain code)
void Curve::compute_worringI_curv(const int &q, const real &sigma /*= 1.0/(16*16)*/)
//-----------------------------------------------------------------------------------
{
  if( sigma == 0.0 ) printf( "compute_worringI_curv warning: sigma = 0.0!\n" ) ;
  real div = sqrt( sigma / (2*M_PI) ) ;
  reset() ;

  _mk = 0 ;
  clock_t begin = clock() ;
  for( int p = q ; p < _n-q ; ++p )
  {
    _k[p] = 0 ;
    real nrm = 0 ;
    for( int i = 1 ; i <= 2*q; ++i )
    {
      int    j = p-q+i ;
      real s = i - q - 0.5 ;
      real m = ( -s * sigma ) ;
      real e = div * exp( -s*s * sigma / 2.0 ) ;
      nrm += e ;
      _k[p] += atan2( (_y[j-1]-_y[j]), (_x[j-1]-_x[j]) ) * m * e ;
    }
    _k[p] *= 180.0 / M_PI / nrm ;
    _mk += fabs(_k[p]) ;
  }
  _mk /= _n - 2*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_worringI_curv: %d points (q=%d,sigma=%g), %.0f msec\n", _n, q, (double)sigma, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
//
int Curve::resample(const int &p, const int &q, real *&xr, real *&yr)
//-----------------------------------------------------------------------------------
{
  real minlen = FLT_MAX ;
  real len = 0;

  int i ;
  for( i=p-q ; i<p ; ++i )
  {
    real d = hypot( _x[i]-_x[i+1], _y[i]-_y[i+1] ) ;
    if(minlen > d && d > 1e-5) minlen = d ;
    len += d ;
  }
  int np = (int) floor((double)(len/minlen));

  len = 0;
  for( i=p ; i<p+q ; ++i )
  {
    real d = hypot( _x[i]-_x[i+1], _y[i]-_y[i+1] ) ;
    if(minlen > d && d > 1e-5) minlen = d ;
    len += d ;
  }
  if( np > (int) floor((double)(len/minlen)) ) np = (int) floor((double)(len/minlen)) ;

  xr = new real [2*np+1];
  yr = new real [2*np+1];

  xr[np]=_x[p];
  yr[np]=_y[p];

  int j = p-1 ;
  for( i=np-1 ; i>=0 ; --i )
  {
    real d = hypot( xr[i+1]-_x[j], yr[i+1]-_y[j] ) ;
    if( d > minlen)
    {
      real r = minlen/d ;
      xr[i] = xr[i+1] + r * (_x[j]-xr[i+1]);
      yr[i] = yr[i+1] + r * (_y[j]-yr[i+1]);
    }
    else
    {
      if( fabs(minlen-d) < (real)FLT_EPSILON )
      {
        xr[i] = _x[j-1];
        yr[i] = _y[j-1];
      }
      else
      {
        real r = minlen/(minlen-d) ;
        xr[i] = _x[j] + r * (_x[j-1]-_x[j]);
        yr[i] = _y[j] + r * (_y[j-1]-_y[j]);
      }
      --j ;
    }
  }

  j = p+1 ;
  for( i=np+1 ; i<=2*np ; ++i )
  {
    real d = hypot( xr[i-1]-_x[j], yr[i-1]-_y[j] ) ;
    if( d > minlen)
    {
      real r = minlen/d ;
      xr[i] = xr[i-1] + r * (_x[j]-xr[i-1]);
      yr[i] = yr[i-1] + r * (_y[j]-yr[i-1]);
    }
    else
    {
      if( fabs(minlen-d) < (real)FLT_EPSILON )
      {
        xr[i] = _x[j+1];
        yr[i] = _y[j+1];
      }
      else
      {
        real r = minlen/(minlen-d) ;
        xr[i] = _x[j] + r * (_x[j+1]-_x[j]);
        yr[i] = _y[j] + r * (_y[j+1]-_y[j]);
      }
      ++j ;
    }
  }
  return np ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
// Worring orientation based (resampling)
void Curve::compute_worringII_curv(const int &q, const real &sigma /*= 1.0/(16*16)*/)
//-----------------------------------------------------------------------------------
{
  if( sigma == 0.0 ) printf( "compute_worringI_curv warning: sigma = 0.0!\n" ) ;
  real div = sqrt( sigma / (2*M_PI) ) ;
  reset() ;

  _mk = 0 ;
  clock_t begin = clock() ;
  for( int p = q ; p < _n-q ; ++p )
  {
    real *xr, *yr ;
    real  nrm = 0 ;
    int     np  = resample( p, q, xr, yr ) ;
    for( int i = 1 ; i <= 2*np; ++i )
    {
      real j = i - np - 0.5 ;
      real m = ( -j * sigma ) ;
      real e = div * exp( -j*j * sigma / 2.0 ) ;
      nrm += e ;
      _k[p] += atan2( (yr[i-1]-yr[i]), (xr[i-1]-xr[i]) ) * m * e ;
    }
    _nx[p] = _ny[p] = _tx[p] = _ty[p] = 0 ;
    _k[p] *= 180.0 / M_PI / 1.107 / nrm ;
    _mk += fabs(_k[p]) ;
    delete [] xr ;
    delete [] yr ;
  }
  _mk /= _n - 2*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_worringII_curv: %d points (q=%d,sigma=%g), %.0f msec\n", _n, q, (double)sigma, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
// Worring based on line fitting
void Curve::compute_worringIII_curv(const int &q, const real &sigma /*= 1.0/(16*16)*/)
//-----------------------------------------------------------------------------------
{
  double div = sqrt( (double)sigma / (2.0*M_PI) ) ;
  reset() ;

  double *xc = new double[2*q-1] ;
  double *yc = new double[2*q-1] ;
  double *w  = new double[2*q-1] ;

  _mk = 0 ;
  clock_t begin = clock() ;
  for( int j = -q+1 ; j < q-1 ; ++j )
  {
    w[q-1+j] = div * exp( -j*j * (double)sigma / 2.0 ) ;
  }

  for( int p = q ; p < _n-q ; ++p )
  {
    int i ;

    // estimate the orientation at p-1
    double cl0, cl1 ; // left
    for( i = 0 ; i < 2*q-1 ; ++i )
    {
      xc[i] = _x[p-q+i] ;
      yc[i] = _y[p-q+i] ;

      double cov00, cov01, cov11, chisq;
      gsl_fit_wlinear (xc, 1, w, 1, yc, 1, 2*q-1, &cl0, &cl1, &cov00, &cov01, &cov11, &chisq) ;
    }

    // estimate the orientation at p+1
    double cr0, cr1 ; // right
    for( i = 0 ; i < 2*q-1 ; ++i )
    {
      xc[i] = _x[p-q+i+2] ;
      yc[i] = _y[p-q+i+2] ;

      double cov00, cov01, cov11, chisq;
      gsl_fit_wlinear (xc, 1, w, 1, yc, 1, 2*q-1, &cr0, &cr1, &cov00, &cov01, &cov11, &chisq) ;
    }

    if( fabs(cr1) > 2 && fabs(cl1) > 2 )
      _k[p] = ( atan2( 1.0,cr1 ) - atan2( 1.0,cl1 ) ) / hypot( _x[p+1]-_x[p-1], _y[p+1]-_y[p-1] ) ;
    else
      _k[p] = ( atan( cr1 ) - atan( cl1 ) ) / hypot( _x[p+1]-_x[p-1], _y[p+1]-_y[p-1] ) ;
    _mk += fabs(_k[p]) ;
  }
  _mk /= _n - 2*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_worringII_curv: %d points (q=%d,sigma=%g), %.0f msec\n", _n, q, (double)sigma, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;

  delete [] xc ;
  delete [] yc ;
  delete [] w  ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
// Worring based on path
void Curve::compute_worringIV_curv(const int &q, const real &sigma /*= 1.0/(16*16)*/)
//-----------------------------------------------------------------------------------
{
  if( sigma == 0.0 ) printf( "compute_worringII_curv warning: sigma = 0.0!\n" ) ;
  reset() ;

  real D, Dx, Dy, DDx, DDy ;

  _mk = 0 ;
  clock_t begin = clock() ;
  for( int p = q ; p < _n-q ; ++p )
  {
    D = Dx = Dy = DDx = DDy = 0 ;
    for( int i = -q ; i <= q; ++i )
    {
      int    r = p+i ;
      real m = ( -(real)i * sigma ) ;
      real e = exp( (real)(-i*i) * sigma / 2.0 ) ;

      D   += e ;
      Dx  += _x[r] * m * e ;
      Dy  += _y[r] * m * e ;
      DDx += _x[r] * (m*m - sigma) * e ;
      DDy += _y[r] * (m*m - sigma) * e ;
    }
    Dx  /= D ;
    Dy  /= D ;
    DDx /= D ;
    DDy /= D ;

    real nrm = hypot( Dx,Dy ) ;
    _k [p] = -( Dx * DDy - Dy * DDx ) / ( nrm * nrm * nrm ) ;
    _mk += fabs(_k[p]) ;
    nrm *= sign_of( _k[p] ) ;
    _nx[p] = ( DDx * Dy * Dy - Dx * Dy * DDy ) / ( nrm * nrm * nrm ) ;
    _ny[p] = ( DDy * Dx * Dx - Dy * Dx * DDx ) / ( nrm * nrm * nrm ) ;

    _tx[p] = Dx / nrm ;
    _ty[p] = Dy / nrm ;
  }
  _mk /= _n - 2*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_worringIV_curv: %d points (q=%d,sigma=%g), %.0f msec\n", _n, q, (double)sigma, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
// Utcke based on circle fitting
void Curve::compute_utcke_curv(const int &q)
//-----------------------------------------------------------------------------------
{
  gsl_set_error_handler_off();

  gsl_matrix *M = gsl_matrix_alloc(2*q+1,4);
  gsl_matrix *A = gsl_matrix_alloc(4,4);
  reset() ;

  _mk = 0 ;
  clock_t begin = clock() ;
  for( int p = q ; p < _n-q ; ++p )
  {
    for( int i = 0 ; i <= 2*q; ++i )
    {
      real Dx = _x[p-q+i] ;
      real Dy = _y[p-q+i] ;
      gsl_matrix_set(M,i,0,1);
      gsl_matrix_set(M,i,1,Dx);
      gsl_matrix_set(M,i,2,Dy);
      gsl_matrix_set(M,i,3,Dx*Dx + Dy*Dy);
    }

    gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0,M, M, 0.0, A );
    if( gsl_linalg_cholesky_decomp( A ) == GSL_EDOM ) { _k[p] = nan() ; return ; }

    real a = gsl_matrix_get(A,0,0) * gsl_matrix_get(A,1,1) * gsl_matrix_get(A,2,2) ;
    real b = gsl_matrix_get(A,0,0) * gsl_matrix_get(A,1,1) * gsl_matrix_get(A,2,3) ;
    real c = gsl_matrix_get(A,0,0) * ( gsl_matrix_get(A,1,2) * gsl_matrix_get(A,2,3) -
                                         gsl_matrix_get(A,2,2) * gsl_matrix_get(A,1,3) ) ;
    real d =
      gsl_matrix_get(A,0,1) * ( gsl_matrix_get(A,1,2) * gsl_matrix_get(A,2,3) -
                                gsl_matrix_get(A,2,2) * gsl_matrix_get(A,1,3) ) +
      gsl_matrix_get(A,1,1) * ( gsl_matrix_get(A,0,3) * gsl_matrix_get(A,2,2) -
                                gsl_matrix_get(A,0,2) * gsl_matrix_get(A,2,3) ) ;

    b /=  2.0*a ;
    c /= -2.0*a ;
    d /=     a  ;
    _k[p] = 1.0 / sqrt( b*b + c*c + d ) ;

    _nx[p] = c - _x[p] ;
    _ny[p] = b - _y[p] ;

    if( (double) (_nx[p] * (_y[p+1]-_y[p]) - _ny[p] * (_x[p+1]-_x[p])) > 0 )
      _k[p] *= -1.0 ;
    _mk += fabs(_k[p]) ;

    d = sign_of(_k[p] ) * hypot( _nx[p],_ny[p] ) ;
    _nx[p] /= d ;
    _ny[p] /= d ;
    _tx[p] =  _ny[p] ;
    _ty[p] = -_nx[p] ;
  }
  _mk /= _n - 2*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_utcke_curv: %d points (q=%d), %.0f msec\n", _n, q, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;

  gsl_matrix_free(M) ;
  gsl_matrix_free(A) ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
// Belyaev via derivatives
void Curve::compute_belyaev_curv(const int &q)
//-----------------------------------------------------------------------------------
{
  reset() ;

  _mk = 0 ;
  clock_t begin = clock() ;
  for( int p = q ; p < _n-q ; ++p )
  {
    real vx = _x[p] - _x[p-q] ;
    real vy = _y[p] - _y[p-q] ;
    real ux = _x[p+q] - _x[p] ;
    real uy = _y[p+q] - _y[p] ;
    real a  = hypot( vx,vy ) ;
    real b  = hypot( ux,uy ) ;

    _tx[p] = ux/b + vx/a - (ux+vx)/(a+b) ;
    _ty[p] = uy/b + vy/a - (uy+vy)/(a+b) ;
    _nx[p] = 2.0 * (_x[p-q]/((a+b)*a) - _x[p]/(b*a) + _x[p+q]/((a+b)*b) ) ;
    _ny[p] = 2.0 * (_y[p-q]/((a+b)*a) - _y[p]/(b*a) + _y[p+q]/((a+b)*b) ) ;

    real nrm = hypot( _tx[p],_ty[p] ) ;
    _k [p] = ( _tx[p] * _ny[p] - _ty[p] * _nx[p] ) / ( nrm * nrm * nrm ) ;
    _mk += fabs(_k[p]) ;

    _tx[p] /= nrm ;
    _ty[p] /= nrm ;
    nrm = sign_of( _k[p] ) * hypot( _nx[p],_ny[p] ) ;
    _nx[p] /= nrm ;
    _ny[p] /= nrm ;
  }
  _mk /= _n - 2*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_belyaev_curv: %d points (q=%d), %.0f msec\n", _n, q, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
// Gumhold based on external angles
void Curve::compute_gumhold_curv(const int &q, const real &alpha /*= 1.0*/)
//-----------------------------------------------------------------------------------
{
  if( alpha == 0.0 ) printf( "compute_gumhold_curv warning: alpha = 0.0!\n" ) ;
  reset() ;

  _mk = 0 ;
  clock_t begin = clock() ;
  for( int p = q ; p < _n-q ; ++p )
  {
    real vx = _x[p] - _x[p-q] ;
    real vy = _y[p] - _y[p-q] ;
    real ux = _x[p+q] - _x[p] ;
    real uy = _y[p+q] - _y[p] ;
    real a  = hypot( vx,vy ) ;
    real b  = hypot( ux,uy ) ;
    real t = acos( (vx*ux + vy*uy) / (a*b) ) ;

    _k [p] = 2.0 * t / (a+b) / cos( t/(2.0+alpha/10.0) ) ;
    _mk += fabs(_k[p]) ;
  }
  _mk /= _n - 2.0*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_gumhold_curv: %d points (q=%d,alpha=%g), %.0f msec\n", _n, q, (double)alpha, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;
}
//___________________________________________________________________________________



//___________________________________________________________________________________
// Estrozi based on FFT
#define REAL(z,i) ((z)[ 2*(i) ])
#define IMAG(z,i) ((z)[2*(i)+1])
void Curve::compute_estrozi_curv(const int &q)
//-----------------------------------------------------------------------------------
{
  reset() ;

  int nq = 2*q+1 ;
  double *u1data = new double[2*nq] ;
  double *u2data = new double[2*nq] ;

  gsl_fft_complex_wavetable *wavetable = gsl_fft_complex_wavetable_alloc( nq ) ;
  gsl_fft_complex_workspace *workspace = gsl_fft_complex_workspace_alloc( nq ) ;

  _mk = 0 ;
  clock_t begin = clock() ;
  int i ;
  for( int p = q ; p < _n-q ; ++p )
  {
    // fills the signals
    for( i = 0 ; i <= 2*q ; ++i )
    {
      REAL(u1data,i) = _x[p-q+i] ;
      IMAG(u1data,i) = _y[p-q+i] ;
    }

    int res = gsl_fft_complex_forward( u1data, 1, nq, wavetable, workspace ) ;
    if( res == GSL_EDOM || res == GSL_EINVAL )  { _k[p] = nan() ; return ; }
    memcpy( u2data, u1data, 2*nq*sizeof(double) ) ;

    real re,im ;
    for( i = 0 ; i <= 2*q ; ++i )
    {
      // multiplication by 2.i.Pi.s
      double d = 2*M_PI*(i-q) ;
      re = REAL(u1data,i) ;
      im = IMAG(u1data,i) ;
      REAL(u1data,i) =  -d * im ;
      IMAG(u1data,i) =   d * re ;

      // multiplication by -(2.Pi.s)^2
      d *= -d ;
      REAL(u2data,i) *= d ;
      IMAG(u2data,i) *= d ;
    }
    res = gsl_fft_complex_backward( u1data, 1, nq, wavetable, workspace ) ;
    if( res == GSL_EDOM || res == GSL_EINVAL )  { _k[p] = nan() ; return ; }
    res = gsl_fft_complex_backward( u2data, 1, nq, wavetable, workspace ) ;
    if( res == GSL_EDOM || res == GSL_EINVAL )  { _k[p] = nan() ; return ; }

    _tx[p] = REAL(u1data,q) ;
    _ty[p] = IMAG(u1data,q) ;
    _nx[p] = REAL(u2data,q) ;
    _ny[p] = IMAG(u2data,q) ;
    real nrm = hypot( _tx[p],_ty[p] ) ;
    _k [p] = -( _tx[p] * _ny[p] - _ty[p] * _nx[p] ) / ( nrm * nrm * nrm ) ;
    _mk += fabs(_k[p]) ;
    _tx[p] /= nrm ;
    _ty[p] /= nrm ;
  }
  _mk /= _n - 2*q ;
  clock_t end = clock() ;
  printf( "Curve::compute_estrozi_curv: %d points (q=%d), %.0f msec\n", _n, q, (double)(end - begin) * 1000.0 / CLOCKS_PER_SEC ) ;

  delete [] u1data ;
  delete [] u2data ;

  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);
}
//___________________________________________________________________________________




//___________________________________________________________________________________
//___________________________________________________________________________________




//___________________________________________________________________________________
//
int Curve::refine()
//-----------------------------------------------------------------------------------
{
  int i,o_n = _n ;
  vector<real> o_x  = _x , o_y  = _y  ;      // sample points
  vector<real> o_tx = _tx, o_ty = _ty ;      // tangents   of sample points
  vector<real> o_nx = _nx, o_ny = _ny ;      // normals    of sample points
  vector<real> o_k  = _k              ;      // curvatures of sample points

  set_n( 2*_n-1 ) ;
  for( i = 0 ; i < o_n ; ++i )
  {
    int j = 2*i ;
    _x [j] = o_x [i] ;
    _y [j] = o_y [i] ;
    _tx[j] = o_tx[i] ;
    _ty[j] = o_ty[i] ;
    _nx[j] = o_nx[i] ;
    _ny[j] = o_ny[i] ;
    _k [j] = o_k [i] ;
  }

  for( i = 0 ; i < o_n-1 ; ++i )
  {
    int j  = 2*i+1 ;
    int i1 = i - (i?1:0) ;
    int i2 = i + (i==o_n-2?1:2) ;
    _x [j] = ( 9.0 * ( o_x [i] + o_x [i+1] ) - o_x [i1] - o_x [i2] ) / 16.0 ;
    _y [j] = ( 9.0 * ( o_y [i] + o_y [i+1] ) - o_y [i1] - o_y [i2] ) / 16.0 ;
    _tx[j] = ( 9.0 * ( o_tx[i] + o_tx[i+1] ) - o_tx[i1] - o_tx[i2] ) / 16.0 ;
    _ty[j] = ( 9.0 * ( o_ty[i] + o_ty[i+1] ) - o_ty[i1] - o_ty[i2] ) / 16.0 ;
    _nx[j] = ( 9.0 * ( o_nx[i] + o_nx[i+1] ) - o_nx[i1] - o_nx[i2] ) / 16.0 ;
    _ny[j] = ( 9.0 * ( o_ny[i] + o_ny[i+1] ) - o_ny[i1] - o_ny[i2] ) / 16.0 ;
    _k [j] = ( 9.0 * ( o_k [i] + o_k [i+1] ) - o_k [i1] - o_k [i2] ) / 16.0 ;
  }

  return _n ;
}
//___________________________________________________________________________________


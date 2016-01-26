// Copyright (C) 2010  Michael Seyfert <m@codesand.org>
/*
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef LAS_MTYPES_H_
#define LAS_MTYPES_H_

#include<gmpxx.h>

namespace las{

// Large integer.
typedef mpz_class bigint;

// Large rational.
typedef mpq_class bigq;

// Large float.
typedef mpf_class bigf;

//============================================================
// Float
// Floating point class with epsilong checking.
//
// TODO: Plan on getting rid of this entirely.
// It is a messy way to handle floating point algorithms and inaccuracies.
//============================================================
static const bigf F_EPS(FLOAT_ZERO_EPS);
class Float
{
public:
	bigf f;

	Float(const char* p):f(p){}
	Float(const bigf &val = bigf()):f(val){}
	Float(int i):f(i){}
	Float(long l):f(l){}
	Float(double d):f(d){}

	bool operator==(const Float &r)const{ return f > r.f-F_EPS && f < r.f+F_EPS; }
	bool operator==(int i)const{return *this == Float(i);}
	bool operator!=(const Float &r)const{ return !(*this == r); }
	bool operator<(const Float &r)const{ return  f < r.f; }
	bool operator<=(const Float &r)const{ return  f <= r.f; }
	bool operator>(const Float &r)const{ return  f > r.f; }
	bool operator>=(const Float &r)const{ return  f >= r.f; }

	Float operator*(const Float &r)const{return bigf(f*r.f);}
	Float operator/(const Float &r)const{return bigf(f/r.f);}
	Float operator+(const Float &r)const{return bigf(f+r.f);}
	Float operator-(const Float &r)const{return bigf(f-r.f);}
	Float operator%(const Float &r)const{return bigf(f - r.f * floor(f / r.f));}
	Float& operator+=(const Float &r){ f += r.f; return *this; }
	Float& operator-=(const Float &r){ f -= r.f; return *this; }
	Float& operator/=(const Float &r){ f /= r.f; return *this; }
	Float& operator*=(const Float &r){ f *= r.f; return *this; }
	Float& operator%=(const Float &r){ f -= r.f * floor(f / r.f); return *this; }
	Float operator-()const{return bigf(-f);}

	void Canonicalize(){
		if(-F_EPS < f && f < F_EPS) f = 0;
	}

private:

	// Do not allow comparison with the base type.
	bool operator==(bigf)const;
	bool operator!=(bigf)const;
};

inline std::istream& operator>>(std::istream &is,Float &t)
{
	return is >> t.f;
}

inline std::ostream& operator<<(std::ostream &os,const Float &t)
{
	return os << t.f;
}

// abs for Float.
inline Float abs(const Float &f){
	return bigf(abs(f.f));
}

// sqrt for Float.
inline Float sqrt(const Float &f){
	return bigf(sqrt(f.f));
}

// operator (int) * (Float)
inline Float operator*(int integer,const Float &f){
	return bigf(integer * f.f);
}
//============================================================

/*
Canonicalize will zero out small floating point values.
For rational types it will remove common factors from the numerator and
denominator (Reduce them, which is necessary for libgmp).
This is useful for input sanitization.
 */
template<class T> void Canonicalize(T&);

template<> void Canonicalize<bigq>(bigq &c)
{
	c.canonicalize();
}

template<> void Canonicalize<Float>(Float &c)
{
	c.Canonicalize();
}

template<class T>
void Canonicalize(std::complex<T> &c)
{
    // Hack to get real/imaginary part as a reference (since c++11).
    T &r = reinterpret_cast<T(&)[2]>(c)[0];
	Canonicalize(r);
    T &i = reinterpret_cast<T(&)[2]>(c)[1];
	Canonicalize(i);
}

// Canonicalize every element of a matrix.
template<class T>
void Canonicalize(boost::numeric::ublas::matrix<T> &m);
template<class T>
void Canonicalize(boost::numeric::ublas::triangular_matrix<T,boost::numeric::ublas::upper> &m);
//============================================================


// === Rational Matrix types.
typedef std::complex<bigq> q_type;
typedef boost::numeric::ublas::matrix<q_type> mxq; // rational w/ imaginary matrix
typedef boost::numeric::ublas::matrix<bigq> mxqr; // rational real matrix.
typedef boost::numeric::ublas::triangular_matrix<q_type,boost::numeric::ublas::upper> mxqtu;
// ===

// === Floating point matrix types.
typedef std::complex<Float> f_type;
typedef boost::numeric::ublas::matrix<f_type> mxf; // float w/ imaginary matrix
typedef boost::numeric::ublas::matrix<Float> mxfr; // float real matrix.
typedef boost::numeric::ublas::triangular_matrix<f_type,boost::numeric::ublas::upper> mxftu;
// ===

// Default matrix.
#if defined(MATRIX_RATIONAL)
typedef mxq mx;
typedef mxqtu mxtu;
typedef q_type mx_value_type;
typedef bigq scalar_type;
#elif defined(MATRIX_FLOAT)
typedef mxf mx;
typedef mxftu mxtu;
typedef f_type mx_value_type;
typedef Float scalar_type;
#endif

} //namespace las

// Define complex type traits.
// So our complex types will work with libboost.
namespace boost{ namespace numeric{ namespace ublas{

template<>
struct type_traits<std::complex<las::Float> > : complex_traits<std::complex<las::Float> >{
	typedef type_traits<std::complex<las::Float> > self_type;
	typedef std::complex<las::Float> value_type;
	typedef const value_type &const_reference;
	typedef value_type &reference;
	typedef las::Float real_type;
	typedef std::complex<las::Float> precision_type;
};
template<>
struct type_traits<std::complex<las::bigq> > : complex_traits<std::complex<las::bigq> >{
	typedef type_traits<std::complex<las::bigq> > self_type;
	typedef std::complex<las::bigq> value_type;
	typedef const value_type &const_reference;
	typedef value_type &reference;
	typedef las::bigq real_type;
	typedef std::complex<las::bigq> precision_type;
};

namespace detail {
// libgmp types do not have a simple constructor.
template<>
struct has_trivial_constructor<std::complex<las::Float> > : public boost::false_type {};
template<>
struct has_trivial_destructor<std::complex<las::Float> > : public boost::false_type {};

template<>
struct has_trivial_constructor<las::Float> : public boost::false_type {};
template<>
struct has_trivial_destructor<las::Float> : public boost::false_type {};

template<>
struct has_trivial_constructor<std::complex<las::bigq> > : public boost::false_type {};
template<>
struct has_trivial_destructor<std::complex<las::bigq> > : public boost::false_type {};

template<>
struct has_trivial_constructor<las::bigq> : public boost::false_type {};
template<>
struct has_trivial_destructor<las::bigq> : public boost::false_type {};

} // namespace detail
}}} // namespace boost, numeric, ublas

#include"mtypes-inl.h"
#endif  //LAS_MTYPES_H_

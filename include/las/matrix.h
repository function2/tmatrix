// Copyright (C) 2012  Michael Seyfert <m@codesand.org>
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
#ifndef LAS_MATRIX_H_
#define LAS_MATRIX_H_

#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/triangular.hpp>

namespace las{

/*
Read a matrix from stream 'in'.
The matrix is entered one row at a time, with each element whitespace separated.
Each row is separated by a newline.
Enter a blank row to finish.
 */
template<class T>
std::istream& MatrixInput(
	std::istream &in
	,boost::numeric::ublas::matrix<T> &m);

/*
Write a matrix to stream 'out'.
 */
template<class M1>
std::ostream& MatrixOutput(
	std::ostream &os,
	const M1 &m);

// Returns the identity matrix of size n.
template<class T>
inline boost::numeric::ublas::identity_matrix<T> IdentityMatrix(size_t n)
{
	return boost::numeric::ublas::identity_matrix<T>(n);
}

// Take A to the nth power. A must be square, n >= 0.
// Return the result.
template<class M1,class M2>
M2& Pow(const M1 &a,long n,M2 &r);

// Compute e^A using terms of the power series expansion.
template<class M1,class M2>
M2& Exp(const M1 &a,long terms,M2 &r);

// Returns true if the matrices are equal
// (they have the same width, height, and the elements are equal).
template<class M1,class M2>
bool MatrixEquality(const M1 &m1,const M2 &m2);

// Return the number of row exchanges needed to convert this permutation
// matrix to the identity.
template<class M1>
int NumRowExchanges(const M1 &perm);

// Returns true if the matrix is square.
template<class M1>
inline bool IsSquare(const M1 &m);

// Returns true if the matrix is symmetric.
template<class M1>
bool IsSymmetric(const M1 &a);

// Swaps rows
template<class M1>
void SwapRows(M1 &m,size_t r1,size_t r2);

// Swaps cols
template<class M1>
void SwapCols(M1 &m,size_t c1,size_t c2);

// Divides a row by a scalar
template<class M1,class T>
void DivRow(M1 &m,size_t row,const T &val);

// Divides a column by a scalar
template<class M1,class T>
void DivCol(M1 &m,size_t col,const T &val);

// Multiply a column by a scalar
template<class M1,class T>
void MultCol(M1 &m,size_t col,const T &val);

// Subtracts val * (src row) from (dest row).
template<class M1,class T>
void MultSubRow(M1 &m,size_t src,size_t dest,const T &val);

// Subtracts val * (src col) from (dest col).
template<class M1,class T>
void MultSubCol(M1 &m,size_t src,size_t dest,const T &val);

// Creates a matrix from the given rows.
template<class M1>
static void GetSubMatrixRows(
	const M1 &m,const std::vector<size_t> &rows,M1 &res);

// Creates a submatrix from [low row, high row]
template<class M1>
static void GetSubMatrixRows(
	const M1 &m,int lowrow,int highrow,M1 &res);

// Creates a matrix from the given columns
template<class M1>
static void GetSubMatrixCols(
	const M1 &m,const std::vector<size_t> &cols,M1 &res);

// Find the Rx = c form of this system Ax = b.
// R is row echelon form.
// b may have multiple columns.
// r,c do not need to be initialized.
template<class M1,class M2,class M3,class M4>
void RRef(const M1 &a,const M2 &b, M3 &r, M4 &c);
template<class M1,class M2>
void RRef(const M1 &a,M2 &r);

// Returns the rank of a matrix
// The input matrix must be in RRef form.
template<class M1>
size_t Rank(const M1 &a);

// Returns a matrix built of every nonzero row.
template<class M1>
void GetNonZeroRows(const M1 &a,M1 &out);

// Returns the nullspace of this matrix.
// The input matrix 'a' must be in RRef form.
template<class M1>
void Nullspace(const M1 &a,M1 &n);
// This will also return the pivot columns of the matrix.
template<class M1>
void Nullspace(const M1 &a,M1 &n,std::vector<size_t> &pivot_cols);
// If you know the rank, you can pass it to speed up this operation.
template<class M1>
void Nullspace(const M1 &a,M1 &n,size_t rank,std::vector<size_t> &pivot_cols);

// Returns pivot columns (must be in rref form)
template<class M1>
void GetPivotCols(const M1 &a,std::vector<size_t> &pivot_cols);

/*
Computes the inverse of A.
Returns false if A is rectangular or there is no inverse.
inv does not need to be initialized.
 */
template<class M1>
bool Inverse(const M1 &a,M1 &inv);

// Finds a PA = LU form of matrix A.
template<class M1>
void LU(const M1 &a,M1 &p,M1 &l,M1 &u);

// Returns PA = LDU form,
// Returns false if no such factorization exists.
template<class M1>
bool LDU(const M1 &a,M1 &p,M1 &l,M1 &d,M1 &u);

// Solve the equation Ax=b. x will contain only pivot column values.
// (If the system has infinite solutions, this will only return one.)
template<class M1>
bool SolvePivots(const M1 &a,const M1 &b,M1 &x);
// Same as above but using Rx = c inputs (see RRef)
template<class M1>
bool SolvePivotsR(const M1 &r,const M1 &c,M1 &x);

/*
A = QR Factorization.
The columns of Q are orthogonal.
R will be upper triangular.

A = UR Factorization.
U is unitary.
Returns false if the columns are not independent.
 */
template<class M1,class M2,class M3>
bool QR_Factor(const M1 &a,M2 &q,M3 &r);
template<class M1,class M2,class M3>
bool UR_Factor(const M1 &a,M2 &q,M3 &r);
template<class M1,class M2,class M3>
bool QR_FactorNormal(const M1 &a,M2 &q,M3 &r); // make Q orthonormal.
template<class M1,class M2,class M3>
bool UR_FactorNormal(const M1 &a,M2 &q,M3 &r); // make U orthonormal.

// Get the determinant of matrix A
// This will compute the PA = LDU factorization and return them if wanted.
// The determinant is the product of
// determinant(p) x (product of diagonals of d).
// INPUT matrix A must be square.
// INPUT p,l,d,u should point to matrices stored by the caller, or null
// FIXME this is slow (numrowexchanges)
template<class M1>
typename M1::value_type Determinant(const M1 &a,M1 *p=NULL,M1 *l=NULL,M1 *d=NULL,M1 *u=NULL);

// Get the minor matrix at location (r,c)
// This is matrix A with row r and column c removed.
template<class M1,class M2>
void MinorMatrix(const M1 &a,M2 &m,size_t row,size_t col);

// Get the cofactor matrix of A
template<class M1,class M2>
void GetCofactorMatrix(const M1 &a,M2 &cofactor);

} //namespace las

#include"matrix-inl.h"

#endif  // LAS_MATRIX_H_

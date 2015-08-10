// Copyright (C) 2008  Michael Seyfert <michael@codesand.org>
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

/*
Compute the 4 Fundamental Subspaces of a matrix.

Rowspace
Nullspace
Left nullspace
Columnspace
 */

#include<iostream>
#include<boost/format.hpp>

#include"las/matrix.h"
#include"las/mtypes.h"

using namespace std;
using namespace las;
using boost::format;

int main(int argc, char**argv)
{
#if defined(MATRIX_FLOAT)
	mpf_set_default_prec(MATRIX_FLOAT_PREC);
#endif  //MATRIX_FLOAT

	mx a;
	printf("Input matrix A:\n");
	MatrixInput(cin,a);
	Canonicalize(a);

	//b starts as the identity matrix of size a.Rows()
	mx b(IdentityMatrix<mx::value_type>(a.size1()));

	//reduce A to R form
	mx r,c;
	RRef(a,b,r,c);

	int rank;
	mx ca,cat,n,natp,nat;
	rank = Rank(r);

	//N(A) = nullspace
	vector<size_t> pivot_cols;
	Nullspace(r,n,rank,pivot_cols);

	//C(A) = column space
	// the pivot columns form a basis for the column space of A.
	GetSubMatrixCols(a,pivot_cols,ca);

	//C(A^T) = Rowspace
	GetNonZeroRows(r,cat);

	//N(A^T) = left nullspace
	GetSubMatrixRows(c,rank,r.size1()-1,natp);
	nat = trans(natp);

#if defined(MATRIX_FLOAT)
	Canonicalize(ca);
	Canonicalize(n);
	Canonicalize(cat);
	Canonicalize(nat);
#endif  //MATRIX_FLOAT

	cout << format("C(A) {rank = %1%} = ") % rank;
	MatrixOutput(cout,ca);

	cout << format("N(A) {rank = %1%} = ") % (a.size2()-rank);
	MatrixOutput(cout,n);

	cout << format("C(A^T) {rank = %1%} = ") % rank;
	MatrixOutput(cout,trans(cat));

	cout << format("N(A^T) {rank = %1%} = ") % (a.size1()-rank);
	MatrixOutput(cout,nat);
}

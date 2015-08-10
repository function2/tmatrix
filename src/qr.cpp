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
A = QR Factorization
For rationals the columns of Q will be orthogonal but not orthonormal.

A = UR Factorization.
U is a unitary matrix.
 */
#include<iostream>

#include"las/matrix.h"
#include"las/mtypes.h"

using namespace std;
using namespace las;

int main(int argc, char**argv)
{
#if defined(MATRIX_FLOAT)
	mpf_set_default_prec(MATRIX_FLOAT_PREC);
#endif  //MATRIX_FLOAT

	mx a;
	cout << "Input matrix A (empty line to finish):\n";
	MatrixInput(cin,a);
	Canonicalize(a);

	mx q;
	mxtu r; // Use a triangular matrix to save space.

	bool check;
#if defined(MATRIX_FLOAT)
	check = QR_FactorNormal(a,q,r);
#else
	check = QR_Factor(a,q,r);
#endif
	if(!check)
		printf("Warning: Columns not independent.\n");

#if defined(MATRIX_FLOAT)
	Canonicalize(q);
	Canonicalize(r);
#endif  //MATRIX_FLOAT

	printf("q = "); MatrixOutput(cout,q);
	printf("r = "); MatrixOutput(cout,r);

	// Check the result.
	printf("Check, Does a = qr? ");
	cout << flush;
	if(MatrixEquality(prod(q,r), a)) {
		printf("Yes\n");
	}else
		printf("NO, ERROR\n");

	// A = UR
	mx u;
#if defined(MATRIX_FLOAT)
	check = UR_FactorNormal(a,u,r);
#else
	check = UR_Factor(a,u,r);
#endif
	if(!check)
		printf("Warning: Columns not independent.\n");

#if defined(MATRIX_FLOAT)
	Canonicalize(u);
	Canonicalize(r);
#endif  //MATRIX_FLOAT

	printf("u = "); MatrixOutput(cout,u);
	printf("r = "); MatrixOutput(cout,r);

	// Check the result.
	printf("Check, Does a = ur? ");
	cout << flush;
	if(MatrixEquality(prod(u,r), a)) {
		printf("Yes\n");
	}else
		printf("NO, ERROR\n");

	// Further check.
	mx c1 = prod(trans(q),q);
	mx c2 = prod(herm(u),u);
#if defined(MATRIX_FLOAT)
	Canonicalize(c1);
	Canonicalize(c2);
#endif  //MATRIX_FLOAT
	printf("trans(q) * q = "); MatrixOutput(cout,c1);
	printf("herm(u) * u = "); MatrixOutput(cout,c2);
}

// Copyright (C) 2009  Michael Seyfert <michael@codesand.org>
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
Find the eigenvalues and eigenvectors of a 2 by 2 matrix.
 */

#include<iostream>
#include<algorithm>

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
	printf("Input matrix A (empty line to finish):\n");
	MatrixInput(cin,a);
//	Canonicalize(a);

	if(!IsSquare(a) || a.size1() != 2){
		printf("Matrix is not square 2 by 2.\n");
		return EXIT_FAILURE;
	}

	printf("a = ");
	MatrixOutput(cout,a);

//========
	typedef mx::value_type var;

	// Solve for the eigenvalues.
	var b = -(a(0,0) + a(1,1)); // b = -trace
	var c = a(0,0) * a(1,1) - a(1,0) * a(0,1); // c = det A
	Canonicalize(b);
	Canonicalize(c);

	var d = -b / var(2);
	var e = sqrt(b*b - var(4)*c) / var(2);

	var lambda1 = d+e;
	var lambda2 = d-e;

	cout << format("lambda_1 = %1%, lambda_2 = %2%\n") % lambda1 % lambda2;

	mx x1,x2;
	mx m1 = a - var(lambda1) * IdentityMatrix<var>(2);
	mx m2 = a - var(lambda2) * IdentityMatrix<var>(2);
	mx r1,r2;
	RRef(m1,r1); Nullspace(r1,x1);
	RRef(m2,r2); Nullspace(r2,x2);

	cout << "x1 = "; MatrixOutput(cout,x1);
	cout << "x2 = "; MatrixOutput(cout,x2);
	cout << format("sum of eigenvalues = %1%\n") % (lambda1 + lambda2);
	cout << format("product of eigenvalues = %1%\n") % (lambda1 * lambda2);

	// S^-1 A S = Lambda factorization.
	mx s(2,2);
	s(0,0) = x1(0,0);
	s(1,0) = x1(1,0);
	s(0,1) = x2(0,0);
	s(1,1) = x2(1,0);

	mx si;

	if(!Inverse(s,si)){
		cout << "No S^-1 A S factorization.\n";
	}else{

		mx pp = prod(si,a);
		mx l = prod(pp,s);
		if(!MatrixEquality(l,boost::numeric::ublas::zero_matrix<mx::value_type>(2))){
			Canonicalize(s);
			Canonicalize(si);
			cout << "s^-1 = "; MatrixOutput(cout,si);
			cout << "s = "; MatrixOutput(cout,s);
			Canonicalize(l);
			cout << "Lambda = "; MatrixOutput(cout,l);
		}
	}
}

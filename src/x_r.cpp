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
compute x_r and x_n
where x = x_r + x_n
x_r is in the rowspace of A
x_n is in the nullspace of A
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

	mx a,x;

	printf("Input matrix A:\n");
	MatrixInput(cin,a);
	Canonicalize(a);

	printf("Enter x:\n");
	MatrixInput(cin,x);
	Canonicalize(x);

	if(x.size2() != 1 || x.size1() != a.size2()){
		printf("Invalid entry.\n");
		exit(EXIT_FAILURE);
	}

	//to find x_r we solve A(A^T v + x_n) = Ax
	//where x_r = A^T * v
	mx aat = prod(a,trans(a));
	mx ax = prod(a,x);

	//solve
	mx v;
	SolvePivots(aat,ax,v);

	mx x_r = prod(trans(a),v);
	mx x_n = x - x_r;

#if defined(MATRIX_FLOAT)
	Canonicalize(x_r);
	Canonicalize(x_n);
#endif  //MATRIX_FLOAT

	printf("x_r =\n");
	MatrixOutput(cout,x_r);

	printf("x_n =\n");
	MatrixOutput(cout,x_n);
}

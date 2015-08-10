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
Convert a matrix to row echelon form.
Ax = b -> Rx = c
b is optional and may contain multiple columns.
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
	printf("Input matrix A (empty line to finish):\n");
	MatrixInput(cin,a);
	Canonicalize(a);

	mx b;
	printf("Input vector (or matrix) b: \n");
	MatrixInput(cin,b);
	Canonicalize(b);

	//reduce A to R form
	mx r,c;
	if(!b.size1())
		RRef(a,r);
	else
		RRef(a,b,r,c);

#if defined(MATRIX_FLOAT)
	Canonicalize(r);
#endif  //MATRIX_FLOAT

	//Print out the resulting Rx = c system
	printf("Rank of A = %zu\n",Rank(r));

	if(!b.size1()){
		printf("R = \n");
		MatrixOutput(cout,r);
	}else{
#if defined(MATRIX_FLOAT)
		Canonicalize(c);
#endif  //MATRIX_FLOAT
		printf("Reducing the system Ax = b to Rx = c\n");
		printf("R = \n");
		MatrixOutput(cout,r);
		printf("c = \n");
		MatrixOutput(cout,c);
	}
}

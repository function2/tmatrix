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
Compute the determinant of a matrix.
Return error if the matrix is not invertible.
 */
#include<iostream>
#include<cstdlib>

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

	if(!IsSquare(a)){
		fprintf(stderr,"Invalid matrix: The matrix must be square.\n");
		return EXIT_FAILURE;
	}

	mx::value_type det = Determinant(a);
	cout << det << '\n';
	return det == mx::value_type(0) ? EXIT_FAILURE : EXIT_SUCCESS;
}

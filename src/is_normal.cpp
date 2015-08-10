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
Return true if a matrix commutes with N^\herm
N^\herm N = N N^\herm
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

	printf("a = ");
	MatrixOutput(cout,a);

	if(!IsSquare(a)){
		printf("Matrix must be square.\n");
		return EXIT_FAILURE;
	}

	if(MatrixEquality(prod(a,herm(a)), prod(herm(a),a))) {
		printf("Is normal.\n");
		return EXIT_SUCCESS;
	}else{
		printf("Is NOT normal.\n");
		return EXIT_FAILURE;
	}
}

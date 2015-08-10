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
Find the inverse and left/right inverse of a matrix.
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

	if(IsSquare(a)){
		mx inv;
		if(!Inverse(a,inv))
			printf("Inverse does not exist.\n");
		else{
#if defined(MATRIX_FLOAT)
			Canonicalize(inv);
#endif  //MATRIX_FLOAT
			printf("Inverse is ");
			MatrixOutput(cout,inv);
			return EXIT_SUCCESS;
		}
	}else
		printf("Rectangular matrix does not have an inverse.\n");

	// Check for left,right inverses.
	mx ata,aat;
	mx atai,aati;
	ata = prod(trans(a),a);
	aat = prod(a,trans(a));

	if(Inverse(ata,atai)){
		mx leftinv = prod(atai,trans(a));
#if defined(MATRIX_FLOAT)
		Canonicalize(leftinv);
#endif  //MATRIX_FLOAT
		printf("Left inverse is ");
		MatrixOutput(cout,leftinv);
	}else
		printf("Left inverse does not exist.\n");

	if(Inverse(aat,aati)){
		mx rightinv = prod(trans(a),aati);
#if defined(MATRIX_FLOAT)
		Canonicalize(rightinv);
#endif  //MATRIX_FLOAT
		printf("Right inverse is ");
		MatrixOutput(cout,rightinv);
	}else
		printf("Right inverse does not exist.\n");
}

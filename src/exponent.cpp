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
Find the exponent of a matrix e^A (using n terms of the power series)
 */
#include<cerrno>
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

	// Default number of terms (if not given on the command line).
	long n = 50;

	if(argc == 2) {
		char *endptr;
		errno = 0;
		n = strtol(argv[1],&endptr,0);
		if(*endptr != '\0' || errno != 0 || n < 1) return EXIT_FAILURE;
	}

	mx a;
	printf("Input matrix A (empty line to finish):\n");
	MatrixInput(cin,a);
	Canonicalize(a);

	printf("a = ");
	MatrixOutput(cout,a);

	mx b;
	Exp(a,n,b);
	printf("exp = ");
	MatrixOutput(cout,b);
}

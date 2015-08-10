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
Take a square matrix to a power n.
 */

#include<iostream>
#include<cerrno>

#include"las/matrix.h"
#include"las/mtypes.h"

using namespace std;
using namespace las;

int main(int argc, char**argv)
{
#if defined(MATRIX_FLOAT)
	mpf_set_default_prec(MATRIX_FLOAT_PREC);
#endif  //MATRIX_FLOAT

	if(argc != 2)  return EXIT_FAILURE;
	char *endptr;
	errno = 0;
	long n = strtol(argv[1],&endptr,0);
	if(*endptr != '\0' || errno != 0) return EXIT_FAILURE;

	mx a;
	printf("Input matrix A (empty line to finish):\n");
	MatrixInput(cin,a);
	Canonicalize(a);

	printf("a = ");
	MatrixOutput(cout,a);

	if(!IsSquare(a)){
		printf("A is not square!\n");
		return EXIT_FAILURE;
	}

	mx r;
	printf("result = ");
	MatrixOutput(cout,Pow(a,n,r));
}

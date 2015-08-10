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
PA = LDU factorization
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
	printf("Input matrix A (empty line to finish):\n");
	MatrixInput(cin,a);
	Canonicalize(a);

	if(!IsSquare(a)){
		fprintf(stderr,"Invalid matrix: The matrix must be square.\n");
		return EXIT_FAILURE;
	}

	mx p,l,d,u;
	mx::value_type det = Determinant(a,&p,&l,&d,&u);

	if(det == mx::value_type(0))
		printf("A is not invertible, factorization failed.\n");

	mx pa = prod(p,a);
	mx ld = prod(l,d);
	mx ldu = prod(ld,u);

#if defined(MATRIX_FLOAT)
	Canonicalize(l);
	Canonicalize(u);
	Canonicalize(ldu);
	// no need to canonicalize 'd', it shouldn't have zeros on the diagonal.
#endif  //MATRIX_FLOAT

	printf("p = "); MatrixOutput(cout,p);
	printf("pa = "); MatrixOutput(cout,pa);
	printf("l = "); MatrixOutput(cout,l);
	printf("d = "); MatrixOutput(cout,d);
	printf("u = "); MatrixOutput(cout,u);
	printf("l*d*u = "); MatrixOutput(cout,ldu);

	cout << format("Determinant = %1%\n") % det;

	// Check the result for accuracy.
	if(MatrixEquality(prod(trans(p),ldu),a)) cout << "CHECK GOOD\n";
	else {
		cerr << "CHECK ERROR\n";
		return EXIT_FAILURE;
	}
}

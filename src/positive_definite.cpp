// Copyright (C) 2014  Michael Seyfert <michael@codesand.org>
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
Determine whether the input matrix is positive definite.
A positive definite matrix has all of the following.
1 x^T A x > 0 for all nonzero real vectors x.
2 All the eigenvalues of A satisfy \lambda_i > 0.
3 All the upper left submatrices A_k have positive determinants.
4 All the pivots (without row exchanges) satisfy d_k > 0.

Any one of these is sufficient, they are either all true or none
are true.
 */

#include <iostream>
#include <boost/format.hpp>

#include "las/matrix.h"
#include "las/mtypes.h"

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

	// The determinant test is the easiest.
	if(!IsSymmetric(a)){
		fprintf(stderr,"Invalid matrix: the matrix is not symmetric (not positive definite).\n");
		return EXIT_FAILURE;
	}

	mx p,l,d,u;
	mx::value_type det = Determinant(a,&p,&l,&d,&u);
	Canonicalize(d);

	// Check all the pivots.
	for(size_t pivot = 0; pivot < d.size1();++pivot){
		mx::value_type pivot_val = d(pivot,pivot);
		cout << format("Pivot %1% = %2%\n") % (pivot+1) % pivot_val;
		if(pivot_val.real() <= mx::value_type(0).real()){
			printf("Is not positive definite.\n");
			return EXIT_FAILURE;
		}
	}

	printf("Is positive definite.\n");
	return EXIT_SUCCESS;
}

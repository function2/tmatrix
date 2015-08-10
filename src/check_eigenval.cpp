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
Given a matrix and some probable eigenvalues, check that they
are actually eigenvalues and print the eigenvectors.
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
	printf("Input matrix A:\n");
	MatrixInput(cin,a);
	Canonicalize(a);

	if(!IsSquare(a)) return 1;

	string line;
	while(printf("Enter probable eigenvalues:"), getline(cin,line)){
		istringstream iss(line);
		vector<mx::value_type> check;
		mx::value_type t;
		while(iss >> t) check.push_back(t);

		for(size_t k = 0;k < check.size();++k){
			Canonicalize(check[k]); // Reduce fractions.
			const mx::value_type &e = check[k];

			// test = A - \lambda I
			mx test = a;
			test -= e * IdentityMatrix<mx::value_type>(a.size1());

			mx r;
			RRef(test,r); // Nullspace() expects A to be in RRef form.

#if defined(MATRIX_FLOAT)
			Canonicalize(r);
#endif  //MATRIX_FLOAT

			// The nullspace should contain the eigenvectors.
			mx n;
			Nullspace(r,n);
			if(n.size2() == 0){
				cout << format("%1% is not an eigenvalue.\n") % e;
			}else{
				if(n.size2() == 1)
					cout << format("Eigenvector for eigenvalue %1%: ") % e;
				else
					cout << format("Eigenvectors for eigenvalue %1%: ") % e;
#if defined(MATRIX_FLOAT)
				Canonicalize(n);
#endif  //MATRIX_FLOAT
				MatrixOutput(cout,n);
			}
		}
	}
}

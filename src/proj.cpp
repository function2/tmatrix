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
Project b onto the column space of A.

If b is not given find a projection matrix P that will project any b. Pb.

(b-Av) is orthogonal to A   ->   AT(b - Av) = 0   ->   ATAv = ATb
Solve for v.
Av is a combination of the columns of A that is the projection of b
onto A.
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
	printf("Input matrix A\n");
	MatrixInput(cin,a);
	Canonicalize(a);
	printf("a =\n");
	MatrixOutput(cout,a);

	mx b;
	printf("Enter b to be projected onto A:\n");
	MatrixInput(cin,b);
	Canonicalize(b);

	mx at,ata;
	at = trans(a);
	ata = prod(at,a);

	if(b.size1()){
		if(a.size1() != b.size1()) return EXIT_FAILURE;

		mx atb;
		atb = prod(at,b);

		mx v;
		if(!SolvePivots(ata,atb,v)){
			printf("No solution?? ERROR.\n");
			return EXIT_FAILURE;
		}

		mx proj(prod(a,v));
		mx error(b - proj);

#if defined(MATRIX_FLOAT)
		Canonicalize(v);
		Canonicalize(proj);
		Canonicalize(error);
#endif  //MATRIX_FLOAT

		printf("v =\n"); MatrixOutput(cout,v);
		printf("Projection = Av =\n"); MatrixOutput(cout,proj);
		printf("Error = (b - Av) = \n"); MatrixOutput(cout,error);
		ostringstream oss;
		oss << prod(trans(error),error)(0,0);
		printf("E^2 = %s\n", oss.str().c_str());
	}

	//find the projection matrix (only works if the columns are independent)
	mx atai;
	if(Inverse(ata,atai)){
//		mx p = a * atai * at;
		mx p = prod(a,atai);
		p = prod(p,at);
		printf("Projection matrix =\n");
		MatrixOutput(cout,p);
	}
}

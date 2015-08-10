/*
Sandbox. Use this file to test code using the library!
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
	MatrixInput(cin,a); Canonicalize(a);

	mx b;
	printf("Input matrix B (empty line to finish):\n");
	MatrixInput(cin,b); Canonicalize(b);

	mx c(4,4);
	c(0,0) = a(0,0) - b(0,0);
	c(0,1) = -b(1,0);
	c(0,2) = a(0,1);
	c(0,3) = 0;

	c(1,0) = a(1,0);
	c(1,1) = 0;
	c(1,2) = a(1,1) - b(0,0);
	c(1,3) = -b(1,0);

	c(2,0) = -b(0,1);
	c(2,1) = a(0,0) - b(1,1);
	c(2,2) = 0;
	c(2,3) = a(0,1);

	c(3,0) = 0;
	c(3,1) = a(1,0);
	c(3,2) = -b(0,1);
	c(3,3) = a(1,1) - b(1,1);

	mx r;
	RRef(c,r);
	mx n;
	Nullspace(r,n);

	mx t(n.size2(),1);
	for(int k = 0;k < (int)n.size2();++k)
		t(k,0) = rand() % 2 + 1;

	mx v = prod(n,t);
	mx m(2,2);
	m(0,0) = v(0,0);
	m(0,1) = v(1,0);
	m(1,0) = v(2,0);
	m(1,1) = v(3,0);
	cout << "M = "; MatrixOutput(cout,m);
	mx mi;
	if(!Inverse(m,mi)) throw int();
	MatrixOutput(cout,mi);

	t = prod(mi,a);
	t = prod(t,m);
	MatrixOutput(cout,t);
	cout << "Equal? " << MatrixEquality(b,t) << endl;
}

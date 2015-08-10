// Copyright (C) 2012  Michael Seyfert <m@codesand.org>
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
#ifndef LAS_MATRIX_INL_H_
#define LAS_MATRIX_INL_H_

#include<vector>

namespace las{

// Create a permutation matrix from a vector permutation
template<class M1>
static void PermutationMatrix(M1 &p,const std::vector<size_t> &pt)
{
	p.resize(pt.size(),pt.size(),false);
	for(size_t k = 0;k < pt.size();++k)
		p(k,pt[k]) = (typename M1::value_type)(1);
}

template<class M1>
static void GetSubMatrixRows(
	const M1 &m,const std::vector<size_t> &rows,M1 &res)
{
	if(!rows.size() || !m.size1()){
		res.clear();
		return;
	}
	res.resize(rows.size(),m.size2(),false);
	for(size_t k = 0;k < rows.size();++k)
		for(size_t j = 0;j < m.size2();++j)
			res(k,j) = m(rows[k],j);
}

template<class M1>
static void GetSubMatrixRows(
	const M1 &m,int low,int high,M1 &res)
{
	if(low > high || !m.size1()){
		res.clear();
		return;
	}

	res.resize(high-low+1,m.size2(),false);
	for(size_t k = low;k <= (size_t)high;++k)
		for(size_t j = 0;j < m.size2();++j)
			res(k-low,j) = m(k,j);
}

template<class M1>
static void GetSubMatrixCols(
	const M1 &m,const std::vector<size_t> &cols,M1 &res)
{
	if(!cols.size() || !m.size1()){
		res.clear();
		return;
	}
	res.resize(m.size1(),cols.size(),false);
	for(size_t k = 0;k < m.size1();++k)
		for(size_t j = 0;j < cols.size();++j)
			res(k,j) = m(k,cols[j]);
}

template<class M1,class M2>
M2& Pow(const M1 &a,long n,M2 &r)
{
	r = IdentityMatrix<typename M2::value_type>(a.size1());

	M1 p = a;
	while(true){
		if(n & 1) r = prod(r,p);
		n >>= 1;
		if(!n) break;
		p = prod(p,p);
	}

	return r;
}

template<class M1,class M2>
M2& Exp(const M1 &a,long terms,M2 &r)
{
	r = IdentityMatrix<typename M2::value_type>(a.size1());

	M2 power = IdentityMatrix<typename M2::value_type>(a.size1());
	typename M2::value_type fact = typename M2::value_type(1); // factorial

	for(long k = 1;k < terms; ++k) {
		fact *= typename M2::value_type(k);
		power = prod(power,a);
		r += power / fact;
	}

	return r;
}

template<class M1>
int NumRowExchanges(const M1 &p)
{
	// Determine how many row swaps there are in this permutation matrix.
	int num_row_swaps = 0;
	std::vector<bool> checked(p.size1(),false);
	for(size_t k = 0;k < p.size1();++k){
		int cur = k;
		int this_loop = 0;
		// Each loop found can subtract one from the number of rows out of order.
		while(!checked[cur]){
			checked[cur] = true;
			if(p(cur,cur) != typename M1::value_type(1)){
				this_loop++;
				for(size_t j = 0;j < p.size2();++j)
					if(p(cur,j) != typename M1::value_type(0)){
						cur = j;
						break;
					}
			}
		}
		if(this_loop) num_row_swaps += this_loop - 1;
	}
	return num_row_swaps;
}

template<class M1>
static bool IsZeroRow(const M1 &m,size_t row)
{
	for(size_t k = 0;k < m.size2();++k)
		if(m(row,k) != (typename M1::value_type)(0)) return false;
	return true;
}


template<class M1,class M2>
bool MatrixEquality(const M1 &m1,const M2 &m2)
{
	if(m1.size1() != m2.size1() || m1.size2() != m2.size2())
		return false;

	for(size_t row = 0;row < m1.size1();++row)
		for(size_t col = 0;col < m1.size2();++col)
			if(m1(row,col) != m2(row,col))
				return false;

	return true;
}

template<class M1>
bool IsSquare(const M1 &m)
{
	return m.size1() == m.size2();
}

template<class M1>
bool IsSymmetric(const M1 &a)
{
	if(!IsSquare(a)) return false;

	for(size_t row = 0; row < a.size1(); ++row)
		for(size_t col = row+1; col < a.size2(); ++col)
			if(a(row,col) != a(col,row))
				return false;

	return true;
}

template<class M1>
void SwapRows(M1 &m,size_t r1,size_t r2)
{
	if(r1 == r2) return;

	for(size_t k = 0;k < m.size2();++k)
		std::swap(m(r1,k),m(r2,k));
}

template<class M1>
void SwapCols(M1 &m,size_t c1,size_t c2)
{
	if(c1 == c2) return;

	for(size_t k = 0;k < m.size1();++k)
		std::swap(m(k,c1),m(k,c2));
}

template<class M1,class T>
void DivRow(M1 &m,size_t row,const T &val)
{
	for(size_t k = 0;k < m.size2();++k)
		m(row,k) /= val;
}

template<class M1,class T>
void DivCol(M1 &m,size_t col,const T &val)
{
	for(size_t k = 0;k < m.size1();++k)
		m(k,col) /= val;
}

template<class M1,class T>
void MultCol(M1 &m,size_t col,const T &val)
{
	for(size_t k = 0;k < m.size1();++k)
		m(k,col) *= val;
}

template<class M1,class T>
void MultSubRow(M1 &m,size_t src,size_t dest,const T &val)
{
	for(size_t k = 0;k < m.size2();++k)
		m(dest,k) -= m(src,k) * val;
}

template<class M1,class T>
void MultSubCol(M1 &m,size_t src,size_t dest,const T &val)
{
	for(size_t k = 0;k < m.size1();++k)
		m(k,dest) -= m(k,src) * val;
}

template<class M1,class M2,class M3,class M4>
void RRef(const M1 &a,const M2 &b, M3 &r, M4 &c)
{
	//Ax = b   ->   Rx = c
	if(b.size1() < a.size1()) throw int(); //TODO
	r = a;
	c = b;

	typedef typename M1::value_type T;

	//triangular form, or the closest thing we can get to that form.
	size_t curr_row = 0;
	for(size_t curr_col = 0; curr_col < a.size2() ; ++curr_col){
		size_t k;
		for(k = curr_row;k < a.size1();++k)
			if(r(k,curr_col) != T()) break;

		if(k == a.size1()) continue; //couldn't find nonzero pivot

		SwapRows(r,k,curr_row); //make this row our current
		SwapRows(c,k,curr_row);

		//make this pivot 1
		T tmp = r(curr_row,curr_col);
		DivRow(r,curr_row,tmp);
		DivRow(c,curr_row,tmp);

		//zero out everything below this pivot, and above it
		for(k = 0;k < a.size1();++k)
			if(k != curr_row){
				if(r(k,curr_col) == T()) continue;

				tmp = r(k,curr_col);
				MultSubRow(r,curr_row,k,tmp);
				MultSubRow(c,curr_row,k,tmp);
			}

		++curr_row;
	}
}

template<class M1,class M2>
void RRef(const M1 &a,M2 &r)
{
	M1 b(a.size1(),1),c;
	RRef(a,b,r,c);
}

template<class M1>
size_t Rank(const M1 &a)
{
	//Simply count the number of nonzero rows.
	// This works because of the input restriction.
	size_t c=0;
	for(size_t k = 0;k < a.size1();++k)
		if(!IsZeroRow(a,k)) ++c;
	return c;
}

template<class M1>
void GetNonZeroRows(const M1 &a,M1 &out)
{
	size_t num_nonzero_rows = Rank(a);
	if(num_nonzero_rows == 0){
		out.clear();
		return;
	}

	out.resize(num_nonzero_rows,a.size2());
	size_t curr_out_row = 0;
	for(size_t row = 0;row < a.size1();++row){
		if(!IsZeroRow(a,row)){
			for(size_t col = 0;col < a.size2();++col)
				out(curr_out_row,col) = a(row,col);
			++curr_out_row;
		}
	}
}

template<class M1>
void Nullspace(const M1 &a,M1 &n)
{
	std::vector<size_t> tmp;
	Nullspace(a,n,tmp);
}

template<class M1>
void Nullspace(const M1 &a,M1 &n,std::vector<size_t> &pivot_cols)
{
	Nullspace(a,n,Rank(a),pivot_cols);
}

template<class M1>
void Nullspace(const M1 &a,M1 &n,size_t rank,std::vector<size_t> &pivot_cols)
{
	typedef typename M1::value_type T;

	if(a.size2() == rank)
		n.clear();
	else
		n.resize(a.size2(),a.size2() - rank,false);

	size_t curr_row=0;
	size_t curr_nullspace_col = 0;
	for(size_t curr_col = 0; curr_col < a.size2(); ++curr_col){

		//found a free column?
		if(curr_row >= a.size1() || a(curr_row,curr_col) == T(0)){

			//put the values in this nullspace column
			n(curr_col, curr_nullspace_col) = T(1);
			for(size_t k = 0;k < curr_row;++k)
				n(pivot_cols[k], curr_nullspace_col) = -a(k,curr_col);

			++curr_nullspace_col;

		}else{
			++curr_row;
			pivot_cols.push_back(curr_col);
		}
	}
}

template<class M1>
void GetPivotCols(const M1 &a,std::vector<size_t> &pivot_cols)
{
	typedef typename M1::value_type T;
	std::vector<size_t> t;

	size_t curr_row=0;
	for(size_t curr_col = 0; curr_col < a.size2(); ++curr_col){
		//found a free column?
		if(curr_row >= a.size1() || a(curr_row,curr_col) == T(0));
		else{
			pivot_cols.push_back(curr_col);
			++curr_row;
		}
	}
}

template<class M1>
bool Inverse(const M1 &a,M1 &inv)
{
	if(!IsSquare(a)) return false;
	M1 I(boost::numeric::ublas::identity_matrix<typename M1::value_type>(a.size1()));
	M1 b(I),r;
	RRef(a,b,r,inv);

	//if R is the identity, it has an inverse
	return MatrixEquality(r,I);
}

template<class M1>
void LU(const M1 &a,M1 &p,M1 &l,M1 &u)
{
	typedef typename M1::value_type T;

	//initialize l,u, and pt (permutation vector)
	l = boost::numeric::ublas::identity_matrix<T>(a.size1());
	u = a;
	std::vector<size_t> pt(a.size1());
	for(size_t i = 0;i < a.size1();++i) pt[i] = i;

	//triangular form
	size_t curr_row = 0;
	for(size_t curr_col = 0; curr_col < a.size2() ; ++curr_col){
		size_t k;
		for(k = curr_row;k < a.size1();++k)
			if(u(k,curr_col) != T(0)) break;

		if(k == a.size1()) continue; //couldn't find nonzero pivot

		SwapRows(u,k,curr_row); //make this row our current
		SwapCols(l,k,curr_row);
		SwapRows(l,k,curr_row);
		std::swap(pt[k],pt[curr_row]);

		//make this pivot 1
		T tmp = u(curr_row,curr_col);
		DivRow(u,curr_row,tmp);
		MultCol(l,curr_row,tmp);

		//zero out everything below this pivot
		for(k = curr_row+1;k < a.size1();++k){
			if(u(k,curr_col) == T(0)) continue;

			tmp = u(k,curr_col);
			MultSubRow(u,curr_row,k,tmp);
			MultSubCol(l,k,curr_row,-tmp);
		}

		++curr_row;
	}

	//Create permutation matrix
	PermutationMatrix(p,pt);
}

template<class M1>
bool LDU(const M1 &a,M1 &p,M1 &l,M1 &d,M1 &u)
{
	LU(a,p,l,u);

	d = boost::numeric::ublas::identity_matrix<typename M1::value_type>(a.size1());

	bool good = true;
	for(size_t k = 0;k < a.size1();++k){
		if(u(k,k) == typename M1::value_type(0))
			good = false;

		typename M1::value_type pivot = l(k,k);
		if(pivot != (typename M1::value_type)(1))
			DivCol(l,k,pivot);
		d(k,k) = pivot;
	}

	return good;
}

template<class T>
std::istream& MatrixInput(
	std::istream &in,
	boost::numeric::ublas::matrix<T> &m)
{
	std::string t;
	std::vector<std::string> v;
	while(getline(in,t) && t != "")
		v.push_back(t);

	if(v.empty()) return in;

	//count columns
	T tmp;
	size_t cols = 0;
	std::istringstream iss(v[0]);
	while(iss >> tmp) ++cols;

	m.resize(v.size(),cols,false);

	for(size_t row = 0;row < v.size();++row){
		std::istringstream iss(v[row]);
		size_t col = 0;
		while(iss >> tmp)
			m(row,col++) = tmp;

		if(col != cols) throw int(); //TODO: throw a real exception.
	}

	return in;
}

template<class M1>
std::ostream& MatrixOutput(
	std::ostream &os,
	const M1 &m)
{
	os << "[" << m.size1() << ',' << m.size2() << "]:\n";
	for(size_t row = 0;row < m.size1();++row){
		for(size_t col = 0;col < m.size2();++col)
			os << ' ' << m(row,col);
		os << '\n';
	}
	return os;
}

template<class M1>
bool SolvePivots(const M1 &a,const M1 &b,M1 &x)
{
	//row echelon reduce first
	M1 r,c;
	RRef(a,b,r,c);

	return SolvePivotsR(r,c,x);
}

template<class M1>
bool SolvePivotsR(const M1 &r,const M1 &c,M1 &x)
{
	//for every row(r) == 0, we must have c == 0 to solve
	for(size_t k = 0;k < r.size1();++k)
		if(IsZeroRow(r,k) && !IsZeroRow(c,k)) return false;

	//find pivot columns
	std::vector<size_t> pivot_cols;
	GetPivotCols(r,pivot_cols);

	x.resize(r.size2(),c.size2(),false);
	//every pivot column has exactly one pivot row associated with it.
	// pivot rows are just 0,1,2,...,rank-1
	for(size_t j = 0;j < c.size2();++j){
		for(size_t k = 0;k < pivot_cols.size();++k){
			x(pivot_cols[k],j) = c(k,j);
		}
	}

	return true;
}

namespace {
// Template to switch between hermitian inner product and regular.
template<class M1,class M2,class M3,bool use_herm>
bool UQR_Factor(const M1 &a,M2 &q,M3 &r)
{
	typedef typename M2::value_type val;
	q = a;

	for(size_t k = 0;k < q.size2();++k){
		for(size_t j = 0;j < k;++j){
			val length2 = use_herm ? inner_prod(herm(column(q,j)),column(q,j)) : inner_prod(column(q,j),column(q,j));
			if(length2 == val(0)) continue;
			MultSubCol(q,j,k,(use_herm ? inner_prod(herm(column(q,j)),column(q,k)) : inner_prod(column(q,j),column(q,k))) / length2);
		}
	}

	bool ret = true;

	r.resize(q.size2(),q.size2(),false);
	for(size_t k = 0;k < r.size2();++k){
		for(size_t j = 0;j <= k;++j){
			if(k == j) r(j,k) = 1;
			val length2 = (use_herm ? inner_prod(herm(column(q,j)),column(q,j)) : inner_prod(column(q,j),column(q,j)) );
			if(length2 == val(0)){
				ret = false;
				continue; //continue, even though we will return false
			}else
				r(j,k) = (use_herm ? inner_prod(herm(column(q,j)),column(a,k)) : inner_prod(column(q,j),column(a,k)) ) / length2;
		}
	}

	return ret;
}

//Template to switch between hermitian inner product and regular.
template<class M1,class M2,class M3,bool use_herm>
bool UQR_FactorNormal(const M1 &a,M2 &q,M3 &r)
{
	typedef typename M2::value_type val;
	q = a;

	bool ret = true;

	for(size_t k = 0;k < q.size2();++k){

		for(size_t j = 0;j < k;++j)
			MultSubCol(q,j,k, (use_herm ? inner_prod(herm(column(q,j)),column(q,k)) : inner_prod(column(q,j),column(q,k))) );

		// Normalize this column.
		val length = (use_herm ? inner_prod(herm(column(q,k)),column(q,k)) : inner_prod(column(q,k),column(q,k)));
		if(length == val(0)){
			ret = false;
			continue; //continue, even though we will return false
		}
		length = sqrt(length);
		DivCol(q,k,length);
	}

	r.resize(q.size2(),q.size2(),false);
	for(size_t k = 0;k < r.size2();++k)
		for(size_t j = 0;j <= k;++j)
			r(j,k) = (use_herm ? inner_prod(herm(column(q,j)),column(a,k)) : inner_prod(column(q,j),column(a,k)));

	return ret;
}
} // namespace <unnamed>

template<class M1,class M2,class M3>
bool QR_Factor(const M1 &a,M2 &q,M3 &r)
{
	return UQR_Factor<M1,M2,M3,false>(a,q,r);
}

template<class M1,class M2,class M3>
bool UR_Factor(const M1 &a,M2 &q,M3 &r)
{
	return UQR_Factor<M1,M2,M3,true>(a,q,r);
}

template<class M1,class M2,class M3>
bool QR_FactorNormal(const M1 &a,M2 &q,M3 &r)
{
	return UQR_FactorNormal<M1,M2,M3,false>(a,q,r);
}

template<class M1,class M2,class M3>
bool UR_FactorNormal(const M1 &a,M2 &q,M3 &r)
{
	return UQR_FactorNormal<M1,M2,M3,true>(a,q,r);
}

template<class M1>
typename M1::value_type Determinant(const M1 &a,M1 *p,M1 *l,M1 *d,M1 *u)
{
	// local scope (if the caller doesn't want these returned)
	M1 pm,lm,dm,um;
	if(p == NULL) p = &pm;
	if(l == NULL) l = &lm;
	if(d == NULL) d = &dm;
	if(u == NULL) u = &um;

	typename M1::value_type det(1);

	if(!LDU(a,*p,*l,*d,*u)){
		det = typename M1::value_type(0);
	}else{
		// Compute the determinant.
		for(size_t k = 0;k < d->size1();++k)
			det *= (*d)(k,k); // Product of pivots.
		if(NumRowExchanges(*p) % 2) //FIXME this is slow
			det *= typename M1::value_type(-1); // Sign.
	}

	return det;
}

template<class M1,class M2>
void MinorMatrix(const M1 &a,M2 &m,size_t row,size_t col)
{
	size_t row_low[] = {0,row+1};
	size_t row_hi[]  = {row,a.size1()};
	size_t col_low[] = {0,col+1};
	size_t col_hi[]  = {col,a.size2()};
	size_t add_val[] = {0,(size_t)(-1)};
	static const size_t M = sizeof(row_low) / sizeof(row_low[0]);

	m.resize(a.size1()-1,a.size2()-1,false);

	for(size_t kr = 0;kr < M;++kr){
		size_t row_add = add_val[kr];
		for(size_t kc = 0;kc < M;++kc){
			size_t col_add = add_val[kc];

			for(size_t r = row_low[kr]; r < row_hi[kr];++r)
				for(size_t c = col_low[kc]; c < col_hi[kc];++c)
					m(r + row_add,c + col_add) = a(r,c);

		}
	}

}

template<class M1,class M2>
void GetCofactorMatrix(const M1 &a,M2 &cofactor)
{
	// FIXME find a way to speed this up.
	typedef typename M1::value_type vtype;

	M1 inverse;
	if(Inverse(a,inverse)){
		//matrix has an inverse.
		cofactor = Determinant(a) * trans(inverse);
		return;
	}

	//matrix does not have an inverse.
	cofactor.resize(a.size1(),a.size2(),false);
	M1 m; //minor

	for(size_t r = 0;r < a.size1();++r){
		for(size_t c = 0;c < a.size2();++c){
			MinorMatrix(a,m,r,c);
			cofactor(r,c) = (((r+c)%2) ? vtype(-1) : vtype(1)) * Determinant(m);
		}
	}
}

} //namespace las

#endif // LAS_MATRIX_INL_H_

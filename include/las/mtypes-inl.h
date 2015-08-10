// Copyright (C) 2008  Michael Seyfert <m@codesand.org>
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

#ifndef MTYPES_INL_H_
#define MTYPES_INL_H_
namespace las{

template<class T>
void Canonicalize(boost::numeric::ublas::matrix<T> &m)
{
	for(size_t row = 0;row < m.size1();++row)
		for(size_t col = 0;col < m.size2();++col)
			Canonicalize(m(row,col));
}

template<class T>
void Canonicalize(boost::numeric::ublas::triangular_matrix<T,boost::numeric::ublas::upper> &m)
{
	for(size_t row = 0;row < m.size1();++row)
		for(size_t col = row;col < m.size2();++col)
			Canonicalize(m(row,col));
}

} //namespace las
#endif  //MTYPES_INL_H_

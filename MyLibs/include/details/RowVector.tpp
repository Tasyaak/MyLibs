#include "RowVector.hpp"
#include "ColumnVector.hpp"
#include "Matrix.hpp"

template <typename T>
T RowVector<T>::operator * (const ColumnVector<T>& v) const
{
	assert(n == v.size() && "RowVector size doesn't match ColumnVector size");

	std::size_t size = n;
	T s = T{};
	for (std::size_t i = 0; i < size; ++i)
		s += data[i] * v.data[i];
	return s;
}
template <typename T>
RowVector<T> RowVector<T>::operator * (const Matrix<T>& A) const
{
	assert(n == A.numRows() && "RowVector size doesn't match Matrix number of rows");

	std::size_t size = n, m = A.numCols();
	RowVector<T> u(m);
	for (std::size_t i = 0; i < m; ++i)
		for (std::size_t j = 0; j < n; ++j)
			u[i] += data[j] * A(j, i);
	return u;
}
template <typename T>
ColumnVector<T> RowVector<T>::operator ~ () const { return ColumnVector<T>(data, n); }
#include "ColumnVector.hpp"
#include "RowVector.hpp"
#include "Matrix.hpp"

template <typename T>
Matrix<T> ColumnVector<T>::operator * (const RowVector<T>& v) const
{
	assert(n == v.size() && "ColumnVector size doesn't match RowVector size");

	std::size_t size = n, m = v.size();
	Matrix<T> A(size, m, mathDetails::NoInitTag{});
	for (std::size_t i = 0; i < size; ++i)
		for (std::size_t j = 0; j < m; ++j)
			A(i, j) = data[i] * v[j];
	return A;
}
template <typename T>
RowVector<T> ColumnVector<T>::operator ~ () const { return RowVector<T>(data, n); }
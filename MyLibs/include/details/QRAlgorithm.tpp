#include <cmath>
#include <tuple>
#include "QRAlgorithm.hpp"
#include "GaussMethod.hpp"
#include "Matrix.hpp"
#include "ColumnVector.hpp"

template<typename T>
T* QRAlgorithm<T>::solveSystem(const Matrix<T>& A, T* b, const double eps)
{
	std::size_t n = A.numRows();
	Matrix<T> Q, R;
	std::tie(Q, R) = A.QRdecomposition(eps);
	T* qb = new T[n]{};
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j)
			qb[i] += Q(i, j) * b[j];

	return GaussMethod<T>::backSubstitution(R, qb, eps);
}

template<typename T>
ColumnVector<T> QRAlgorithm<T>::solveSystem(const Matrix<T>& A, const ColumnVector<T>& b, const double eps)
{
	std::size_t n = A.numRows();
	assert(n == b.size() && "System matrix order doesn't match ColumnVector size");
	Matrix<T> Q, R;
	std::tie(Q, R) = A.QRdecomposition(eps );
	ColumnVector<T> B = Q * b;

	return GaussMethod<T>::backSubstitution(R, B, eps);
}
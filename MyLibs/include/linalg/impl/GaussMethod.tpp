#include <cmath>
#include "../GaussMethod.hpp"
#include "../Matrix.hpp"
#include "../ColumnVector.hpp"

template<typename T>
void GaussMethod<T>::forwardElimination(Matrix<T>& A, const double eps)
{
	std::size_t n = A.numRows();
	assert(n == A.numCols() && "System matrix isn't square");

	for (std::size_t i = 0; i < n; ++i)
	{
		if (abs(A(i, i) < eps))
		{
			short ind = -1;
			for (std::size_t j = i + 1; j < n; ++j)
				if (abs(A(j, i)) >= eps)
				{
					ind = j;
					break;
				}
			if (ind == -1)
				continue;
			A.swapRows(i, ind);
		}
		for (std::size_t j = i + 1; j < n; ++j)
			if (abs(A(j, i)) >= eps)
			{
				T k = A(j, i) / A(i, i);
				A(j, i) = T{};
				for (std::size_t l = i + 1; l < n; ++l)
					A(j, l) -= A(i, l) * k;
			}
	}
}

template<typename T>
void GaussMethod<T>::forwardElimination(Matrix<T>& A, T* b, const double eps)
{
	std::size_t n = A.numRows();
	assert(n == A.numCols() && "System matrix isn't square");

	for (std::size_t i = 0; i < n; ++i)
	{
		if (abs(A(i, i)) < eps)
		{
			short ind = -1;
			for (std::size_t j = i + 1; j < n; ++j)
				if (abs(A(j, i)) >= eps)
				{
					ind = j;
					break;
				}
			if (ind == -1)
				continue;
			A.swapRows(i, ind);
			std::swap(b[i], b[ind]);
		}
		for (std::size_t j = i + 1; j < n; ++j)
			if (abs(A(j, i)) >= eps)
			{
				T k = A(j, i) / A(i, i);
				A(j, i) = T{};
				for (std::size_t l = i + 1; l < n; ++l)
					A(j, l) -= A(i, l) * k;
				b[j] -= b[i] * k;
			}
	}
}

template<typename T>
void GaussMethod<T>::forwardElimination(Matrix<T>& A, ColumnVector<T>& b, const double eps)
{
	std::size_t n = A.numRows();
	assert(n == A.numCols() && "System matrix isn't square");
	assert(n == b.size() && "System matrix order doesn't match ColumnVector size");

	for (std::size_t i = 0; i < n; ++i)
	{
		if (abs(A(i, i)) < eps)
		{
			short ind = -1;
			for (std::size_t j = i + 1; j < n; ++j)
				if (abs(A(j, i)) >= eps)
				{
					ind = j;
					break;
				}
			if (ind == -1)
				continue;
			A.swapRows(i, ind);
			std::swap(b[i], b[ind]);
		}
		for (std::size_t j = i + 1; j < n; ++j)
			if (abs(A(j, i)) >= eps)
			{
				T k = A(j, i) / A(i, i);
				A(j, i) = T{};
				for (std::size_t l = i + 1; l < n; ++l)
					A(j, l) -= A(i, l) * k;
				b[j] -= b[i] * k;
			}
	}
}

template<typename T>
T* GaussMethod<T>::backSubstitution(const Matrix<T>& A, T* b, const double eps)
{
	std::size_t n = A.numRows();
	assert(n == A.numCols() && "System matrix isn't square");
	for (std::size_t i = 0; i < n; ++i)
		assert(std::abs(A(i, i)) >= eps && "SLAE has more than one solution");

	T* x = new T[n];
	for (std::size_t i = n; i-- > 0;)
	{
		T s = T{};
		for (std::size_t j = i + 1; j < n; ++j)
			s += A(i, j) * x[j];
		x[i] = (b[i] - s) / A(i, i);
	}
	return x;
}

template<typename T>
ColumnVector<T> GaussMethod<T>::backSubstitution(const Matrix<T>& A, const ColumnVector<T>& b, const double eps)
{
	std::size_t n = A.numRows();
	assert(n == A.numCols() && "System matrix isn't square");
	assert(n == b.size() && "System matrix order doesn't match ColumnVector size");
	for (std::size_t i = 0; i < n; ++i)
		assert(std::abs(A(i, i)) >= eps && "SLAE has more than one solution");

	ColumnVector<T> x(n, mathDetails::NoInitTag{});
	for (std::size_t i = n; i-- > 0;)
	{
		T s = T{};
		for (std::size_t j = i + 1; j < n; ++j)
			s += A(i, j) * x[j];
		x[i] = (b[i] - s) / A(i, i);
	}
	return x;
}

template<typename T>
T* GaussMethod<T>::solveSystem(Matrix<T> A, T* b, const double eps)
{
	std::size_t n = A.numRows();
	T* bt = new T[n];
	std::copy(b, b + n, bt);
	forwardElimination(A, bt, eps);

	T* res = backSubstitution(A, bt, eps);
	delete[] bt;

	return res;
}

template<typename T>
ColumnVector<T> GaussMethod<T>::solveSystem(Matrix<T> A, ColumnVector<T> b, const double eps)
{
	forwardElimination(A, b, eps);

	return backSubstitution(A, b, eps);
}
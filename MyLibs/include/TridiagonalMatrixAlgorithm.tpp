#include "TridiagonalMatrixAlgorithm.hpp"
#include "Matrix.hpp"
#include "ColumnVector.hpp"
#include "mathDetails.hpp"

template<typename T>
T* TridiagonalMatrixAlgorithm<T>::solveSystem(const Matrix<T>& A, T* b)
{
	std::size_t n = A.numRows();

	T* alp = new T[n - 1];
	T* bet = new T[n];
	alp[0] = -A(0, 1) / A(0, 0);
	bet[0] = b[0] / A(0, 0);

	for (std::size_t i = 1; i < n - 1; ++i)
	{
		T temp = A(i, i) + A(i, i - 1) * alp[i - 1];
		alp[i] = -A(i, i + 1) / temp;
		bet[i] = (b[i] - A(i, i - 1) * bet[i - 1]) / temp;
	}
	bet[n - 1] = (b[n - 1] - A(n - 1, n - 2) * bet[n - 2]) / (A(n - 1, n - 1) + A(n - 1, n - 2) * alp[n - 2]);

	T* x = new T[n];
	x[n - 1] = bet[n - 1];
	for (std::size_t i = n - 2; i >= 0; --i)
		x[i] = alp[i] * x[i + 1] + bet[i];
	delete[] alp;
	delete[] bet;

	return x;
}

template<typename T>
ColumnVector<T> TridiagonalMatrixAlgorithm<T>::solveSystem(const Matrix<T>& A, const ColumnVector<T>& b)
{
	std::size_t n = A.numRows();

	T* alp = new T[n - 1];
	T* bet = new T[n];
	alp[0] = -A(0, 1) / A(0, 0);
	bet[0] = b[0] / A(0, 0);

	for (std::size_t i = 1; i < n - 1; ++i)
	{
		T temp = A(i, i) + A(i, i - 1) * alp[i - 1];
		alp[i] = -A(i, i + 1) / temp;
		bet[i] = (b[i] - A(i, i - 1) * bet[i - 1]) / temp;
	}
	bet[n - 1] = (b[n - 1] - A(n - 1, n - 2) * bet[n - 2]) / (A(n - 1, n - 1) + A(n - 1, n - 2) * alp[n - 2]);

	ColumnVector<T> x(n, mathDetails::NoInitTag{})
	x[n - 1] = bet[n - 1];
	for (std::size_t i = n - 2; i >= 0; --i)
		x[i] = alp[i] * x[i + 1] + bet[i];
	delete[] alp;
	delete[] bet;

	return x;
}
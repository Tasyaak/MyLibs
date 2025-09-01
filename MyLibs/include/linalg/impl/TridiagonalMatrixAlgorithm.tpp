#include "../TridiagonalMatrixAlgorithm.hpp"
#include "../Matrix.hpp"
#include "../ColumnVector.hpp"

template<typename T>
T* TridiagonalMatrixAlgorithm<T>::solveSystem(const Matrix<T>& A, T* b)
{
	std::size_t n = A.numRows();

	T* alpha = new T[n - 1];
	T* beta = new T[n];
	alpha[0] = -A(0, 1) / A(0, 0);
	beta[0] = b[0] / A(0, 0);

	for (std::size_t i = 1; i < n - 1; ++i)
	{
		T temp = A(i, i) + A(i, i - 1) * alpha[i - 1];
		alpha[i] = -A(i, i + 1) / temp;
		beta[i] = (b[i] - A(i, i - 1) * beta[i - 1]) / temp;
	}
	beta[n - 1] = (b[n - 1] - A(n - 1, n - 2) * beta[n - 2]) / (A(n - 1, n - 1) + A(n - 1, n - 2) * alpha[n - 2]);

	T* x = new T[n];
	x[n - 1] = beta[n - 1];
	for (std::size_t i = n - 1; i-- > 0;)
		x[i] = alpha[i] * x[i + 1] + beta[i];
	delete[] alpha;
	delete[] beta;

	return x;
}

template<typename T>
ColumnVector<T> TridiagonalMatrixAlgorithm<T>::solveSystem(const Matrix<T>& A, const ColumnVector<T>& b)
{
	std::size_t n = A.numRows();

	T* alpha = new T[n - 1];
	T* beta = new T[n];
	alpha[0] = -A(0, 1) / A(0, 0);
	beta[0] = b[0] / A(0, 0);

	for (std::size_t i = 1; i < n - 1; ++i)
	{
		T temp = A(i, i) + A(i, i - 1) * alpha[i - 1];
		alpha[i] = -A(i, i + 1) / temp;
		beta[i] = (b[i] - A(i, i - 1) * beta[i - 1]) / temp;
	}
	beta[n - 1] = (b[n - 1] - A(n - 1, n - 2) * beta[n - 2]) / (A(n - 1, n - 1) + A(n - 1, n - 2) * alpha[n - 2]);

	ColumnVector<T> x(n, mathDetails::NoInitTag{});
	x[n - 1] = beta[n - 1];
	for (std::size_t i = n - 1; i-- > 0;)
		x[i] = alpha[i] * x[i + 1] + beta[i];
	delete[] alpha;
	delete[] beta;

	return x;
}
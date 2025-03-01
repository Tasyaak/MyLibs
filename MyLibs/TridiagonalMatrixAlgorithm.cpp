#include "TridiagonalMatrixAlgorithm.h"
#include "Matrix.h"

double* TridiagonalMatrixAlgorithm::solveSystem(const Matrix& A, const double* b)
{
	int n = A.getN();

	double* alp = new double[n - 1];
	double* bet = new double[n];
	alp[0] = -A[0][1] / A[0][0];
	bet[0] = b[0] / A[0][0];

	for (int i = 1; i < n - 1; ++i)
	{
		double temp = A[i][i] + A[i][i - 1] * alp[i - 1];
		alp[i] = -A[i][i + 1] / temp;
		bet[i] = (b[i] - A[i][i - 1] * bet[i - 1]) / temp;
	}
	bet[n - 1] = (b[n - 1] - A[n - 1][n - 2] * bet[n - 2]) / (A[n - 1][n - 1] + A[n - 1][n - 2] * alp[n - 2]);

	double* x = new double[n];
	x[n - 1] = bet[n - 1];
	for (int i = n - 2; i >= 0; --i)
		x[i] = alp[i] * x[i + 1] + bet[i];

	return x;
}
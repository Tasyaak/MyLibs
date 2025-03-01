#include <iostream>
#include <stdexcept>
#include "GaussMethod.h"
#include "Matrix.h"
#include "ColumnVector.h"

const double eps = 1e-10;

Matrix GaussMethod::forwardElimination(Matrix A)
{
	int n = A.getN();
	for (int i = 0; i < n; ++i)
	{
		if (abs(A[i][i]) < eps)
		{
			int ind = -1;
			for (int j = i + 1; j < n; ++j)
				if (abs(A[j][i]) >= eps)
				{
					ind = j;
					break;
				}
			if (ind == -1)
				continue;
			A.swapRows(i, ind);
		}
		for (int j = i + 1; j < n; ++j)
			if (abs(A[j][i]) >= eps)
			{
				double k = A[j][i] / A[i][i];
				A[j][i] = 0;
				for (int l = i + 1; l < n; ++l)
				{
					A[j][l] -= A[i][l] * k;
					if (abs(A[j][l]) < eps)
						A[j][l] = 0;
				}
			}
	}
	return A;
}

Matrix GaussMethod::forwardElimination(Matrix A, double* b)
{
	int n = A.getN();
	for (int i = 0; i < n; ++i)
	{
		if (abs(A[i][i]) < eps)
		{
			int ind = -1;
			for (int j = i + 1; j < n; ++j)
				if (abs(A[j][i]) >= eps)
				{
					ind = j;
					break;
				}
			if (ind == -1)
				continue;
			A.swapRows(i, ind);
			std::swap(b[i], b[ind]);
		}
		for (int j = i + 1; j < n; ++j)
			if (abs(A[j][i]) >= eps)
			{
				double k = A[j][i] / A[i][i];
				A[j][i] = 0;
				for (int l = i + 1; l < n; ++l)
					A[j][l] -= A[i][l] * k;
				b[j] -= b[i] * k;
			}
	}
	return A;
}

double* GaussMethod::backSubstitution(const Matrix& A, double* b)
{
	int n = A.getN();
	double* x = new double[n];
	for (int i = n - 1; i >= 0; --i)
	{
		double s = 0;
		for (int j = i + 1; j < n; ++j)
			s += A[i][j] * x[j];
		x[i] = (b[i] - s) / A[i][i];
		if (abs(x[i]) < eps)
			x[i] = 0;
	}
	return x;
}

ColumnVector GaussMethod::backSubstitution(const Matrix& A, const ColumnVector& b)
{
	int n = A.getN();
	for (int i = 0; i < n; ++i)
		if (A[i][i] == 0)
			throw std::invalid_argument("СЛАУ имеет не единственное решение");

	ColumnVector x(n);
	for (int i = n - 1; i >= 0; --i)
	{
		double s = 0;
		for (int j = i + 1; j < n; ++j)
			s += A[i][j] * x[j];
		x[i] = (b[i] - s) / A[i][i];
		if (abs(x[i]) < eps)
			x[i] = 0;
	}
	return x;
}

double* GaussMethod::solveSystem(const Matrix& A, double* b)
{
	int n = A.getN();
	double* bt = new double[n];
	for (int i = 0; i < n; ++i)
		bt[i] = b[i];

	Matrix At = forwardElimination(A, bt);

	return backSubstitution(At, bt);
}
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include "Matrix.h"
#include "ColumnVector.h"
#include "RowVector.h"
using namespace std;

const double Matrix::eps = 1e-12;

Matrix::Matrix(double** B, int n) : A(nullptr), n(n), m(n)
{
	A = new double* [n];
	for (int i = 0; i < n; ++i)
	{
		A[i] = new double[n];
		for (int j = 0; j < n; ++j)
			A[i][j] = B[i][j];
	}
}
Matrix::Matrix(double** B, int n, int m) : A(nullptr), n(n), m(m)
{
	A = new double* [n];
	for (int i = 0; i < n; ++i)
	{
		A[i] = new double[m];
		for (int j = 0; j < m; ++j)
			A[i][j] = B[i][j];
	}
}
Matrix::Matrix(int n) : n(n), m(n)
{
	A = new double* [n];
	for (int i = 0; i < n; ++i)
	{
		A[i] = new double[n];
		for (int j = 0; j < n; ++j)
			A[i][j] = 0;
	}
	for (int i = 0; i < n; ++i)
		A[i][i] = 1;
}
Matrix::Matrix(int n, int m) : n(n), m(m)
{
	A = new double* [n];
	for (int i = 0; i < n; ++i)
	{
		A[i] = new double[m];
		for (int j = 0; j < m; ++j)
			A[i][j] = 0;
	}
}
Matrix::Matrix(const char* fname, bool isSquare) : A(nullptr), n(0), m(0)
{
	ifstream fin{ fname };
	if (!fin)
		throw runtime_error("Cannot open file");
	fin >> n;
	if (isSquare)
		m = n;
	else
		fin >> m;
	A = new double* [n];
	for (int i = 0; i < n; ++i)
	{
		A[i] = new double[m];
		for (int j = 0; j < m; ++j)
			fin >> A[i][j];
	}
	fin.close();
}
Matrix::Matrix(const char* fname, double*& b, bool isSquare) : A(nullptr), n(0), m(0)
{
	ifstream fin{ fname };
	if (!fin)
		throw runtime_error("Cannot open file");
	fin >> n;
	if (isSquare)
		m = n;
	else
		fin >> m;
	A = new double* [n];
	for (int i = 0; i < n; ++i)
	{
		A[i] = new double[m];
		for (int j = 0; j < m; ++j)
			fin >> A[i][j];
	}
	b = new double[n];
	for (int i = 0; i < n; ++i)
		fin >> b[i];
	fin.close();
}
Matrix::Matrix(const Matrix& B) : n(B.n), m(B.m)
{
	if (B.empty())
		A = nullptr;
	else
	{
		A = new double* [n];
		for (int i = 0; i < n; ++i)
		{
			A[i] = new double[m];
			for (int j = 0; j < m; ++j)
				A[i][j] = B[i][j];
		}
	}
}
Matrix::Matrix(Matrix&& B) noexcept : A(B.A), n(B.n), m(B.m)
{
	B.A = nullptr;
	B.n = 0;
	B.m = 0;
}

double* Matrix::operator [] (int i)
{
	if (i >= 0 && i < n)
		return A[i];
	throw out_of_range("Index out of range");
}
const double* Matrix::operator [] (int i) const
{
	if (i >= 0 && i < n)
		return A[i];
	throw out_of_range("Index out of range");
}
Matrix& Matrix::operator = (const Matrix& B)
{
	if (this == &B)
		return *this;
	del();
	n = B.n;
	m = B.m;
	if (B.empty())
		A = nullptr;
	else
	{
		A = new double* [n];
		for (int i = 0; i < n; ++i)
		{
			A[i] = new double[m];
			for (int j = 0; j < m; ++j)
				A[i][j] = B[i][j];
		}
	}
	return *this;
}
Matrix& Matrix::operator = (Matrix&& B) noexcept 
{
	if (this == &B)
		return *this;
	del();
	A = B.A;
	n = B.n;
	m = B.m;

	B.A = nullptr;
	B.n = 0;
	B.m = 0;

	return *this;
}
Matrix Matrix::operator * (double k) const
{
	Matrix B(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			B[i][j] = A[i][j] * k;
	return B;
}
Matrix Matrix::operator * (const Matrix& B) const
{
	if (m != B.n)
		throw invalid_argument("Matrix m doesn't match Matrix n");

	int bm = B.m;
	Matrix C(n, bm);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < bm; ++j)
		{
			for (int k = 0; k < m; ++k)
				C[i][j] += A[i][k] * B[k][j];
			if (abs(C[i][j]) < eps)
				C[i][j] = 0;
		}
	}
	return C;
}
ColumnVector Matrix::operator * (const ColumnVector& v) const
{
	if (m != v.getN())
		throw invalid_argument("Matrix m doesn't match ColumnVector n");

	ColumnVector u(n);
	for (int i = 0; i < n; ++i)
	{
		u[i] = 0;
		for (int j = 0; j < m; ++j)
			u[i] += A[i][j] * v[j];
		if (abs(u[i]) < eps)
			u[i] = 0;
	}
	return u;
}
Matrix Matrix::operator + (const Matrix& B) const
{
	if (n != B.n || m != B.m)
		throw invalid_argument("Matrix n and m doesn't match Matrix n and m");

	Matrix C(n, m);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			C[i][j] = A[i][j] + B[i][j];
			if (abs(C[i][j]) < eps)
				C[i][j] = 0;
		}
	}
	return C;
}
Matrix Matrix::operator - (const Matrix& B) const
{
	if (n != B.n || m != B.m)
		throw invalid_argument("Matrix n and m doesn't match Matrix n and m");

	Matrix C(n, m);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			C[i][j] = A[i][j] - B[i][j];
			if (abs(C[i][j]) < eps)
				C[i][j] = 0;
		}
	}
	return C;
}
Matrix Matrix::operator ~ () const
{
	Matrix B(m, n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			B[j][i] = A[i][j];
	return B;
}

ColumnVector Matrix::getColumn(int i) const
{
	if (i < 0 || i >= m)
		throw out_of_range("Column index out of range");

	ColumnVector v(n);
	for (int k = 0; k < n; ++k)
		v[k] = A[k][i];
	return v;
}
RowVector Matrix::getRow(int i) const
{
	if (i >= 0 && i < n)
		return RowVector(A[i], n);
	throw out_of_range("Row index out of range");
}
void Matrix::swapRows(int i, int j)
{
	if (i >= 0 && i < n && j >= 0 && j < n)
		swap(A[i], A[j]);
	else
		throw std::out_of_range("Row index out of range");
}
void Matrix::print() const
{
	if (empty())
		throw invalid_argument("Matrix is empty");

	for (int i = 0; i < n; ++i)
	{
		cout << '(';
		for (int j = 0; j < m - 1; ++j)
			cout << A[i][j] << '\t';
		cout << A[i][m - 1] << ')' << endl;
	}
	cout << endl;
}
void Matrix::print(double* b) const
{
	if (empty())
		throw invalid_argument("Matrix is empty");
	if(b == nullptr)
		throw invalid_argument("b is empty");

	for (int i = 0; i < n; ++i)
	{
		cout << '(';
		for (int j = 0; j < m; ++j)
			cout << A[i][j] << '\t';
		cout << "|\t" << b[i] << ')' << endl;
	}
	cout << endl;
}
void Matrix::print(const ColumnVector& b) const
{
	if (empty())
		throw invalid_argument("Matrix is empty");
	if (b.empty())
		throw invalid_argument("ColumnVector b is empty");

	for (int i = 0; i < n; ++i)
	{
		cout << '(';
		for (int j = 0; j < m; ++j)
			cout << A[i][j] << '\t';
		cout << "|\t" << b[i] << ')' << endl;
	}
	cout << endl;
}

Matrix Matrix::LUdecomposition() const
{
	if (n != m)
		throw invalid_argument("Matrix isn't square => LU-decomposition doesn't exist");

	Matrix LU(n);
	for (int ind = 0; ind < n; ++ind)
	{
		// L
		for (int i = ind; i < n; ++i)
		{
			double s = 0;
			for (int k = 0; k <= ind - 1; ++k)
				s += LU[i][k] * LU[k][ind];
			LU[i][ind] = A[i][ind] - s;
			if (abs(LU[i][ind]) < eps)
				LU[i][ind] = 0;
		}
		if (LU[ind][ind] == 0)
			throw invalid_argument("Matrix has a zero principal minor => LU-decomposition doesn't exist");
		// U
		for (int j = ind + 1; j < n; ++j)
		{
			double s = 0;
			for (int k = 0; k <= ind - 1; ++k)
				s += LU[ind][k] * LU[k][j];
			LU[ind][j] = (A[ind][j] - s) / LU[ind][ind];
			if (abs(LU[ind][j]) < eps)
				LU[ind][j] = 0;
		}
	}
	return LU;
}
pair<Matrix, Matrix> Matrix::QRdecomposition() const
{
	if (n != m)
		throw invalid_argument("Matrix isn't square => QR-decomposition doesn't exist");

	Matrix R(A, n);
	Matrix Q(n), P, Ptmp;
	ColumnVector u;
	for (int i = 0; i < n - 1; ++i)
	{
		int mni = n - i;
		u = ColumnVector(mni);
		bool b = true;
		for (int j = 0; j < mni; ++j)
		{
			u[j] = R[j + i][i];
			if (b && j && abs(u[j]) >= eps)
				b = false;
		}
		if (b)
			continue;

		if (u[0] < 0)
			u[0] += u.norm2();
		else
			u[0] -= u.norm2();
		Ptmp = Matrix(mni) - (u * ~u) * (2.0 / u.scalarSquare());
		P = Matrix(n);
		for (int j = i; j < n; ++j)
			for (int k = i; k < n; ++k)
				P[j][k] = Ptmp[j - i][k - i];
		R = P * R;
		Q = P * Q;
	}
	return make_pair(Q, R);
}
#pragma once
#include <utility>

class ColumnVector;
class RowVector;

class Matrix
{
private:
	static const double eps;
	double** A;
	int n, m;

	inline void del()
	{
		for (int i = 0; i < n; ++i)
			delete[] A[i];
		delete[] A;
		A = nullptr;
		n = m = 0;
	}
public:
	inline Matrix() : A(nullptr), n(0), m(0) {}
	Matrix(double** A, int n);
	Matrix(double** A, int n, int m);
	Matrix(int n);
	Matrix(int n, int m);
	Matrix(const char* fname, bool isSquare = true);
	Matrix(const char* fname, double*& b, bool isSquare = true);
	Matrix(const Matrix& B);
	Matrix(Matrix&& B) noexcept;

	double* operator [] (int i);
	const double* operator [] (int i) const;
	Matrix& operator = (const Matrix& B);
	Matrix& operator = (Matrix&& B) noexcept;
	Matrix operator * (double k) const;
	Matrix operator * (const Matrix& B) const;
	ColumnVector operator * (const ColumnVector& v) const;
	Matrix operator + (const Matrix& B) const;
	Matrix operator - (const Matrix& B) const;
	Matrix operator ~ () const;

	ColumnVector getColumn(int i) const;
	RowVector getRow(int i) const;
	void swapRows(int i, int j);
	inline int getN() const { return n; }
	inline int getM() const { return m; }
	inline bool empty() const { return A == nullptr; }
	void print() const;
	void print(double* b) const;
	void print(const ColumnVector& b) const;

	Matrix LUdecomposition() const; // единицы на главной диагонали матрицы U
	std::pair<Matrix, Matrix> QRdecomposition() const;

	~Matrix() { del(); }
};
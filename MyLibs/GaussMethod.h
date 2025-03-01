#pragma once

class Matrix;
class ColumnVector;

class GaussMethod
{
public:
	static Matrix forwardElimination(Matrix A);
	static Matrix forwardElimination(Matrix A, double* b);

	static double* backSubstitution(const Matrix& A, double* b);
	static ColumnVector backSubstitution(const Matrix& A, const ColumnVector& b);

	static double* solveSystem(const Matrix& A, double* b);
};
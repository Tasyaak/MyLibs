#pragma once

class Matrix;

class TridiagonalMatrixAlgorithm
{
public:
	static double* solveSystem(const Matrix& A, const double* b);
};
#pragma once
#include <utility>

class Matrix;
class ColumnVector;

class QRAlgorithm
{
public:
	static ColumnVector solveSystem(const Matrix& A, double* b);
};
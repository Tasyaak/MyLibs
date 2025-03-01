#include <stdexcept>
#include "ColumnVector.h"
#include "RowVector.h"
#include "Matrix.h"

Matrix ColumnVector::operator * (const RowVector& u) const
{
	int m = u.getN();
	Matrix A(n, m);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < m; ++j)
			A[i][j] = v[i] * u[j];
	return A;
}
RowVector ColumnVector::operator ~ () const { return RowVector(v, n); }
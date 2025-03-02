#include <stdexcept>
#include "RowVector.h"
#include "ColumnVector.h"
#include "Matrix.h"

double RowVector::operator * (const ColumnVector& u) const
{
	if (n != u.size())
		throw std::invalid_argument("ColumnVector n doesn't match RowVector n");

	double s = 0;
	for (int i = 0; i < n; ++i)
		s += v[i] * u[i];
	if (abs(s) < eps)
		s = 0;
	return s;
}
RowVector RowVector::operator * (const Matrix& A) const
{
	if (n != A.getN())
		throw std::invalid_argument("RowVector n doesn't match Matrix n");

	int m = A.getM();
	RowVector u(m);
	for (int i = 0; i < m; ++i)
	{
		u[i] = 0;
		for (int j = 0; j < n; ++j)
			u[i] += v[j] * A[j][i];
		if (abs(u[i]) < eps)
			u[i] = 0;
	}
	return u;
}
ColumnVector RowVector::operator ~ () const { return ColumnVector(v, n); }
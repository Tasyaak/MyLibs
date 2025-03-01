#pragma once
#include "Vector.h"

class Matrix;
class ColumnVector;

class RowVector : public Vector
{
public:
	using Vector::Vector;
	using Vector::operator [];
	using Vector::operator =;

	inline RowVector operator + (const RowVector& u) const { return Vector::operator + <RowVector>(u); }
	inline RowVector operator - (const RowVector& u) const { return Vector::operator - <RowVector>(u); }
	inline RowVector operator * (double k) const { return Vector::operator * <RowVector>(k); }
	inline RowVector operator / (double k) const { return Vector::operator / <RowVector>(k); }

	double operator * (const ColumnVector& u) const;
	RowVector operator * (const Matrix& A) const;
	ColumnVector operator ~ () const;

	~RowVector() { del(); }
};
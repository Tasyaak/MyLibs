#pragma once
#include "Vector.h"

class RowVector;
class Matrix;

class ColumnVector : public Vector
{
public:
	using Vector::Vector;
	using Vector::operator [];
	using Vector::operator =;

	inline ColumnVector operator + (const ColumnVector& u) const { return Vector::operator + <ColumnVector>(u); }
	inline ColumnVector operator - (const ColumnVector& u) const { return Vector::operator - <ColumnVector>(u); }
	inline ColumnVector operator * (double k) const { return Vector::operator * <ColumnVector>(k); }
	inline ColumnVector operator / (double k) const { return Vector::operator / <ColumnVector>(k); }

	Matrix operator * (const RowVector& u) const;
	RowVector operator ~ () const;

	~ColumnVector() { del(); }
};
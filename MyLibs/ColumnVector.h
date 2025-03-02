#pragma once
#include "Vector.h"

class RowVector;
class Matrix;

template <typename T>
class ColumnVector : public Vector<ColumnVector<T>, T>
{
public:
	using Vector<ColumnVector<T>, T>::Vector;

	Matrix operator * (const RowVector& u) const;
	RowVector operator ~ () const;
};
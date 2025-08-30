#pragma once
#include "Vector.hpp"

template <typename T>
class RowVector;
template <typename T>
class Matrix;

template <typename T>
class ColumnVector : public Vector<ColumnVector<T>, T>
{
protected:
	using Vector<ColumnVector<T>, T>::n;
	using Vector<ColumnVector<T>, T>::data;

public:
	using Vector<ColumnVector<T>, T>::Vector;

	Matrix<T> operator * (const RowVector<T>& u) const;
	RowVector<T> operator ~ () const;
};

#include "impl/ColumnVector.tpp"
#pragma once
#include "Vector.hpp"

template <typename T>
class Matrix;
template <typename T>
class ColumnVector;

template <typename T>
class RowVector : public Vector<RowVector<T>, T>
{
protected:
	using Vector<RowVector<T>, T>::n;
	using Vector<RowVector<T>, T>::data;

public:
	using Vector<RowVector<T>, T>::Vector;

	T operator * (const ColumnVector<T>& u) const;
	RowVector<T> operator * (const Matrix<T>& A) const;
	ColumnVector<T> operator ~ () const;
};

#include "RowVector.tpp"
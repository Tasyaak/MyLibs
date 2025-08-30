#pragma once

template<typename T>
class Matrix;
template<typename T>
class ColumnVector;

template<typename T>
class QRAlgorithm
{
public:
	static T* solveSystem(const Matrix<T>& A, T* b, const double eps = 1e-9);
	static ColumnVector<T> solveSystem(const Matrix<T>& A, const ColumnVector<T>& b, const double eps = 1e-9);
};

#include "impl/QRAlgorithm.tpp"
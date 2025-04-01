#pragma once

template<typename T>
class Matrix;
template<typename T>
class ColumnVector;

template<typename T>
class GaussMethod
{
public:
	static void forwardElimination(Matrix<T>& A, const double eps = 1e-9);
	static void forwardElimination(Matrix<T>& A, T* b, const double eps = 1e-9);
	static void forwardElimination(Matrix<T>& A, ColumnVector<T>& b, const double eps = 1e-9);

	static T* backSubstitution(const Matrix<T>& A, T* b, const double eps = 1e-9);
	static ColumnVector<T> backSubstitution(const Matrix<T>& A, const ColumnVector<T>& b, const double eps = 1e-9);

	static T* solveSystem(Matrix<T> A, T* b, const double eps = 1e-9);
	static ColumnVector<T> solveSystem(Matrix<T> A, ColumnVector<T> b, const double eps = 1e-9);
};

#include "details/GaussMethod.tpp"
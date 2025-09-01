#pragma once
#include "core/config.hpp"

template<typename T>
class Matrix;
template<typename T>
class ColumnVector;

template<typename T>
class TridiagonalMatrixAlgorithm
{
public:
	static T* solveSystem(const Matrix<T>& A, T* b);
	static ColumnVector<T> solveSystem(const Matrix<T>& A, const ColumnVector<T>& b);
};

#include "impl/TridiagonalMatrixAlgorithm.tpp"
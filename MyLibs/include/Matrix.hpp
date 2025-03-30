#pragma once
#include <string>
#include <initializer_list>
#include <cassert>
#include <utility>
#include "underlyingType.hpp"

template<typename T>
class ColumnVector;
template<typename T>
class RowVector;

template<typename T>
class Matrix
{
private:
	T* data;
	std::size_t n, m;

public:
	using matrix_tag = void;

	explicit Matrix(std::size_t n, std::size_t m, mathDetails::NoInitTag no) : data(new T[n * m]), n(n), m(m) {}
	explicit Matrix() : data(nullptr), n(0), m(0) {}
	explicit Matrix(T** ptrMatrix, std::size_t n);
	explicit Matrix(T** ptrMatrix, std::size_t n, std::size_t m);
	explicit Matrix(T* ptrMatrix, std::size_t n, std::size_t m);
	explicit Matrix(std::size_t n);
	explicit Matrix(std::size_t n, std::size_t m);
	explicit Matrix(std::initializer_list<std::initializer_list<T>> init);
	Matrix(const Matrix& B);
	Matrix(Matrix&& B) noexcept;
	template<typename Expr, typename = std::enable_if_t<is_matrix_expression<Expr>::value>>
	Matrix(const Expr& expr);

	T* operator [] (std::size_t i);
	const T* operator [] (std::size_t i) const;
	T& operator () (std::size_t i, std::size_t j);
	const T& operator () (std::size_t i, std::size_t j) const;
	Matrix& operator = (const Matrix& B);
	Matrix& operator = (Matrix&& B) noexcept;
	template<typename Expr, typename = std::enable_if_t<is_matrix_expression<Expr>::value>>
	Matrix& operator = (const Expr& expr);
	Matrix& operator += (const Matrix& A);
	Matrix& operator -= (const Matrix& A);
	Matrix& operator *= (const T& scalar);
	Matrix operator * (const Matrix& B) const;
	ColumnVector<T> operator * (const ColumnVector<T>& v) const;
	Matrix operator ~ () const;

	RowVector<T> getRow(std::size_t i) const;
	ColumnVector<T> getCol(std::size_t j) const;
	void swapRows(std::size_t k, std::size_t l);
	void swapCols(std::size_t k, std::size_t l);
	std::size_t numRows() const { return n; }
	std::size_t numCols() const { return m; }
	bool empty() const { return n == 0 || m == 0; }
	std::string toString() const;
	void print() const;
	void print(T* b) const;
	void print(const ColumnVector<T>& b) const;

	Matrix LUdecomposition() const; // единицы на главной диагонали матрицы U
	std::pair<Matrix, Matrix> QRdecomposition(const double eps = 1e-9) const;

	~Matrix();
};

template<typename LHS, typename RHS,
	typename = std::enable_if_t<is_matrix_or_expression<LHS>::value>,
	typename = std::enable_if_t<is_matrix_or_expression<RHS>::value>>
auto operator + (const LHS& lhs, const RHS& rhs) -> mathDetails::MatBinaryOp<LHS, RHS, mathDetails::Add>
{
	using LeftUnderlying = underlying_matrix_t<LHS>;
	using RightUnderlying = underlying_matrix_t<RHS>;

	static_assert(std::is_same<LeftUnderlying, RightUnderlying>::value, "Matrix addition error: operands must have the same underlying type");
	assert(!lhs.empty() && "Matrix addition error: left operand is empty");
	assert(!rhs.empty() && "Matrix addition error: right operand is empty");
	assert(lhs.numRows() == rhs.numRows() && lhs.numCols() == rhs.numCols() && "Matrix addition error: mismatched sizes");

	return mathDetails::MatBinaryOp<LHS, RHS, mathDetails::Add>(lhs, rhs);
}

template<typename LHS, typename RHS,
	typename = std::enable_if_t<is_matrix_or_expression<LHS>::value>,
	typename = std::enable_if_t<is_matrix_or_expression<RHS>::value>>
auto operator - (const LHS& lhs, const RHS& rhs) -> mathDetails::MatBinaryOp<LHS, RHS, mathDetails::Subtract>
{
	using LeftUnderlying = underlying_matrix_t<LHS>;
	using RightUnderlying = underlying_matrix_t<RHS>;

	static_assert(std::is_same<LeftUnderlying, RightUnderlying>::value, "Matrix subtract error: operands must have the same underlying type");
	assert(!lhs.empty() && "Matrix subtract error: left operand is empty");
	assert(!rhs.empty() && "Matrix subtract error: right operand is empty");
	assert(lhs.numRows() == rhs.numRows() && lhs.numCols() == rhs.numCols() && "Matrix subtract error: mismatched sizes");

	return mathDetails::MatBinaryOp<LHS, RHS, mathDetails::Subtract>(lhs, rhs);
}

template<typename Expr, typename T,
	typename = std::enable_if_t<is_matrix_or_expression<Expr>::value>>
auto operator * (const Expr& expr, const T& scalar) -> mathDetails::MatScalarOp<Expr, T, mathDetails::Multiply>
{
	assert(!expr.empty() && "Matrix is empty");

	return mathDetails::MatScalarOp<Expr, T, mathDetails::Multiply>(expr, scalar);
}
template<typename Expr, typename T,
	typename = std::enable_if_t<is_matrix_or_expression<Expr>::value>>
	auto operator * (const T& scalar, const Expr& expr) -> mathDetails::MatScalarOp<Expr, T, mathDetails::Multiply>
{
	assert(!expr.empty() && "Matrix is empty");

	return mathDetails::MatScalarOp<Expr, T, mathDetails::Multiply>(expr, scalar);
}

template<typename Expr, typename T,
	typename = std::enable_if_t<is_matrix_or_expression<Expr>::value>>
auto operator / (const Expr& expr, const T& scalar) -> mathDetails::MatScalarOp<Expr, T, mathDetails::Divide>
{
	assert(!expr.empty() && "Matrix is empty");

	return mathDetails::MatScalarOp<Expr, T, mathDetails::Divide>(expr, scalar);
}

#include "Matrix.tpp"
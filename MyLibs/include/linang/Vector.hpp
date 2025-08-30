#pragma once
#include <string>
#include <initializer_list>
#include <cassert>
#include "details/underlyingType.hpp"
#include "details/expr.hpp"
#include "details/traits.hpp"


template<typename Derived, typename T>
class Vector
{
protected:
	T* data;
	std::size_t n;

public:
	using vector_tag = void;

	explicit Vector(std::size_t n, mathDetails::NoInitTag no) : data(new T[n]), n(n) {}
	explicit Vector() : data(nullptr), n(0) {}
	explicit Vector(T* ptrArray, std::size_t n);
	explicit Vector(std::size_t n);
	explicit Vector(std::initializer_list<T> init);
	Vector(const Vector& v);
	Vector(Vector&& v) noexcept;
	template<typename VecExpr, typename = std::enable_if_t<is_vector_expression<VecExpr>::value>>
	Vector(const VecExpr& expr);

	T& operator [] (std::size_t i);
	const T& operator [] (std::size_t i) const;
	Vector& operator = (const Vector& v);
	Vector& operator = (Vector&& v) noexcept;
	template<typename VecExpr, typename = std::enable_if_t<is_vector_expression<VecExpr>::value>>
	Vector& operator = (const VecExpr& expr);
	Derived operator - () const;
	Derived& operator += (const Derived& v);
	Derived& operator -= (const Derived& v);
	template<typename U, typename = std::enable_if_t<std::is_convertible_v<U, T>>>
	Derived& operator *= (const U& scalar);

	double normInfinity() const;
	double norm1() const;
	double norm2() const;
	double scalarSquare() const;
	T dotProduct(const Vector& v) const;
	bool isZero() const;
	bool empty() const { return n == 0; }
	std::size_t size() const { return n; }
	std::string toString() const;
	void print() const;

	~Vector();
};

template<typename VecExpr>
double normInfinity(const VecExpr& expr);
template<typename VecExpr>
double norm1(const VecExpr& expr);
template<typename VecExpr>
double norm2(const VecExpr& expr);

template<typename LHS, typename RHS,
	typename = std::enable_if_t<is_vector_or_expression<LHS>::value>,
	typename = std::enable_if_t<is_vector_or_expression<RHS>::value>>
auto operator + (const LHS& lhs, const RHS& rhs) -> mathDetails::VecBinaryOp<LHS, RHS, mathDetails::Add>
{
	using LeftUnderlying = underlying_vector_t<LHS>;
	using RightUnderlying = underlying_vector_t<RHS>;

	static_assert(std::is_same<LeftUnderlying, RightUnderlying>::value, "Vector addition error: operands must have the same underlying type");
	assert(!lhs.empty() && "Vector addition error: left operand is empty");
	assert(!rhs.empty() && "Vector addition error: right operand is empty");
	assert(lhs.size() == rhs.size() && "Vector addition error: mismatched sizes");

	return mathDetails::VecBinaryOp<LHS, RHS, mathDetails::Add>(lhs, rhs);
}

template<typename LHS, typename RHS,
	typename = std::enable_if_t<is_vector_or_expression<LHS>::value>,
	typename = std::enable_if_t<is_vector_or_expression<RHS>::value>>
auto operator - (const LHS& lhs, const RHS& rhs) -> mathDetails::VecBinaryOp<LHS, RHS, mathDetails::Subtract>
{
	using LeftUnderlying = underlying_vector_t<LHS>;
	using RightUnderlying = underlying_vector_t<RHS>;

	static_assert(std::is_same<LeftUnderlying, RightUnderlying>::value, "Vector subtract error: operands must have the same underlying type");
	assert(!lhs.empty() && "Vector subtract error: left operand is empty");
	assert(!rhs.empty() && "Vector subtract error: right operand is empty");
	assert(lhs.size() == rhs.size() && "Vector subtract error: mismatched sizes");

	return mathDetails::VecBinaryOp<LHS, RHS, mathDetails::Subtract>(lhs, rhs);
}

template<typename VecExpr, typename T,
	typename = std::enable_if_t<is_vector_or_expression<VecExpr>::value>>
auto operator * (const VecExpr& expr, const T& scalar) -> mathDetails::VecScalarOp<VecExpr, T, mathDetails::Multiply>
{
	assert(!expr.empty() && "Vector is empty");

	return mathDetails::VecScalarOp<VecExpr, T, mathDetails::Multiply>(expr, scalar);
}
template<typename VecExpr, typename T,
	typename = std::enable_if_t<is_vector_or_expression<VecExpr>::value>>
auto operator * (const T& scalar, const VecExpr& expr) -> mathDetails::VecScalarOp<VecExpr, T, mathDetails::Multiply>
{
	assert(!expr.empty() && "Vector is empty");

	return mathDetails::VecScalarOp<VecExpr, T, mathDetails::Multiply>(expr, scalar);
}

template<typename VecExpr, typename T,
	typename = std::enable_if_t<is_vector_or_expression<VecExpr>::value>>
auto operator / (const VecExpr& expr, const T& scalar) -> mathDetails::VecScalarOp<VecExpr, T, mathDetails::Divide>
{
	assert(!expr.empty() && "Vector is empty");

	return mathDetails::VecScalarOp<VecExpr, T, mathDetails::Divide>(expr, scalar);
}

#include "impl/Vector.tpp"
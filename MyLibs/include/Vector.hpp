#pragma once
#include <string>
#include <initializer_list>
#include <cassert>
#include "underlyingType.hpp"

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
	template<typename Expr, typename = std::enable_if_t<is_vector_expression<Expr>::value>>
	Vector(const Expr& expr);

	T& operator [] (std::size_t i);
	const T& operator [] (std::size_t i) const;
	Vector& operator = (const Vector& v);
	Vector& operator = (Vector&& v) noexcept;
	template<typename Expr, typename = std::enable_if_t<is_vector_expression<Expr>::value>>
	Vector& operator = (const Expr& expr);
	Derived operator - () const;
	Derived& operator += (const Derived& v);
	Derived& operator -= (const Derived& v);
	Derived& operator *= (const T& scalar);

	double normInfinity() const;
	double norm1() const;
	double norm2() const;
	T scalarSquare() const;
	T dotProduct(const Vector& v) const;
	bool isZero() const;
	bool empty() const { return n == 0; }
	std::size_t size() const { return n; }
	std::string toString() const;
	void print() const;

	~Vector();
};

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

template<typename Expr, typename T,
	typename = std::enable_if_t<is_vector_or_expression<Expr>::value>>
auto operator * (const Expr& expr, const T& scalar) -> mathDetails::VecScalarOp<Expr, T, mathDetails::Multiply>
{
	assert(!expr.empty() && "Vector is empty");

	return mathDetails::VecScalarOp<Expr, T, mathDetails::Multiply>(expr, scalar);
}
template<typename Expr, typename T,
	typename = std::enable_if_t<is_vector_or_expression<Expr>::value>>
auto operator * (const T& scalar, const Expr& expr) -> mathDetails::VecScalarOp<Expr, T, mathDetails::Multiply>
{
	assert(!expr.empty() && "Vector is empty");

	return mathDetails::VecScalarOp<Expr, T, mathDetails::Multiply>(expr, scalar);
}

template<typename Expr, typename T,
	typename = std::enable_if_t<is_vector_or_expression<Expr>::value>>
auto operator / (const Expr& expr, const T& scalar) -> mathDetails::VecScalarOp<Expr, T, mathDetails::Divide>
{
	assert(!expr.empty() && "Vector is empty");

	return mathDetails::VecScalarOp<Expr, T, mathDetails::Divide>(expr, scalar);
}

#include "Vector.tpp"
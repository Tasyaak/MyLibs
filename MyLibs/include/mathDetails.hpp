#pragma once
namespace mathDetails
{
	struct NoInitTag {};

	struct Add
	{
		template<typename T>
		static T apply(const T& a, const T& b) { return a + b; }
	};

	struct Subtract
	{
		template<typename T>
		static T apply(const T& a, const T& b) { return a - b; }
	};

	struct Multiply
	{
		template<typename T>
		static T apply(const T& a, const T& scalar) { return a * scalar; }
	};

	struct Divide
	{
		template<typename T>
		static T apply(const T& a, const T& scalar) { return a / scalar; }
	};

	template<typename LHS, typename RHS, typename Op>
	class VecBinaryOp
	{
	private:
		const LHS& lhs;
		const RHS& rhs;

	public:
		using vector_expr_tag = void;

		VecBinaryOp(const LHS& lhs, const RHS& rhs) : lhs(lhs), rhs(rhs) {}

		auto operator[](std::size_t i) const -> decltype(Op::apply(lhs[i], rhs[i]))
		{
			return Op::apply(lhs[i], rhs[i]);
		}

		bool empty() const { return lhs.size() == 0; }
		std::size_t size() const { return lhs.size(); }
	};

	template<typename Expr, typename T, typename Op>
	class VecScalarOp
	{
	private:
		const Expr& expr;
		T scalar;

	public:
		using vector_expr_tag = void;

		VecScalarOp(const Expr& expr, T scalar) : expr(expr), scalar(scalar) {}

		auto operator[](std::size_t i) const -> decltype(Op::apply(expr[i], scalar))
		{
			return Op::apply(expr[i], scalar);
		}

		bool empty() const { return expr.size() == 0; }
		std::size_t size() const { return expr.size(); }
	};

	template<typename LHS, typename RHS, typename Op>
	class MatBinaryOp
	{
	private:
		const LHS& lhs;
		const RHS& rhs;

	public:
		using matrix_expr_tag = void;

		MatBinaryOp(const LHS& lhs, const RHS& rhs) : lhs(lhs), rhs(rhs) {}

		auto operator()(std::size_t i, std::size_t j) const -> decltype(Op::apply(lhs(i, j), rhs(i, j)))
		{
			return Op::apply(lhs(i, j), rhs(i, j));
		}

		bool empty() const { return lhs.numRows() == 0 || lhs.numCols() == 0; }
		std::size_t numRows() const { return lhs.numRows(); }
		std::size_t numCols() const { return lhs.numCols(); }
	};

	template<typename Expr, typename T, typename Op>
	class MatScalarOp
	{
	private:
		const Expr& expr;
		T scalar;

	public:
		using matrix_expr_tag = void;

		MatScalarOp(const Expr& expr, T scalar) : expr(expr), scalar(scalar) {}

		auto operator()(std::size_t i, std::size_t j) const -> decltype(Op::apply(expr(i, j), scalar))
		{
			return Op::apply(expr(i, j), scalar);
		}

		bool empty() const { return expr.numRows() == 0 || expr.numCols() == 0; }
		std::size_t numRows() const { return expr.numRows(); }
		std::size_t numCols() const { return expr.numCols(); }
	};
}
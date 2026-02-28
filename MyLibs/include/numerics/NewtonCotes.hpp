#pragma once
#include <cstddef>      // std::size_t
#include <type_traits>  // std::enable_if_t, std::is_invocable_r_v
#include <functional>   // std::invoke
#include <utility>      // std::forward

// The order of accuracy equals degree_ + 1 when degree_ is even, and degree_ + 2 when degree_ is odd
template <class Real = double>
class NewtonCotes
{
	static_assert(std::is_floating_point_v<Real>, "NewtonCotes<Real>: Real must be a floating-point type");
	using real_type = Real;

public:
	static constexpr std::size_t min_degree = 1;
	/*
		* @param degree  n >= 1: rule uses (n+1) equally spaced nodes per panel
		* @param panels  M >= 1: number of panels in composite rule
	*/
	NewtonCotes() = delete;
	NewtonCotes(const NewtonCotes&) = delete;
	NewtonCotes& operator=(const NewtonCotes&) = delete;
	NewtonCotes(NewtonCotes&& other) noexcept
		: nodes_unit_(other.nodes_unit_), weights_unit_(other.weights_unit_),
		degree_(other.degree_), panels_(other.panels_)
	{
		other.nodes_unit_ = nullptr;
		other.weights_unit_ = nullptr;
		other.degree_ = 0;
		other.panels_ = 0;
	}
	NewtonCotes& operator=(NewtonCotes&& other) noexcept
	{
		if (this == &other) return *this;
		delete[] nodes_unit_;
		delete[] weights_unit_;
		nodes_unit_ = other.nodes_unit_;
		weights_unit_ = other.weights_unit_;
		degree_ = other.degree_;
		panels_ = other.panels_;
		other.nodes_unit_ = nullptr;
		other.weights_unit_ = nullptr;
		other.degree_ = 0;
		other.panels_ = 0;
		return *this;
	}
	explicit NewtonCotes(std::size_t degree, std::size_t panels = 1) : degree_(degree), panels_(panels)
	{
		set_degree(degree);
		set_panels(panels);
	}

	[[nodiscard]] std::size_t degree() const noexcept { return degree_; }
	[[nodiscard]] std::size_t panels() const noexcept { return panels_; }

	void set_panels(std::size_t panels); // веса не пересчитываются
	void set_degree(std::size_t degree); // пересчитываем базовые узлы/веса

	/// Доступ к узлам/весам базового правила на [0;1] (узлы равномерные)
	[[nodiscard]] const Real* nodes_unit() const noexcept { return nodes_unit_; }		// size = degree_+1
	[[nodiscard]] const Real* weights_unit() const noexcept { return weights_unit_; }	// size = degree_+1

	template <class F>
	[[nodiscard]] std::enable_if_t<std::is_invocable_r_v<Real, F, Real>, Real>
		integrate(F&& f, Real a, Real b) const;

	// Удобный вызов как функтора
	template <class F>
	[[nodiscard]] std::enable_if_t<std::is_invocable_r_v<Real, F, Real>, Real>
		operator()(F&& f, Real a, Real b) const
	{
		return integrate(std::forward<F>(f), a, b);
	}

	~NewtonCotes()
	{
		delete[] nodes_unit_;
		delete[] weights_unit_;
	}

private:
	Real* nodes_unit_ = nullptr;
	Real* weights_unit_ = nullptr;
	std::size_t degree_ = 0;
	std::size_t panels_ = 0;

	void build_rule_(); // Строит nodes_unit_ и weights_unit_ для closed Newton–Cotes на [0;1]
};

#include "impl/NewtonCotes.tpp"
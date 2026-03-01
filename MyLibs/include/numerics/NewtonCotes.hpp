#pragma once
#include <cstddef>      // std::size_t
#include <type_traits>  // std::enable_if_t, std::is_invocable_r_v
#include <utility>      // std::forward

template <class Real = double>
class NewtonCotes
{
	static_assert(std::is_floating_point_v<Real>, "NewtonCotes<Real>: Real must be a floating-point type");
	using real_type = Real;

public:
	static constexpr std::size_t min_degree = 1;

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

	void set_panels(std::size_t panels);
	void set_degree(std::size_t degree);

	[[nodiscard]] const Real* nodes_unit() const noexcept { return nodes_unit_; }		// size = degree_+1
	[[nodiscard]] const Real* weights_unit() const noexcept { return weights_unit_; }	// size = degree_+1

	template <class F>
	[[nodiscard]] std::enable_if_t<std::is_invocable_r_v<Real, F, Real>, Real>
		integrate_fixed(F&& f, Real a, Real b) const;

	template <class F>
	[[nodiscard]] std::enable_if_t<std::is_invocable_r_v<Real, F, Real>, Real>
		integrate(F&& f, Real a, Real b, Real tol, std::size_t max_depth = 25) const;


	template <class F>
	[[nodiscard]] std::enable_if_t<std::is_invocable_r_v<Real, F, Real>, Real>
		operator()(F&& f, Real a, Real b, Real tol, std::size_t max_depth = 25) const
	{
		return integrate(std::forward<F>(f), a, b, tol, max_depth);
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
	std::size_t p = 0;	// order of accuracy p = degree_ + 1 when degree_ is even, p = degree_ + 2 when degree_ is odd

	void build_rule_(); // nodes_unit_ è weights_unit_ for closed Newton–Cotes on [0;1]

	template <class F>
	[[nodiscard]] std::enable_if_t<std::is_invocable_r_v<Real, F, Real>, Real>
		integrate_panel(F&& f, Real a, Real b, std::size_t panels_local) const;

	template <class F>
	[[nodiscard]] std::enable_if_t<std::is_invocable_r_v<Real, F, Real>, Real>
		adaptive_rec(F&& f, Real a, Real b, Real tol, Real Q1, std::size_t depth, std::size_t panels_local) const;
};

#include "impl/NewtonCotes.tpp"
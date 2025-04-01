#pragma once

template<typename T>
class ColumnVector;
template<typename T>
class Matrix;

template<typename StateType>
class AdamsMoultonMethod
{
private:
	StateType(*F)(double, const StateType&);

	StateType G(const StateType& y, const StateType& yn, double tn, double h, const StateType& f0, const StateType& f1, const StateType& f2, const StateType& f3) const;

	using ReturnType = std::conditional_t<std::is_same_v<StateType, double>, double, Matrix<double>>;
	ReturnType dGdy(const StateType& y, double tn, double h, double H = 1e-6) const;

	void NewtonMethod(StateType& y, const StateType& yn, double tn, double h, const StateType& f0, const StateType& f1, const StateType& f2, const StateType& f3) const;
public:
	double t0, T;

	AdamsMoultonMethod(StateType(*F)(double, const StateType&), double t0, double T) : F(F), t0(t0), T(T) {}

	std::vector<StateType> solve(double h, const StateType& y0) const;
};

#include "details/AdamsMoultonMethod.tpp"
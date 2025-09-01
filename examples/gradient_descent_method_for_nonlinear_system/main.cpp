#include <iostream>
#include "linalg/Matrix.hpp"
#include "linalg/ColumnVector.hpp"

// The gradient descent method for nonlinear systems finds a solution by reducing the problem to minimizing a nonnegative function
// The global minimum of this function is zero, which corresponds to the solution of the system
// It works by moving step by step in the direction of the negative gradient to decrease the function value

std::size_t n;

ColumnVector<double> F(const ColumnVector<double>& x)
{
	ColumnVector<double> f(n);
	// Examples of left sides of the algebraic system of nonlinear equations
	switch (n)
	{
	case 1:
		f[0] = x[0] * x[0] - 4;
		break;
	case 2:
		f[0] = 1.0 / (x[0] - exp(sin(x[0]))) - x[1] * x[1] + 26;
		f[1] = exp(1 - cos(x[1]-5)) - 1 - x[0];
		break;
	case 3:
		f[0] = x[0] * x[1] + sin(x[2] - 3) - 2;
		f[1] = x[0] - 0.25 * x[1] * x[1] - x[2] + 3;
		f[2] = 1.0 / x[0] + x[2] - 4;
		break;
	}
	return f;
}

inline double Phi(const ColumnVector<double>& x) { return F(x).scalarSquare(); }

Matrix<double> Jacobi(ColumnVector<double>(&F)(const ColumnVector<double>&), const ColumnVector<double>& x, const double h = 1e-7)
{
	Matrix<double> J(n);
	ColumnVector<double> xh1, xh2, f1, f2;
	xh1 = xh2 = x;

	for (int i = 0; i < n; ++i)
	{
		xh1[i] += h;
		xh2[i] -= h;

		f1 = F(xh1);
		f2 = F(xh2);
		for (int j = 0; j < n; ++j)
			J(i, j) = (f1[j] - f2[j]) / (2 * h);

		xh1[i] -= h;
		xh2[i] += h;
	}
	return J;
}

int main()
{
	std::cout << "Enter the system dimension: ";
	std::cin >> n;
	std::cout << "Enter the initial guess: ";
	ColumnVector<double> x(n);
	for (std::size_t i = 0; i < n; ++i)
		std::cin >> x[i];

	const double eps = 1e-12, delta = 1e-6;
	ColumnVector<double> xprev, grad;
	double alpha, phiprev, phi, c, gradsqr;
	bool b;
	phi = Phi(x);
	while (true)
	{
		xprev = x;
		grad = Jacobi(F, xprev) * F(xprev) * 2;
		gradsqr = grad.scalarSquare();
		phiprev = phi;

		if (phiprev >= 100)
			c = 0.1;
		else
			c = 0.001;
		alpha = 2;
		do
		{
			x = xprev - grad * alpha;
			alpha *= 0.6;
			phi = Phi(x);
		} while (phi > phiprev - c * alpha * gradsqr);

		if (phiprev - phi < eps)
		{
			if (phi < delta)
				b = true;
			else
				b = false;
			break;
		}
	}
	if (b)
	{
		std::cout << "Solution: ";
		x.print();
		std::cout << "Residual = " << F(x).norm2();
	}
	else
		std::cout << "No solution found for the given initial guess";
}
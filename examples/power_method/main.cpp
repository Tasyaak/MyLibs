#include <iostream>
#include "linalg/ColumnVector.hpp"
#include "linalg/Matrix.hpp"

// Power Method (Power iterations) finds the eigenvalue with the largest magnitude
// and its corresponding eigenvector

int main()
{
	srand(time(0));
	const double eps = 1e-11;
	std::size_t n, count;
	std::cout << "Enter the order and the square matrix: ";
	std::cin >> n;
	Matrix<double> A(n);
	for (std::size_t i = 0; i < n; ++i)
		for (std::size_t j = 0; j < n; ++j)
			std::cin >> A[i][j];
	ColumnVector<double> x, xprev, y(n);
	for (std::size_t i = 0; i < n; ++i)
		y[i] = 15000 - rand();
	std::cout << "Random initial vector: ";
	y.print();
	x = y * (1.0 / y.normInfinity());
	double l = 0, lprev;
	do
	{
		xprev = x;
		lprev = l;
		y = A * xprev;
		x = y * (1.0 / y.normInfinity());
		l = count = 0;
		for (std::size_t i = 0; i < n; ++i)
			if (abs(xprev[i]) >= eps)
			{
				l += y[i] / xprev[i];
				++count;
			}
		l /= count;
	} while (abs(l - lprev) >= eps);

	std::cout << "The eigenvalue with maximum modulus = " << l << ",\n";
	std::cout << "and the corresponding eigenvector: ";
	x.print();
	std::cout << "Residual = " << norm2(A * x - l * x);
}
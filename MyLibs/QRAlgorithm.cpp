#include <iostream>
#include <cmath>
#include <stdexcept>
#include "QRAlgorithm.h"
#include "GaussMethod.h"
#include "Matrix.h"
#include "ColumnVector.h"
using namespace std;

const double eps = 1e-14;

ColumnVector QRAlgorithm::solveSystem(const Matrix& A, double* b)
{
	int n = A.getN();

	pair<Matrix, Matrix> QR = A.QRdecomposition();
	Matrix Q = QR.first, R = QR.second;
	ColumnVector B(b, n);
	B = Q * B;
	try
	{
		return GaussMethod::backSubstitution(R, B);
	}
	catch (const invalid_argument& e)
	{
		cout << "A:\n";
		A.print();
		cout << "Q:\n";
		Q.print();
		cout << "R:\n";
		R.print();
		cout << "\n\n";
		return nullptr;
	}
}
#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
using namespace std;

const double Vector::eps = 1e-12;

Vector::Vector(double* u, int n) : v(nullptr), n(n)
{
	v = new double[n];
	for (int i = 0; i < n; ++i)
		v[i] = u[i];
}
Vector::Vector(int n) : n(n)
{
	v = new double[n];
	for (int i = 0; i < n; ++i)
		v[i] = 0;
}
Vector::Vector(char* fname) : v(nullptr), n(0)
{
	ifstream fin{ fname };
	if (!fin)
		throw runtime_error("Cannot open file");
	fin >> n;
	v = new double[n];
	for (int i = 0; i < n; ++i)
		fin >> v[i];
	fin.close();
}
Vector::Vector(const Vector& u) : v(nullptr), n(u.n)
{
	if (u.empty())
		throw invalid_argument("Second vector is empty");
	v = new double[n];
	for (int i = 0; i < n; ++i)
		v[i] = u[i];
}
Vector::Vector(Vector&& u) noexcept : v(u.v), n(u.n)
{
	u.v = nullptr;
	u.n = 0;
}

double& Vector::operator [] (int i)
{
	if (i >= 0 && i < n)
		return v[i];
	throw out_of_range("Index out of range");
}
const double& Vector::operator [] (int i) const
{
	if (i >= 0 && i < n)
		return v[i];
	throw out_of_range("Index out of range");
}

Vector& Vector::operator = (const Vector& u)
{
	if (this != &u)
	{
		delete[] v;
		n = u.n;
		if (u.empty())
			v = nullptr;
		else
		{
			v = new double[n];
			for (int i = 0; i < n; ++i)
				v[i] = u[i];
		}
	}
	return *this;
}
Vector& Vector::operator = (Vector&& u) noexcept
{
	if (this != &u)
	{
		delete[] v;
		v = u.v;
		n = u.n;

		u.v = nullptr;
		u.n = 0;
	}
	return *this;
}
Vector Vector::operator - () const
{
	Vector u(n);
	for (int i = 0; i < n; ++i)
		u[i] = -v[i];
	return u;
}
//Vector Vector::operator + (const Vector& u) const
//{
//	if (u.empty())
//		throw invalid_argument("Second vector is empty");
//	if (n != u.n)
//		throw invalid_argument("Vector n doesn't match Vector n");
//
//	Vector w(n);
//	for (int i = 0; i < n; ++i)
//	{
//		w[i] = v[i] + u[i];
//		if (abs(w[i]) < eps)
//			w[i] = 0;
//	}
//	return w;
//}
//Vector Vector::operator * (double k) const
//{
//	Vector u(n);
//	for (int i = 0; i < n; ++i)
//		u[i] = k * v[i];
//	return u;
//}
//Vector Vector::operator / (double k) const
//{
//	Vector u(n);
//	for (int i = 0; i < n; ++i)
//		u[i] = v[i] / k;
//	return u;
//}

double Vector::normInfinity() const
{
	if (empty())
		throw invalid_argument("Vector is empty");

	double mx = abs(v[0]);
	for (int i = 1; i < n; ++i)
		if (mx < abs(v[i]))
			mx = abs(v[i]);
	return mx;
}
double Vector::norm1() const
{
	if (empty())
		throw invalid_argument("Vector is empty");

	double sum = 0;
	for (int i = 0; i < n; ++i)
		sum += abs(v[i]);
	return sum;
}
double Vector::scalarSquare() const
{
	if (empty())
		throw invalid_argument("Vector is empty");

	double s = 0;
	for (int i = 0; i < n; ++i)
		s += v[i] * v[i];
	return s;
}
double Vector::dotProduct(const Vector& u) const
{
	if (u.empty())
		throw invalid_argument("Second vector is empty");
	if (n != u.n)
		throw invalid_argument("Vector n doesn't match Vector n");

	double sum = 0;
	for (int i = 0; i < n; ++i)
		sum += v[i] * u[i];
	return sum;
}
bool Vector::isZero() const
{
	bool flag = true;
	for (int i = 0; i < n; ++i)
		if (v[i] != 0)
		{
			flag = false;
			break;
		}
	return flag;
}
void Vector::print() const
{
	if (empty())
		throw invalid_argument("Vector is empty");

	cout << '(';
	for (int i = 0; i < n - 1; ++i)
		cout << v[i] << ' ';
	cout << v[n - 1] << ')';
}
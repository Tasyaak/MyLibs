#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <complex>
#include "Vector.hpp"

template<typename Derived, typename T>
Vector<Derived, T>::Vector(T* ptrArray, std::size_t n) : data(nullptr), n(n)
{
	assert(n > 0 && "Vector size must be greater than 0");
	assert(ptrArray != nullptr && "Pointer cannot be equal nullptr");

	data = new T[n];
	std::copy(ptrArray, ptrArray + n, data);
}
template<typename Derived, typename T>
Vector<Derived, T>::Vector(std::size_t n) : data(nullptr), n(n)
{
	assert(n > 0 && "Vector size must be greater than 0");

	data = new T[n]{};
}
template<typename Derived, typename T>
Vector<Derived, T>::Vector(std::initializer_list<T> init) : data(nullptr), n(init.size())
{
	assert(n > 0 && "Cannot construct vector from an empty initializer list");

	data = new T[n];
	std::copy(init.begin(), init.end(), data);
}
template<typename Derived, typename T>
Vector<Derived, T>::Vector(const Vector& v) : data(nullptr), n(v.n)
{
	assert(!v.empty() && "Vector is empty");

	data = new T[n];
	std::copy(v.data, v.data + n, data);
}
template<typename Derived, typename T>
Vector<Derived, T>::Vector(Vector&& v) noexcept : data(v.data), n(v.n)
{
	v.data = nullptr;
	v.n = 0;
}
template<typename Derived, typename T>
template<typename VecExpr, typename /* SFINAE-параметр omitted */>
Vector<Derived, T>::Vector(const VecExpr& expr) : data(nullptr), n(0)
{
	assert(!expr.empty() && "Cannot construct vector from an empty expression");

	n = expr.size();
	std::size_t size = n;
	data = new T[size];
	for (std::size_t i = 0; i < size; ++i)
		data[i] = expr[i];
}
template<typename Derived, typename T>
Vector<Derived, T>::~Vector()
{
	delete[] data;
	data = nullptr;
	n = 0;
}

template<typename Derived, typename T>
T& Vector<Derived, T>::operator [] (std::size_t i)
{
	assert(!empty() && "Vector is empty");
	assert(i < n && "Index vector out of range");

	return data[i];
}
template<typename Derived, typename T>
const T& Vector<Derived, T>::operator [] (std::size_t i) const
{
	assert(!empty() && "Vector is empty");
	assert(i < n && "Index vector out of range");

	return data[i];
}
template<typename Derived, typename T>
Vector<Derived, T>& Vector<Derived, T>::operator = (const Vector& v)
{
	if (this != &v)
	{
		delete[] data;
		n = v.n;
		if (v.empty())
			data = nullptr;
		else
		{
			data = new T[n];
			std::copy(v.data, v.data + n, data);
		}
	}
	return *this;
}
template<typename Derived, typename T>
Vector<Derived, T>& Vector<Derived, T>::operator = (Vector&& v) noexcept
{
	if (this != &v)
	{
		delete[] data;
		data = v.data;
		n = v.n;
		v.data = nullptr;
		v.n = 0;
	}
	return *this;
}
template<typename Derived, typename T>
template<typename VecExpr, typename /* SFINAE-параметр omitted */>
Vector<Derived, T>& Vector<Derived, T>::operator = (const VecExpr& expr)
{
	assert(!expr.empty() && "Cannot copy vector from an empty expression");
	assert(n == expr.size() && "Vector size doesn't match expression size");

	std::size_t size = n;
	for (std::size_t i = 0; i < size; ++i)
		data[i] = expr[i];
	return *this;
}
template<typename Derived, typename T>
Derived Vector<Derived, T>::operator - () const
{
	assert(!empty() && "Vector is empty");

	std::size_t size = n;
	Derived v(size, mathDetails::NoInitTag{});
	for (std::size_t i = 0; i < size; ++i)
		v.data[i] = -data[i];
	return v;
}
template<typename Derived, typename T>
Derived& Vector<Derived, T>::operator += (const Derived& v)
{
	assert(!empty() && "First vector is empty");
	assert(!v.empty() && "Second vector is empty");
	assert(n == v.n && "First vector size doesn't match second vector size");

	std::size_t size = n;
	for (std::size_t i = 0; i < size; ++i)
		data[i] += v.data[i];
	return static_cast<Derived&>(*this);
}
template<typename Derived, typename T>
Derived& Vector<Derived, T>::operator -= (const Derived& v)
{
	assert(!empty() && "First vector is empty");
	assert(!v.empty() && "Second vector is empty");
	assert(n == v.n && "First vector size doesn't match second vector size");

	std::size_t size = n;
	for (std::size_t i = 0; i < size; ++i)
		data[i] -= v.data[i];
	return static_cast<Derived&>(*this);
}
template<typename Derived, typename T>
template<typename U, typename /* SFINAE-параметр omitted */>
Derived& Vector<Derived, T>::operator *= (const U& scalar)
{
	assert(!empty() && "Vector is empty");

	std::size_t size = n;
	T scalarT = static_cast<T>(scalar);
	for (std::size_t i = 0; i < size; ++i)
		data[i] *= scalarT;
	return static_cast<Derived&>(*this);
}
template<typename Derived, typename T>
double Vector<Derived, T>::normInfinity() const
{
	assert(!empty() && "Vector is empty");

	std::size_t size = n;
	double mx = static_cast<double>(std::abs(data[0]));
	for (std::size_t i = 1; i < size; ++i)
		mx = std::max(mx, static_cast<double>(std::abs(data[i])));
	return mx;
}
template<typename Derived, typename T>
double Vector<Derived, T>::norm1() const
{
	assert(!empty() && "Vector is empty");

	std::size_t size = n;
	double sum = 0.0;
	for (std::size_t i = 0; i < size; ++i)
		sum += static_cast<double>(std::abs(data[i]));
	return sum;
}
template<typename Derived, typename T>
double Vector<Derived, T>::norm2() const
{
	assert(!empty() && "Vector is empty");

	std::size_t size = n;
	double sum = 0.0;
	for (std::size_t i = 0; i < size; ++i)
	{
		double a = static_cast<double>(std::abs(data[i]));
		sum += a * a;
	}
	return std::sqrt(sum);
}
template<typename Derived, typename T>
double Vector<Derived, T>::scalarSquare() const
{
	assert(!empty() && "Vector is empty");

	std::size_t size = n;
	double sum = 0.0;
	for (std::size_t i = 0; i < size; ++i)
	{
		double a = static_cast<double>(std::abs(data[i]));
		sum += a * a;
	}
	return sum;
}

template<typename U>
inline auto myConj(const U& value) -> std::enable_if_t<std::is_floating_point_v<U>, U>
{
	return value;
}

template<typename U>
inline auto myConj(const std::complex<U>& value) -> std::complex<U>
{
	return std::conj(value);
}

template<typename Derived, typename T>
T Vector<Derived, T>::dotProduct(const Vector& v) const
{
	assert(!empty() && "First vector is empty");
	assert(!v.empty() && "Second vector is empty");
	assert(n == v.n && "First vector size doesn't match second vector size");

	std::size_t size = n;
	T sum = T{};
	for (std::size_t i = 0; i < size; ++i)
		sum += data[i] * myConj(v.data[i]);
	return sum;
}
template<typename Derived, typename T>
bool Vector<Derived, T>::isZero() const
{
	assert(!empty() && "Vector is empty");

	std::size_t size = n;
	bool flag = true;
	T temp = T{};
	for (std::size_t i = 0; i < size; ++i)
		if (data[i] != temp)
		{
			flag = false;
			break;
		}
	return flag;
}
template<typename Derived, typename T>
std::string Vector<Derived, T>::toString() const
{
	if (empty())
		return "[]";

	std::size_t size = n;
	std::ostringstream oss;
	oss << '[' << data[0];
	for (std::size_t i = 1; i < size; ++i)
		oss << ' ' << data[i];
	oss << ']';
	return oss.str();
}
template<typename Derived, typename T>
void Vector<Derived, T>::print() const
{
	std::cout << toString() << '\n';
}

template<typename VecExpr>
double normInfinity(const VecExpr& expr)
{
	assert(!expr.empty() && "Vector expression is empty");

	std::size_t size = expr.size();
	double mx = static_cast<double>(std::abs(expr[0]));
	for (std::size_t i = 1; i < size; ++i)
		mx = std::max(mx, static_cast<double>(std::abs(expr[i])));
	return mx;
}
template<typename VecExpr>
double norm1(const VecExpr& expr)
{
	assert(!expr.empty() && "Vector expression is empty");

	std::size_t size = expr.size();
	double sum = 0.0;
	for (std::size_t i = 0; i < size(); ++i)
		sum += static_cast<double>(std::abs(expr[i]));
	return sum;
}
template<typename VecExpr>
double norm2(const VecExpr& expr)
{
	assert(!expr.empty() && "Vector expression is empty");

	std::size_t size = expr.size();
	double sum = 0.0;
	for (std::size_t i = 0; i < size; ++i)
	{
		double a = static_cast<double>(std::abs(expr[i]));
		sum += a * a;
	}
	return std::sqrt(sum);
}
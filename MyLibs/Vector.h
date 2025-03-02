#pragma once
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

template<typename Derived, typename T>
class Vector
{
private:
	struct NoInitTag {};

protected:
	static const double eps = 1e-12;
	T* data;
	int n;

	explicit Vector(int n, NoInitTag no) : data(new T[n]), n(n) {}

public:
	explicit Vector() : data(nullptr), n(0) {}
	explicit Vector(T* v, int n) : data(nullptr), n(n)
	{
		if (n <= 0)
			throw std::invalid_argument("Vector size n cannot be equal " + std::to_string(n));

		data = new T[n];
		std::copy(v, v + n, data);
	}
	explicit Vector(int n) : data(nullptr), n(n)
	{
		if (n <= 0)
			throw std::invalid_argument("Vector size n cannot be equal " + std::to_string(n));

		data = new T[n]();
	}
	Vector(std::initializer_list<T> init) : n(static_cast<int>(init.size())), data(new T[n])
	{
		std::copy(init.begin(), init.end(), data);
	}
	Vector(const Vector& v) : data(nullptr), n(v.n)
	{
		if (empty())
			throw std::invalid_argument("First vector is empty");
		if (v.empty())
			throw std::invalid_argument("Second vector is empty");

		data = new T[n];
		std::copy(v, v + n, data);
	}
	Vector(Vector&& v) noexcept : data(v.data), n(v.n)
	{
		v.data = nullptr;
		v.n = 0;
	}

	T& operator [] (int i)
	{
		if (empty())
			throw std::invalid_argument("Vector is empty");
		if (i < 0 || i >= n)
			throw std::out_of_range("Index vector out of range, ind = " + std::to_string(i));

		return data[i];
	}
	const T& operator [] (int i) const
	{
		if (empty())
			throw std::invalid_argument("Vector is empty");
		if (i < 0 || i >= n)
			throw std::out_of_range("Index vector out of range, ind = " + std::to_string(i));

		return data[i];
	}
	Vector& operator = (const Vector& v)
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
				std::copy(v, v + n, data);
			}
		}
		return *this;
	}
	Vector& operator = (Vector&& v) noexcept
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
	Derived operator - () const
	{
		Derived v(n, NoInitTag{});
		for (int i = 0; i < n; ++i)
			v.data[i] = -data[i];
		return v;
	}
	Derived operator + (const Derived& v) const
	{
		if (empty())
			throw std::invalid_argument("First vector is empty");
		if (v.empty())
			throw std::invalid_argument("Second vector is empty");
		if (n != v.n)
			throw std::invalid_argument("First vector size " + std::to_string(n) + " doesn't match second vector size " + std::to_string(v.n));

		Derived u(n, NoInitTag{});
		for (int i = 0; i < n; ++i)
		{
			u.data[i] = data[i] + v.data[i];
			if (abs(u.data[i]) < eps)
				u.data[i] = T();
		}
		return u;
	}
	Derived& operator += (const Derived& v)
	{
		if (empty())
			throw std::invalid_argument("First vector is empty");
		if (v.empty())
			throw std::invalid_argument("Second vector is empty");
		if (n != v.n)
			throw std::invalid_argument("First vector size " + std::to_string(n) + " doesn't match second vector size " + std::to_string(v.n));

		for (int i = 0; i < n; ++i)
			data[i] += v.data[i];
		return static_cast<Derived&>(*this);
	}
	Derived operator - (const Derived& v) const
	{
		if (empty())
			throw std::invalid_argument("First vector is empty");
		if (v.empty())
			throw std::invalid_argument("Second vector is empty");
		if (n != v.n)
			throw std::invalid_argument("First vector size " + std::to_string(n) + " doesn't match second vector size " + std::to_string(v.n));

		Derived u(n, NoInitTag{});
		for (int i = 0; i < n; ++i)
		{
			u.data[i] = data[i] - v.data[i];
			if (abs(u.data[i]) < eps)
				u.data[i] = T();
		}
		return u;
	}
	Derived& operator -= (const Derived& v)
	{
		if (empty())
			throw std::invalid_argument("First vector is empty");
		if (v.empty())
			throw std::invalid_argument("Second vector is empty");
		if (n != v.n)
			throw std::invalid_argument("First vector size " + std::to_string(n) + " doesn't match second vector size " + std::to_string(v.n));

		for (int i = 0; i < n; ++i)
			data[i] -= v.data[i];
		return static_cast<Derived&>(*this);
	}
	Derived operator * (const T& scalar) const
	{
		if (empty())
			throw std::invalid_argument("Vector is empty");

		Derived v(n, NoInitTag{});
		for (int i = 0; i < n; ++i)
			v.data[i] = scalar * data[i];
		return v;
	}
	Derived operator * (const T& scalar, const Derived& v) { return v * scalar; }
	Derived& operator *= (const T& scalar)
	{
		if (empty())
			throw std::invalid_argument("Vector is empty");

		for (int i = 0; i < n; ++i)
			data[i] *= scalar;
		return static_cast<Derived&>(*this);
	}
	Derived operator / (const T& scalar) const
	{
		if (empty())
			throw std::invalid_argument("Vector is empty");

		Derived v(n, NoInitTag{});
		for (int i = 0; i < n; ++i)
			v.data[i] = data[i] / scalar;
		return v;
	}
	Derived& operator /= (const T& scalar)
	{
		if (empty())
			throw std::invalid_argument("Vector is empty");

		for (int i = 0; i < n; ++i)
			data[i] /= scalar;
		return static_cast<Derived&>(*this);
	}

	double normInfinity() const
	{
		if (empty())
			throw std::invalid_argument("Vector is empty");

		double mx = abs(data[0]);
		for (int i = 1; i < n; ++i)
			if (mx < abs(data[i]))
				mx = abs(data[i]);
		return mx;
	}
	double norm1() const
	{
		if (empty())
			throw std::invalid_argument("Vector is empty");

		double sum = 0;
		for (int i = 0; i < n; ++i)
			sum += abs(data[i]);
		return sum;
	}
	double norm2() const
	{
		if (empty())
			throw std::invalid_argument("Vector is empty");

		double sum = 0;
		for (int i = 0; i < n; ++i)
			sum += abs(data[i]) * abs(data[i]);
		return sqrt(sum);
	}
	T scalarSquare() const
	{
		if (empty())
			throw std::invalid_argument("Vector is empty");

		T sum = T();
		for (int i = 0; i < n; ++i)
			sum += data[i] * data[i];
		return sum;
	}
	T dotProduct(const Vector& v) const
	{
		if (empty())
			throw std::invalid_argument("First vector is empty");
		if (v.empty())
			throw std::invalid_argument("Second vector is empty");
		if (n != v.n)
			throw std::invalid_argument("First vector size " + std::to_string(n) + " doesn't match second vector size " + std::to_string(v.n));

		T sum = T();
		for (int i = 0; i < n; ++i)
			sum += data[i] * v.data[i];
		return sum;
	}
	bool isZero() const
	{
		if (empty())
			throw std::invalid_argument("Vector is empty");

		bool flag = true;
		T temp + T();
		for (int i = 0; i < n; ++i)
			if (v[i] != temp)
			{
				flag = false;
				break;
			}
		return flag;
	}
	bool empty() const { return data == nullptr; }
	int size() const { return n; }
	void print() const { cout << toString() << '\n'; }
	std::string toString() const
	{
		if (empty())
			return "[]";

		std::ostringstream oss;
		oss << '[' << data[0];
		for (int i = 1; i < n; ++i)
			oss << ' ' << data[i];
		oss << ']';
		return oss.str();
	}

	virtiual ~Vector()
	{
		delete[] data;
		data = nullptr;
		n = 0;
	}
};
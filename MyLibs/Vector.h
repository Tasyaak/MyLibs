#pragma once
#include <stdexcept>

class Vector
{
protected:
	static const double eps;
	double* v;
	int n;

	void del()
	{
		delete[] v;
		v = nullptr;
		n = 0;
	}
public:
	inline Vector() : v(nullptr), n(0) {}
	Vector(double* u, int n);
	Vector(int n);
	Vector(char* fname);
	Vector(const Vector& u);
	Vector(Vector&& u) noexcept;

	double& operator [] (int i);
	const double& operator[](int i) const;
	Vector& operator = (const Vector& u);
	Vector& operator = (Vector&& u) noexcept;
	Vector operator - () const;

	template <typename Derived>
	inline Derived operator + (const Derived& u) const
	{
		if (u.empty())
			throw std::invalid_argument("Second vector is empty");
		if (n != u.n)
			throw std::invalid_argument("Vector n doesn't match Vector n");

		Derived w(n);
		for (int i = 0; i < n; ++i)
		{
			w[i] = v[i] + u[i];
			if (abs(w[i]) < eps)
				w[i] = 0;
		}
		return w;
	}

	template <typename Derived>
	inline Derived operator - (const Derived& u) const
	{
		if (u.empty())
			throw std::invalid_argument("Second vector is empty");
		if (n != u.n)
			throw std::invalid_argument("Vector n doesn't match Vector n");

		Derived w(n);
		for (int i = 0; i < n; ++i)
		{
			w[i] = v[i] - u[i];
			if (abs(w[i]) < eps)
				w[i] = 0;
		}
		return w;
	}

	template <typename Derived>
	inline Derived operator * (double k) const
	{
		Derived u(n);
		for (int i = 0; i < n; ++i)
			u[i] = k * v[i];
		return u;
	}

	template <typename Derived>
	inline Derived operator / (double k) const
	{
		Derived u(n);
		for (int i = 0; i < n; ++i)
			u[i] = v[i] / k;
		return u;
	}

	double normInfinity() const;
	double norm1() const;
	inline double norm2() const { return sqrt(scalarSquare()); }
	double scalarSquare() const;
	double dotProduct(const Vector& v) const;
	bool isZero() const;
	inline bool empty() const { return v == nullptr; }
	inline int getN() const { return n; }
	void print() const;

	~Vector() { del(); }
};
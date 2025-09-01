#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cmath>
#include "../Matrix.hpp"
#include "../ColumnVector.hpp"
#include "../RowVector.hpp"
#include "core/complexUtils.hpp"

template<typename T>
Matrix<T>::Matrix(T** ptrMatrix, std::size_t n) : data(nullptr), n(n), m(n)
{
	assert(n > 0 && "Matrix number of rows must be greater than 0");
	assert(ptrMatrix != nullptr && "Pointer cannot be equal nullptr");

	data = new T[n * n];
	for (size_t i = 0; i < n; ++i)
		std::copy(ptrMatrix[i], ptrMatrix[i] + n, data + i * n);
}
template<typename T>
Matrix<T>::Matrix(T** ptrMatrix, std::size_t n, std::size_t m) : data(nullptr), n(n), m(m)
{
	assert(n > 0 && "Matrix number of rows must be greater than 0");
	assert(m > 0 && "Matrix number of columns must be greater than 0");
	assert(ptrMatrix != nullptr && "Pointer cannot be equal nullptr");

	data = new T[n * m];
	for (size_t i = 0; i < n; ++i)
		std::copy(ptrMatrix[i], ptrMatrix[i] + m, data + i * m);
}
template<typename T>
Matrix<T>::Matrix(T* ptrMatrix, std::size_t n) : data(nullptr), n(n), m(n)
{
	assert(n > 0 && "Matrix number of rows must be greater than 0");
	assert(ptrMatrix != nullptr && "Pointer cannot be equal nullptr");

	data = new T[n * n];
	std::copy(ptrMatrix, ptrMatrix + n * n, data);
}
template<typename T>
Matrix<T>::Matrix(T* ptrMatrix, std::size_t n, std::size_t m) : data(nullptr), n(n), m(m)
{
	assert(n > 0 && "Matrix number of rows must be greater than 0");
	assert(m > 0 && "Matrix number of columns must be greater than 0");
	assert(ptrMatrix != nullptr && "Pointer cannot be equal nullptr");

	data = new T[n * m];
	std::copy(ptrMatrix, ptrMatrix + n * m, data);
}
template<typename T>
Matrix<T>::Matrix(std::size_t n) : data(nullptr), n(n), m(n)
{
	assert(n > 0 && "Matrix number of rows must be greater than 0");
	assert(m > 0 && "Matrix number of columns must be greater than 0");

	data = new T[n * n]{};
	for (std::size_t i = 0; i < n; ++i)
		data[i * n + i] = T{ 1 };
}
template<typename T>
Matrix<T>::Matrix(std::size_t n, std::size_t m) : data(nullptr), n(n), m(m)
{
	assert(n > 0 && "Matrix number of rows must be greater than 0");
	assert(m > 0 && "Matrix number of columns must be greater than 0");

	data = new T[n * m]{};
}
template<typename T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> init) : data(nullptr), n(init.size()), m(0)
{
	if (n == 0)
		return;

	m = init.begin()->size();
	for (const auto& row : init)
		assert(row.size() == m && "All rows must have the same size");

	data = new T[n * m];
	size_t i = 0;
	for (auto row : init)
	{
		std::copy(row.begin(), row.begin() + m, data + i * m);
		++i;
	}
}
template<typename T>
Matrix<T>::Matrix(const Matrix& A) : data(nullptr), n(A.n), m(A.m)
{
	assert(!A.empty() && "Matrix is empty");

	data = new T[n * m];
	std::copy(A.data, A.data + n * m, data);
}
template<typename T>
Matrix<T>::Matrix(Matrix&& A) noexcept : data(A.data), n(A.n), m(A.m)
{
	A.data = nullptr;
	A.n = A.m = 0;
}
template<typename T>
template<typename MatExpr, typename /* SFINAE-параметр omitted */>
Matrix<T>::Matrix(const MatExpr& expr) : data(nullptr), n(0), m(0)
{
	assert(!expr.empty() && "Cannot construct matrix from an empty expression");

	n = expr.numRows();
	m = expr.numCols();
	data = new T[n * m];
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		for (std::size_t j = 0; j < m; ++j)
			data[im + j] = expr(i, j);
	}
}
template<typename T>
Matrix<T>::~Matrix()
{
	delete[] data;
	data = nullptr;
	n = m = 0;
}

template<typename T>
T* Matrix<T>::operator [] (std::size_t i)
{
	assert(!empty() && "Matrix is empty");
	assert(i < n && "Index Matrix rows out of range");

	return data + i * m;
}
template<typename T>
const T* Matrix<T>::operator [] (std::size_t i) const
{
	assert(!empty() && "Matrix is empty");
	assert(i < n && "Index Matrix rows out of range");
	
	return data + i * m;
}
template<typename T>
T& Matrix<T>::operator () (std::size_t i, std::size_t j)
{
	assert(!empty() && "Matrix is empty");
	assert(i < n && "Index Matrix rows out of range");
	assert(j < m && "Index Matrix columns out of range");

	return data[i * m + j];
}
template<typename T>
const T& Matrix<T>::operator () (std::size_t i, std::size_t j) const
{
	assert(!empty() && "Matrix is empty");
	assert(i < n && "Index Matrix rows out of range");
	assert(j < m && "Index Matrix columns out of range");

	return data[i * m + j];
}
template<typename T>
Matrix<T>& Matrix<T>::operator = (const Matrix& A)
{
	if (this != &A)
	{
		delete[] data;
		n = A.n;
		m = A.m;
		if (A.empty())
			data = nullptr;
		else
		{
			data = new T[n * m];
			std::copy(A.data, A.data + n * m, data);
		}
	}
	return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::operator = (Matrix&& A) noexcept
{
	if (this != &A)
	{
		delete[] data;
		data = A.data;
		n = A.n;
		m = A.m;
		A.data = nullptr;
		A.n = A.m = 0;
	}
	return *this;
}
template<typename T>
template<typename MatExpr, typename /* SFINAE-параметр omitted */>
Matrix<T>& Matrix<T>::operator = (const MatExpr& expr)
{
	assert(!expr.empty() && "Cannot copy matrix from an empty expression");
	assert(n == expr.numRows() && m == expr.numCols() && "Matrix size doesn't match expression size");

	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		for (std::size_t j = 0; j < m; ++j)
			data[im + j] = expr(i, j);
	}
	return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::operator += (const Matrix& A)
{
	assert(!empty() && "First matrix is empty");
	assert(!A.empty() && "Second matrix is empty");
	assert(n == A.n && m == A.m && "First matrix size doesn't match second matrix size");

	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		for (std::size_t j = 0; j < m; ++j)
			data[im + j] += A.data[im + j];
	}
	return *this;
}
template<typename T>
Matrix<T>& Matrix<T>::operator -= (const Matrix& A)
{
	assert(!empty() && "First matrix is empty");
	assert(!A.empty() && "Second matrix is empty");
	assert(n == A.n && m == A.m && "First matrix size doesn't match second matrix size");

	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		for (std::size_t j = 0; j < m; ++j)
			data[im + j] -= A.data[im + j];
	}
	return *this;
}
template<typename T>
template<typename U, typename /* SFINAE-параметр omitted */>
Matrix<T>& Matrix<T>::operator *= (const U& scalar)
{
	assert(!empty() && "Matrix is empty");

	T scalarT = static_cast<T>(scalar);
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		for (std::size_t j = 0; j < m; ++j)
			data[im + j] *= scalarT;
	}
	return *this;
}
template<typename T>
Matrix<T> Matrix<T>::operator * (const Matrix& A) const
{
	assert(!empty() && "First matrix is empty");
	assert(!A.empty() && "Second matrix is empty");
	assert(m == A.n && "First matrix number of columns doesn't match second matrix number of rows");

	std::size_t newM = A.m;
	Matrix C(n, newM);
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		std::size_t inewM = i * newM;
		for (std::size_t j = 0; j < newM; ++j)
			for (std::size_t k = 0; k < m; ++k)
			{
				std::size_t knewM = k * newM;
				C.data[inewM + j] += data[im + k] * A.data[knewM + j];
			}
	}
	return C;
}
template<typename T>
ColumnVector<T> Matrix<T>::operator * (const ColumnVector<T>& v) const
{
	assert(!empty() && "Matrix is empty");
	assert(!v.empty() && "ColumnVector is empty");
	assert(m == v.size() && "Matrix number of columns doesn't match ColumnVector size");

	ColumnVector<T> u(n);
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		for (std::size_t j = 0; j < m; ++j)
			u[i] += data[im + j] * v[j];
	}
	return u;
}
template<typename T>
Matrix<T> Matrix<T>::operator ~ () const
{
	assert(!empty() && "Matrix is empty");

	Matrix<T> A(m, n, mathDetails::NoInitTag{});
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		for (std::size_t j = 0; j < m; ++j)
		{
			std::size_t jn = j * n;
			A.data[jn + i] = data[im + j];
		}
	}
	return A;
}

template<typename T>
RowVector<T> Matrix<T>::getRow(std::size_t i) const
{
	assert(!empty() && "Matrix is empty");
	assert(i < n && "Index Matrix rows out of range");

	RowVector<T> v(m, mathDetails::NoInitTag{});
	for (std::size_t j = 0; j < m; ++j)
	{
		std::size_t im = i * m;
		v[j] = data[im + j];
	}
	return v;
}
template<typename T>
ColumnVector<T> Matrix<T>::getCol(std::size_t j) const
{
	assert(!empty() && "Matrix is empty");
	assert(j < m && "Index Matrix columns out of range");

	ColumnVector<T> v(n, mathDetails::NoInitTag{});
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		v[i] = data[im + j];
	}
	return v;
}
template<typename T>
void Matrix<T>::swapRows(std::size_t k, std::size_t l)
{
	assert(!empty() && "Matrix is empty");
	assert(k < n && "Index Matrix rows out of range");
	assert(l < n && "Index Matrix rows out of range");

	if (k == l)
		return;

	T* rowK = data + k * m;
	T* rowL = data + l * m;
	std::swap_ranges(rowK, rowK + m, rowL);
}
template<typename T>
void Matrix<T>::swapCols(std::size_t k, std::size_t l)
{
	assert(!empty() && "Matrix is empty");
	assert(k < m && "Index Matrix columns out of range");
	assert(l < m && "Index Matrix columns out of range");

	if (k == l)
		return;

	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		std::swap(data[im + k], data[im + l]);
	}
}
template<typename T>
T Matrix<T>::det(const double eps) const
{
	assert(n == m && "Matrix isn't square for determinant");

	Matrix<T> A(data, n);
	T res = T{ 1 };
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t in = i * n;
		std::size_t ind = i;
		for (std::size_t j = i + 1; j < n; ++j)
			if (abs(A.data[ind * n + i]) < abs(A.data[j * n + i]))
				ind = j;

		if (abs(A.data[ind * n + i]) < eps)
			return T{};

		if (ind != i)
		{
			res *= T{ -1 };
			A.swapRows(i, ind);
		}
		res *= A.data[in + i];
		for (std::size_t j = i + 1; j < n; ++j)
		{
			std::size_t jn = j * n;
			if (abs(A.data[jn + i]) >= eps)
			{
				T k = A.data[jn + i] / A.data[in + i];
				A.data[jn + i] = T{};
				for (std::size_t l = i + 1; l < n; ++l)
					A.data[jn + l] -= A.data[in + l] * k;
			}
		}
	}
	return res;
}
template<typename T>
Matrix<T> Matrix<T>::inverse(const double eps) const
{
	assert(n == m && "Matrix isn't square for determinant");

	Matrix<T> A(data, n);
	Matrix<T> res(n);
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t in = i * n;
		std::size_t ind = i;
		for (std::size_t j = i + 1; j < n; ++j)
			if (abs(A.data[ind * n + i]) < abs(A.data[j * n + i]))
				ind = j;

		assert(abs(A.data[ind * n + i]) >= eps && "Matrix is singular => inverse doesn't exist");

		if (ind != i)
		{
			res.swapRows(i, ind);
			A.swapRows(i, ind);
		}
		for (std::size_t j = i + 1; j < n; ++j)
		{
			std::size_t jn = j * n;
			if (abs(A.data[jn + i]) >= eps)
			{
				T k = A.data[jn + i] / A.data[in + i];
				A.data[jn + i] = T{};

				for (std::size_t l = 0; l <= i; ++l)
					res.data[jn + l] -= res.data[in + l] * k;
				for (std::size_t l = i + 1; l < n; ++l)
				{
					res.data[jn + l] -= res.data[in + l] * k;
					A.data[jn + l] -= A.data[in + l] * k;
				}
			}
		}
	}
	for (std::size_t i = n; i-- > 0;)
	{
		std::size_t in = i * n;

		for (std::size_t j = 0; j < n; ++j)
			res.data[in + j] /= A.data[in + i];
		A.data[in + i] = T{ 1 };

		for (std::size_t j = 0; j < i; ++j)
		{
			std::size_t jn = j * n;
			if (abs(A.data[jn + i]) >= eps)
			{
				for (std::size_t l = 0; l < n; ++l)
					res.data[jn + l] -= res.data[in + l] * A.data[jn + i];

				A.data[jn + i] = T{};
			}
		}
	}
	return res;
}
template<typename T>
std::string Matrix<T>::toString() const
{
	if (empty())
		return "()";

	std::size_t* colWidths = new std::size_t[m]{};
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		for (std::size_t j = 0; j < m; ++j)
		{
			std::ostringstream temp;
			temp << data[im + j];
			colWidths[j] = std::max(colWidths[j], temp.str().size());
		}
	}
	std::ostringstream oss;
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		oss << "(";
		for (std::size_t j = 0; j < m - 1; ++j)
			oss << std::left << std::setw(colWidths[j] + 4) << data[im + j];
		oss << std::left << std::setw(colWidths[m - 1]) << data[im + m - 1];
		oss << ")\n";
	}
	delete[] colWidths;

	return oss.str();
}
template<typename T>
void Matrix<T>::print() const
{
	std::cout << toString() << '\n';
}
template<typename T>
void Matrix<T>::print(T* b) const
{
	assert(!empty() && "Matrix is empty");
	assert(b != nullptr && "Pointer is empty");

	std::size_t* colWidths = new std::size_t[m + 1]{};
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		for (std::size_t j = 0; j < m; ++j)
		{
			std::ostringstream temp;
			temp << data[im + j];
			colWidths[j] = std::max(colWidths[j], temp.str().size());
		}
		{
			std::ostringstream temp;
			temp << b[i];
			colWidths[m] = std::max(colWidths[m], temp.str().size());
		}
	}
	std::ostringstream oss;
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		oss << "(";
		for (std::size_t j = 0; j < m - 1; ++j)
			oss << std::left << std::setw(colWidths[j] + 4) << data[im + j];
		oss << std::left << std::setw(colWidths[m - 1]) << data[im + m - 1] << " | ";
		oss << std::left << std::setw(colWidths[m]) << b[i];
		oss << ")\n";
	}
	delete[] colWidths;

	std::cout << oss.str() << '\n';
}
template<typename T>
void Matrix<T>::print(const ColumnVector<T>& b) const
{
	assert(!empty() && "Matrix is empty");
	assert(!b.empty() && "ColumnVector is empty");
	assert(n == b.size() && "Matrix number of rows doesn't match ColumnVector size");

	std::size_t* colWidths = new std::size_t[m + 1]{};
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		for (std::size_t j = 0; j < m; ++j)
		{
			std::ostringstream temp;
			temp << data[im + j];
			colWidths[j] = std::max(colWidths[j], temp.str().size());
		}
		{
			std::ostringstream temp;
			temp << b[i];
			colWidths[m] = std::max(colWidths[m], temp.str().size());
		}
	}
	std::ostringstream oss;
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t im = i * m;
		oss << "(";
		for (std::size_t j = 0; j < m - 1; ++j)
			oss << std::left << std::setw(colWidths[j] + 4) << data[im + j];
		oss << std::left << std::setw(colWidths[m - 1]) << data[im + m - 1] << " | ";
		oss << std::left << std::setw(colWidths[m]) << b[i];
		oss << ")\n";
	}
	delete[] colWidths;

	std::cout << oss.str() << '\n';
}

template<typename T>
Matrix<T> Matrix<T>::adjoint() const
{
	assert(!empty() && "Matrix is empty");

	Matrix A(data, n, m);
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t in = i * n;
		for (std::size_t j = 0; j < m; ++j)
		{
			std::size_t jn = j * n;
			A.data[jn + i] = myConj(data[in + j]);
		}
	}
	return A;
}
template<typename T>
bool Matrix<T>::isHermitian(const double eps) const
{
	double amax = 0.0;
	for (std::size_t k = 0; k < n; ++k)
		amax = std::max(amax, std::abs(data[k * n + k]));

	double thresh = eps * std::max(1.0, amax);

	for (std::size_t i = 1; i < n; ++i)
	{
		std::size_t in = i * n;
		for (std::size_t j = 0; j < i; ++j)
			if (std::abs(data[in + j] - myConj(data[j * n + i])) > thresh)
				return false;
	}
	return true;
}
template<typename T>
bool Matrix<T>::isPositiveDefinite(const double eps) const
{
	assert(!empty() && "Matrix is empty");
	assert(n == m && "Matrix isn't square => checking for positive definiteness is invalid");
	assert(isHermitian(eps) && "Matrix isn't hermitian => checking for positive definiteness is invalid");

	Matrix<T> A(data, n, n);
	for (std::size_t k = 0; k < n; ++k)
	{
		std::size_t kn = k * n;
		T s = A.data[kn + k];
		for (std::size_t j = 0; j < k; ++j)
			s -= A.data[kn + j] * myConj(A.data[kn + j]);
		if (std::real(s) <= 0 || std::abs(std::imag(s)) > eps * std::real(A.data[kn + k]))
			return false;
		
		A.data[kn + k] = std::sqrt(std::real(s));
		T invAkk = 1 / A.data[kn + k];
		for (std::size_t i = k + 1; i < n; ++i)
		{
			std::size_t in = i * n;
			s = T{};
			for (std::size_t j = 0; j < k; ++j)
				s += A.data[in + j] * myConj(A.data[kn + j]);
			A.data[in + k] = (A.data[in + k] - s) * invAkk;
		}
	}
	return true;
}
template<typename T>
Matrix<T> Matrix<T>::LUdecomposition(const double eps) const
{
	assert(!empty() && "Matrix is empty");
	assert(n == m && "Matrix isn't square => LU-decomposition doesn't exist");

	Matrix<T> LU(n, mathDetails::NoInitTag{});
	for (std::size_t ind = 0; ind < n; ++ind)
	{
		// L
		for (std::size_t i = ind; i < n; ++i)
		{
			std::size_t im = i * m;
			T s = T{};
			for (std::size_t k = 0; k < ind; ++k)
			{
				std::size_t km = k * m;
				s += LU.data[im + k] * LU.data[km + ind];
			}
			LU.data[im + ind] = data[im + ind] - s;
		}
		std::size_t indm = ind * m;
		assert(std::abs(LU.data[indm + ind]) >= eps && "Matrix has a zero principal minor => LU-decomposition doesn't exist");
		// U
		for (std::size_t j = ind + 1; j < n; ++j)
		{
			T s = T{};
			for (std::size_t k = 0; k < ind; ++k)
			{
				std::size_t km = k * m;
				s += LU.data[indm + k] * LU.data[km + j];
			}
			LU.data[indm + j] = (data[indm + j] - s) / LU.data[indm + ind];
		}
	}
	return LU;
}
template<typename T>
std::pair<Matrix<T>, Matrix<T>> Matrix<T>::QRdecomposition(const double eps) const
{
	assert(!empty() && "Matrix is empty");
	assert(n == m && "Matrix isn't square => QR-decomposition doesn't exist");

	Matrix<T> R(data, n);
	Matrix<T> Q(n);
	for (std::size_t i = 0; i < n - 1; ++i)
	{
		std::size_t mni = n - i;
		ColumnVector<T> u(mni, mathDetails::NoInitTag{});
		bool b = true;
		for (std::size_t j = 0; j < mni; ++j)
		{
			std::size_t ijn = (i + j) * n;
			u[j] = R.data[ijn + i];
			if (b && j && std::abs(u[j]) >= eps)
				b = false;
		}
		if (b)
			continue;

		if (u[0] < T{})
			u[0] += u.norm2();
		else
			u[0] -= u.norm2();
		Matrix<T> Ptmp = Matrix<T>(mni) - (u * ~u) * (2.0 / u.scalarSquare());
		Matrix<T> P(n);
		for (std::size_t j = i; j < n; ++j)
		{
			std::size_t jn = j * n, jimni = (j - i) * mni;
			for (std::size_t k = i; k < n; ++k)
				P.data[jn + k] = Ptmp.data[jimni + k - i];
		}
		R = P * R;
		Q = P * Q;
	}
	return std::make_pair(Q, R);
}
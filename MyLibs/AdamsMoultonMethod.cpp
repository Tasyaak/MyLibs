#include <vector>
#include <tuple>
#include "AdamsMoultonMethod.h"
#include "ColumnVector.h"
#include "Matrix.h"
#include "GaussMethod.h"

template<typename StateType>
StateType AdamsMoultonMethod<StateType>::G(const StateType& y, const StateType& yn, double tn, double h, const StateType& f0, const StateType& f1, const StateType& f2, const StateType& f3) const
{
	return y - yn - F(tn, y) * (251 * h / 720) - f3 * (646 * h / 720) + f2 * (264 * h / 720) - f1 * (106 * h / 720) + f0 * (19 * h / 720);
}
template<typename StateType>
typename AdamsMoultonMethod<StateType>::ReturnType AdamsMoultonMethod<StateType>::dGdy(const StateType& y, double tn, double h, double H) const
{
    if constexpr (std::is_same_v<StateType, double>)
        return 1 - h * 251 / 720 * (F(tn, y + H) - F(tn, y - H)) / (2 * H);
	else if constexpr (std::is_same_v<StateType, ColumnVector>)
	{
		int n = y.size();
		Matrix f_prime(n);
		ColumnVector f1, f2, y1 = y, y2 = y;
		for (int j = 0; j < n; ++j)
		{
			y1[j] += H;
			y2[j] -= H;
			f1 = F(tn, y1);
			f2 = F(tn, y2);
			for (int i = 0; i < n; ++i)
				f_prime[i][j] = (f1[i] - f2[i]) / (2 * H);
			y1[j] -= H;
			y2[j] += H;
		}
		return Matrix(n) - f_prime * (251 * h / 720);
	}
}
template<typename StateType>
void AdamsMoultonMethod<StateType>::NewtonMethod(StateType& y, const StateType& yn, double tn, double h, const StateType& f0, const StateType& f1, const StateType& f2, const StateType& f3) const
{
	const double eps = 1e-10;
	StateType prevY, temp;
	if constexpr (std::is_same_v<StateType, double>)
		do
		{
			prevY = y;
			temp = G(y, yn, tn, h, f0, f1, f2, f3);
			y = y - temp / dGdy(y, tn, h);
		} while (abs(y - prevY) >= eps || abs(temp) >= eps);
	else
	{
		Matrix Q, R;
		ColumnVector b1, b2, z;

		temp = G(y, yn, tn, h, f0, f1, f2, f3);
		do
		{
			prevY = y;
			std::tie(Q, R) = dGdy(y, tn, h).QRdecomposition();
			Q = ~Q;

			b1 = R * y - Q * temp;
			z = GaussMethod::backSubstitution(R, b1);

			b2 = R * z - Q * G(z, yn, tn, h, f0, f1, f2, f3);
			y = GaussMethod::backSubstitution(R, b2);
			temp = G(y, yn, tn, h, f0, f1, f2, f3);
		} while ((y - prevY).norm2() >= eps || temp.norm2() >= eps);
	}
}

template<typename StateType>
std::vector<StateType> AdamsMoultonMethod<StateType>::solve(double h, const StateType& y0) const
{
	int n = (T - t0) / h + 1;
	double* t = new double[n];
	std::vector<StateType> y(n);
	std::vector<StateType> f(n);
	t[0] = t0;
	for (int i = 1; i < n; ++i)
		t[i] = t[i - 1] + h;
	y[0] = y0;
	f[0] = F(t0, y0);

	StateType k1, k2, k3, k4, k5, k6;
	for (int i = 0; i < 4; ++i)
	{
		k1 = F(t[i], y[i]) * h;
		k2 = F(t[i] + h / 4, y[i] + k1 / 4) * h;
		k3 = F(t[i] + h / 4, y[i] + (k1 + k2) / 8) * h;
		k4 = F(t[i] + h / 2, y[i] + k3 / 2) * h;
		k5 = F(t[i] + 3 * h / 4, y[i] + k1 * (3.0 / 16) - k2 * (3.0 / 8) + k3 * (3.0 / 8) + k4 * (9.0 / 16)) * h;
		k6 = F(t[i] + h, y[i] - k1 * (3.0 / 7) + k2 * (8.0 / 7) + k3 * (6.0 / 7) - k4 * (12.0 / 7) + k5 * (8.0 / 7)) * h;

		y[i + 1] = y[i] + k1 * (7.0 / 90) + k3 * (16.0 / 45) + k4 * (2.0 / 15) + k5 * (16.0 / 45) + k6 * (7.0 / 90);
		f[i + 1] = F(t[i + 1], y[i + 1]);
	}

	for (int i = 5; i < n; ++i)
	{
		y[i] = y[i - 1] + f[i - 1] * (1901 * h / 720) - f[i - 2] * (2774 * h / 720) + f[i - 3] * (2616 * h / 720) - f[i - 4] * (1274 * h / 720) + f[i - 5] * (251 * h / 720);
		NewtonMethod(y[i], y[i - 1], t[i], h, f[i - 4], f[i - 3], f[i - 2], f[i - 1]);
		f[i] = F(t[i], y[i]);
	}
	return y;
}

template class AdamsMoultonMethod<double>;
template class AdamsMoultonMethod<ColumnVector>;
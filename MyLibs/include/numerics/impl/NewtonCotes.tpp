#pragma once
#include <cassert>
#ifdef __INTELLISENSE__
#include "../NewtonCotes.hpp"
#endif

namespace newton_cotes_detail
{
    template <class Real>
    inline Real* poly_mul_linear(Real* p, std::size_t& n, Real c)
    {
        // multiply p(t) by (t - c)
        ++n;
        Real* r = new Real[n]();
        for (std::size_t k = 0; k < n - 1; ++k)
        {
            r[k + 1] += p[k];    // * t
            r[k] -= c * p[k];    // * (-c)
        }
        delete[] p;
        return r;
    }

    template <class Real>
    inline Real* poly_div_linear_root(Real* P, std::size_t P_size, Real a)
    {
        // P in increasing powers, degree m
        std::size_t m = P_size - 1;
        Real* Q = new Real[m]();
        Real b = P[m];     // leading coefficient
        Q[m - 1] = b;

        for (std::size_t k = m; k-- > 1; )
        {
            b = P[k] + a * b;
            Q[k - 1] = b;
        }
        // remainder would be P[0] + a*b ~ 0
        return Q;
    }

    // Âåñà closed Newton–Cotes íà [0,1] äëÿ degree=n (óçëû i/n)
    template <class Real>
    inline Real* newton_cotes_closed_weights_unit(const std::size_t n)
    {
        // óçëû íà [0..n], ïîòîì ìàñøòàáèðóåì h=1/n ê [0..1]
        const Real h = Real(1) / static_cast<Real>(n);

        // P(t)=\prod_{m=0}^n (t - m)
        std::size_t size = 1;
        Real* P = new Real[size] { Real(1) };
        for (std::size_t m = 0; m <= n; ++m)
            P = poly_mul_linear(P, size, static_cast<Real>(m));

        Real* w = new Real[n + 1]();
        const Real N = static_cast<Real>(n);

        Real* fact = new Real[n + 1];
        fact[0] = Real(1);
        for (size_t k = 1; k <= n; ++k)
            fact[k] = fact[k - 1] * Real(k);

        for (std::size_t i = 0; i <= n; ++i)
        {
            Real* Q = poly_div_linear_root(P, size, static_cast<Real>(i)); // degree n

            // I = \int_0^n Q(t) dt
            Real I = Real(0);

            Real power = N; // N^(k+1) for k=0
            for (std::size_t k = 0; k < size - 1; ++k)
            {
                I += Q[k] * power / static_cast<Real>(k + 1);
                power *= N;
            }
            delete[] Q;

            // P'(i)=(-1)^{n-i} i! (n-i)!
            Real denom = fact[i] * fact[n - i];
            if (((n - i) & 1u) == 1u) denom = -denom;

            const Real omega = I / denom; // weight on [0..n]
            w[i] = h * omega;             // scale to [0..1]
        }
        delete[] P;
        delete[] fact;
        return w;
    }
} // namespace newton_cotes_detail

template <class Real>
void NewtonCotes<Real>::build_rule_()
{
    const std::size_t n = degree_;
    assert ((n <= 12) && "NewtonCotes: degree too large for stable closed Newton–Cotes");

    delete[] nodes_unit_;
    nodes_unit_ = new Real[n + 1];
    for (std::size_t i = 0; i <= n; ++i)
        nodes_unit_[i] = static_cast<Real>(i) / static_cast<Real>(n);

    Real* w_ld = newton_cotes_detail::newton_cotes_closed_weights_unit<Real>(n);
    delete[] weights_unit_;
    weights_unit_ = new Real[n + 1];
    for (std::size_t i = 0; i <= n; ++i)
        weights_unit_[i] = static_cast<Real>(w_ld[i]);
    delete[] w_ld;
}

template <class Real>
void NewtonCotes<Real>::set_panels(std::size_t panels)
{
    assert((panels >= 1) && "NewtonCotes: panels must be >= 1");
    panels_ = panels;
}

template <class Real>
void NewtonCotes<Real>::set_degree(std::size_t degree)
{
    assert((degree >= min_degree) && "NewtonCotes: degree must be >= 1 (closed rule)");
    degree_ = degree;
    build_rule_();
}

template <class Real>
template <class F>
std::enable_if_t<std::is_invocable_r_v<Real, F, Real>, Real>
NewtonCotes<Real>::integrate(F&& f, Real a, Real b) const
{
    assert((a < b) && "NewtonCotes::integrate: expected a < b");

    Real H = (b - a) / static_cast<Real>(panels_);
    Real sum = Real(0);

    for (std::size_t p = 0; p < panels_; ++p)
    {
        const Real left = a + static_cast<Real>(p) * H;

        Real panel_sum = Real(0);
        for (std::size_t i = 0; i <= degree_; ++i) {
            const Real x = left + H * nodes_unit_[i];
            panel_sum += weights_unit_[i] * static_cast<Real>(std::invoke(f, x));
        }
        sum += H * panel_sum;
    }
    return sum;
}
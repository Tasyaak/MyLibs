#pragma once
#include <type_traits>
#include <complex>

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
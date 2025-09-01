#pragma once
#include <type_traits>

namespace mathDetails
{
	struct NoInitTag {};

	struct Add
	{
		template<typename T>
		static T apply(const T& a, const T& b) { return a + b; }
	};

	struct Subtract
	{
		template<typename T>
		static T apply(const T& a, const T& b) { return a - b; }
	};

	struct Multiply
	{
		template<typename T, typename U, typename = std::enable_if_t<std::is_convertible_v<U, T>>>
		static T apply(const T& a, const U& scalar)
		{
			return a * static_cast<T>(scalar);
		}
	};

	struct Divide
	{
		template<typename T, typename U, typename = std::enable_if_t<std::is_convertible_v<U, T>>>
		static T apply(const T& a, const U& scalar)
		{
			return a / static_cast<T>(scalar);
		}
	};
}
#pragma once
#include <type_traits>

// Traits для вектора: по умолчанию false, если нету vector_tag
template<typename T, typename = void>
struct is_vector : std::false_type {};

template<typename T>
struct is_vector<T, std::void_t<typename T::vector_tag>> : std::true_type {};

// Traits для ленивых векторных выражений: по умолчанию false, если нету vector_expr_tag
template<typename T, typename = void>
struct is_vector_expression : std::false_type {};

template<typename T>
struct is_vector_expression<T, std::void_t<typename T::vector_expr_tag>> : std::true_type {};

// Универсальный trait для всех векторных типов (настоящих и выражений)
template<typename T>
struct is_vector_or_expression : std::integral_constant<bool, is_vector<T>::value || is_vector_expression<T>::value> {};

// Traits для матрицы: по умолчанию false, если нету matrix_tag
template<typename T, typename = void>
struct is_matrix : std::false_type {};

template<typename T>
struct is_matrix<T, std::void_t<typename T::matrix_tag>> : std::true_type {};

// Traits для ленивых матричных выражений: по умолчанию false, если нету matrix_expr_tag
template<typename T, typename = void>
struct is_matrix_expression : std::false_type {};

template<typename T>
struct is_matrix_expression<T, std::void_t<typename T::matrix_expr_tag>> : std::true_type {};

// Универсальный trait для всех матричных типов (настоящих и выражений)
template<typename T>
struct is_matrix_or_expression : std::integral_constant<bool, is_matrix<T>::value || is_matrix_expression<T>::value> {};
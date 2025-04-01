#pragma once
#include "traits.hpp"
#include "mathDetails.hpp"

// ������� ������: ��� ���������� ������� ��� ������� ��� ��
template<typename T>
struct underlying_vector
{
    using type = T;
};

// ��� ������� �������� ��������� � ���� underlying �� ������� ��������:
template<typename LHS, typename RHS, typename Op>
struct underlying_vector<mathDetails::VecBinaryOp<LHS, RHS, Op>>
{
    using type = typename underlying_vector<LHS>::type;
};

// ��� ������� �������� � �������� � ����������:
template<typename Expr, typename T, typename Op>
struct underlying_vector<mathDetails::VecScalarOp<Expr, T, Op>>
{
    using type = typename underlying_vector<Expr>::type;
};

template<typename T>
using underlying_vector_t = typename underlying_vector<T>::type;

// ������� ������: ��� ��������� ������� ��� ������� ��� ��
template<typename T>
struct underlying_matrix
{
    using type = T;
};

// ��� ������� �������� ��������� � ���� underlying �� ������� ��������:
template<typename LHS, typename RHS, typename Op>
struct underlying_matrix<mathDetails::MatBinaryOp<LHS, RHS, Op>>
{
    using type = typename underlying_matrix<LHS>::type;
};

// ��� ������� �������� � �������� � ����������:
template<typename Expr, typename T, typename Op>
struct underlying_matrix<mathDetails::MatScalarOp<Expr, T, Op>>
{
    using type = typename underlying_matrix<Expr>::type;
};

template<typename T>
using underlying_matrix_t = typename underlying_matrix<T>::type;
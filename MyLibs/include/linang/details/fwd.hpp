#pragma once
#include <cstddef>

namespace mathDetails
{
    struct NoInitTag;

    template<class LHS, class RHS, class Op> class VecBinaryOp;
    template<class VecExpr, class T, class Op> class VecScalarOp;

    template<class LHS, class RHS, class Op> class MatBinaryOp;
    template<class MatExpr, class T, class Op> class MatScalarOp;
}
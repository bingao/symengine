#include <symengine/basic.h>
#include <symengine/matrices/matrix_expr.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/matrix_derivative.h>

namespace SymEngine
{

MatrixDerivative::MatrixDerivative(const RCP<const MatrixExpr> &arg,
                                   const multiset_basic &x)
: arg_{arg}, x_{x}
{
    SYMENGINE_ASSIGN_TYPEID()
    SYMENGINE_ASSERT(is_canonical(arg, x))
}

hash_t MatrixDerivative::__hash__() const
{
    hash_t seed = SYMENGINE_MATRIXDERIVATIVE;
    hash_combine<Basic>(seed, *arg_);
    for (auto &p : x_) {
        hash_combine<Basic>(seed, *p);
    }
    return seed;
}

bool MatrixDerivative::__eq__(const Basic &o) const
{
    if (is_a<MatrixDerivative>(o)
        and eq(*arg_, *(down_cast<const MatrixDerivative &>(o).arg_))
        and unified_eq(x_, down_cast<const MatrixDerivative &>(o).x_))
        return true;
    return false;
}

int MatrixDerivative::compare(const Basic &o) const
{
    SYMENGINE_ASSERT(is_a<MatrixDerivative>(o))
    const MatrixDerivative &s = down_cast<const MatrixDerivative &>(o);
    int cmp = arg_->__cmp__(*(s.arg_));
    if (cmp != 0)
        return cmp;
    cmp = unified_compare(x_, s.x_);
    return cmp;
}

bool MatrixDerivative::is_canonical(const RCP<const MatrixExpr> &arg,
                                    const multiset_basic &x) const
{
    // Check that 'x' are Symbols:
    for (const auto &a : x)
        if (not is_a<Symbol>(*a))
            return false;
    // Check that 'arg' is MatrixSymbol:
    if (is_a<MatrixSymbol>(*arg))
        return true;
    return false;
}

} // namespace SymEngine

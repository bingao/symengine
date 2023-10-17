#include <symengine/basic.h>
#include <symengine/symengine_exception.h>
#include <symengine/symengine_rcp.h>
#include <symengine/matrices/matrix_expr.h>
#include <symengine/matrices/matrix_symbol.h>
#include <symengine/matrices/matrix_derivative.h>

namespace SymEngine
{

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

RCP<const MatrixExpr> matrix_derivative(const RCP<const MatrixExpr> &arg,
                                        const multiset_basic &x)
{
    if (is_a<MatrixSymbol>(*arg)) {
        return make_rcp<const MatrixDerivative>(arg, x);
    }
    else {
        RCP<const MatrixExpr> result = arg;
        for (const auto &a : x) {
            if (is_a<Symbol>(*a)) {
                result = rcp_dynamic_cast<const MatrixExpr>(
                    result->diff(rcp_dynamic_cast<const Symbol>(a))
                );
            }
            else {
                throw DomainError("Invalid variable type for differentiation.");
            }
        }
        return result;
    }
}

} // namespace SymEngine

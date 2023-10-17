#ifndef SYMENGINE_MATRICES_MATRIX_DERIVATIVE_H
#define SYMENGINE_MATRICES_MATRIX_DERIVATIVE_H

#include <symengine/basic.h>
#include <symengine/matrices/matrix_expr.h>

namespace SymEngine
{

// The class MatrixDerivative is modified from the class Derivative, such that
// it can be used as an argument for classes ConjugateMatrix, Transpose, etc.
class MatrixDerivative : public MatrixExpr
{
private:
    RCP<const MatrixExpr> arg_; //! The matrix expression to be differentiated
    multiset_basic x_; //! x, y, ...

public:
    MatrixDerivative(const RCP<const MatrixExpr> &arg, const multiset_basic &x):
        arg_{arg}, x_{x}
    {
        SYMENGINE_ASSIGN_TYPEID();
        SYMENGINE_ASSERT(is_canonical(arg, x));
    }

    IMPLEMENT_TYPEID(SYMENGINE_MATRIXDERIVATIVE)
    hash_t __hash__() const override;
    bool __eq__(const Basic &o) const override;
    int compare(const Basic &o) const override;
    bool is_canonical(const RCP<const MatrixExpr> &arg,
                      const multiset_basic &x) const;
    inline RCP<const MatrixExpr> get_arg() const
    {
        return arg_;
    }
    inline const multiset_basic &get_symbols() const
    {
        return x_;
    }
    vec_basic get_args() const override
    {
        vec_basic args = {arg_};
        args.insert(args.end(), x_.begin(), x_.end());
        return args;
    }
};

// Helper function to create `MatrixDerivative` that should be always used
// instead of the constructor of `MatrixDerivative`
RCP<const MatrixExpr> matrix_derivative(const RCP<const MatrixExpr> &arg,
                                        const multiset_basic &x);

} // namespace SymEngine

#endif

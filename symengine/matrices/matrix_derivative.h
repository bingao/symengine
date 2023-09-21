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
    IMPLEMENT_TYPEID(SYMENGINE_MATRIXDERIVATIVE)
    MatrixDerivative(const RCP<const MatrixExpr> &arg, const multiset_basic &x);

    static RCP<const MatrixDerivative> create(const RCP<const MatrixExpr> &arg,
                                              const multiset_basic &x)
    {
        return make_rcp<const MatrixDerivative>(arg, x);
    }

    hash_t __hash__() const override;
    bool __eq__(const Basic &o) const override;
    int compare(const Basic &o) const override;
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
    bool is_canonical(const RCP<const MatrixExpr> &arg,
                      const multiset_basic &x) const;
};

} // namespace SymEngine

#endif

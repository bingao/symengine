#include <algorithm>  // std::min_element
#include <cstddef>    // std::size_t
#include <vector>

#include <symengine/basic.h>
#include <symengine/matrices/matrix_expr.h>
#include <symengine/matrices/trace.h>
#include <symengine/visitor.h>

namespace SymEngine
{

hash_t Trace::__hash__() const
{
    hash_t seed = SYMENGINE_TRACE;
    hash_combine<Basic>(seed, *arg_);
    return seed;
}

bool Trace::__eq__(const Basic &o) const
{
    return (is_a<Trace>(o) && arg_->__eq__(*down_cast<const Trace &>(o).arg_));
}

int Trace::compare(const Basic &o) const
{
    SYMENGINE_ASSERT(is_a<Trace>(o));

    return arg_->compare(*down_cast<const Trace &>(o).arg_);
}

vec_basic Trace::get_args() const
{
    return {arg_};
}

class MatrixTraceVisitor : public BaseVisitor<MatrixTraceVisitor>
{
private:
    RCP<const Basic> trace_;

    void trace_error()
    {
        throw DomainError("Trace is only valid for square matrices");
    }

public:
    MatrixTraceVisitor() {}

    void bvisit(const Basic &x){};

    void bvisit(const MatrixExpr &x)
    {
        auto arg = rcp_static_cast<const MatrixExpr>(x.rcp_from_this());
        trace_ = make_rcp<const Trace>(arg);
    }

    void bvisit(const IdentityMatrix &x)
    {
        trace_ = x.size();
    }

    void bvisit(const ZeroMatrix &x)
    {
        tribool sq = is_square(x);
        if (is_true(sq)) {
            trace_ = zero;
        } else if (is_false(sq)) {
            trace_error();
        } else {
            auto arg = rcp_static_cast<const MatrixExpr>(x.rcp_from_this());
            trace_ = make_rcp<const Trace>(arg);
        }
    }

    void bvisit(const DiagonalMatrix &x)
    {
        trace_ = add(x.get_container());
    }

    void bvisit(const ImmutableDenseMatrix &x)
    {
        if (x.nrows() != x.ncols()) {
            trace_error();
        }
        vec_basic diag;
        for (size_t i = 0; i < x.nrows(); i++) {
            diag.push_back(x.get(i, i));
        }
        trace_ = add(diag);
    }

    void bvisit(const MatrixAdd &x)
    {
        // Trace is a linear function so trace(A + B) = trace(A) + trace(B)
        RCP<const Basic> sum = zero;
        for (auto &e : x.get_terms()) {
            e->accept(*this);
            sum = add(sum, trace_);
        }
        trace_ = sum;
    }

    // We need to consider (i) the trace is invariant under circular shifts,
    // (ii) tr(c*A) = c*tr(A) where c is a scalar
    void bvisit(const MatrixMul &x)
    {
        auto scalar = x.get_scalar();
        if (eq(*scalar, *zero)) {
            trace_ = zero;
        } else {
            auto factors = x.get_factors();
            // Rearrange the factors by using the cyclic property
            auto min_factor = std::min_element(
                factors.begin(), factors.end(), RCPBasicKeyLess()
            );
            RCP<const MatrixExpr> product;
            if (min_factor == factors.begin()) {
                product = matrix_mul(factors);
            } else {
                auto sorted_factors = vec_basic(min_factor, factors.end());
                sorted_factors.insert(sorted_factors.end(), factors.begin(), min_factor);
                product = matrix_mul(sorted_factors);
            }
            if (is_a<const MatrixMul>(*product)) {
                trace_ = make_rcp<const Trace>(product);
            } else {
                // Result is saved in `trace_`
                product->accept(*this);
            }
            if (neq(*scalar, *one)) trace_ = mul(trace_, scalar);
        }
    }

    RCP<const Basic> apply(const MatrixExpr &s)
    {
        s.accept(*this);
        return trace_;
    }
};

RCP<const Basic> trace(const RCP<const MatrixExpr> &arg)
{
    MatrixTraceVisitor visitor;
    return visitor.apply(*arg);
}
} // namespace SymEngine

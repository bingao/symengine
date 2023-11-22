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
    // (ii) tr(c*A) = c*tr(A) where c is a scalar, and (iii) tr((A+B)*C) =
    // tr(A*C) + tr(B*C)
    void bvisit(const MatrixMul &x)
    {
        auto scalar = x.get_scalar();
        if (eq(*scalar, *zero)) {
            trace_ = zero;
        } else {
            std::vector<vec_basic> args;
            // We first need to find out if there is MatrixAdd as a factor
            for (auto &factor : x.get_factors()) {
                if (is_a<const MatrixAdd>(*factor)) {
                    // Expand `MatrixAdd`
                    auto terms = down_cast<const MatrixAdd &>(*factor).get_args();
                    auto size_terms = terms.size();
                    auto size_args = args.size();
                    if (size_args == 0) {
                        args.reserve(size_terms);
                        for (std::size_t i = 0; i < size_terms; ++i)
                            args.push_back(vec_basic({terms[i]}));
                    } else {
                        args.reserve(size_args * size_terms);
                        for (std::size_t i = 1; i < size_terms; ++i)
                            for (std::size_t j = 0; j < size_args; ++j)
                                args.push_back(args[j]);
                        auto arg = args.begin();
                        for (std::size_t i = 0; i < size_terms; ++i)
                            for (std::size_t j = 0; j < size_args; ++j, ++arg)
                                arg->push_back(terms[i]);
                    }
                } else {
                    if (args.empty()) {
                        args.push_back(vec_basic({factor}));
                    } else {
                        for (auto &arg : args) arg.push_back(factor);
                    }
                }
            }
            vec_basic terms;
            for (auto &arg: args) {
                // For each `MatrixMul`, we rearrange the factors by using the
                // cyclic property
                auto min_factor = std::min_element(
                    arg.begin(), arg.end(), RCPBasicKeyLess()
                );
                if (min_factor == arg.begin()) {
                    terms.push_back(make_rcp<const Trace>(matrix_mul(arg)));
                } else {
                    auto factors = vec_basic(min_factor, arg.end());
                    factors.insert(factors.end(), arg.begin(), min_factor);
                    terms.push_back(make_rcp<const Trace>(matrix_mul(factors)));
                }
            }
            trace_ = add(terms);
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

#include <symengine/add.h>
#include <symengine/integer.h>
#include <symengine/rational.h>
#include <symengine/matrices/matrix_add.h>
#include <symengine/matrices/matrix_mul.h>
#include <symengine/matrices/zero_matrix.h>
#include <symengine/matrices/diagonal_matrix.h>
#include <symengine/matrices/immutable_dense_matrix.h>
#include <symengine/matrices/trace.h>

namespace SymEngine
{

hash_t MatrixAdd::__hash__() const
{
    hash_t seed = SYMENGINE_MATRIXADD;
    for (const auto &a : terms_) {
        hash_combine<Basic>(seed, *a);
    }
    return seed;
}

bool MatrixAdd::__eq__(const Basic &o) const
{
    if (is_a<MatrixAdd>(o)) {
        const MatrixAdd &other = down_cast<const MatrixAdd &>(o);
        // Because matrix addition satisfies commutative property and
        // associative property, so the order of terms does not matter.
        //
        // Converting `terms_` to `multiset_basic` is also problematic. Because
        // if a term contains `MatrixAdd` as one of its arguments, `__hash__()`
        // may be different if the `MatrixAdd` argument has `terms_` in
        // different order.
        if (terms_.size() != other.terms_.size()) return false;
        for (const auto& p: terms_) {
            auto found = false;
            for (const auto& q: other.terms_) {
                // Do not use `eq(*p, *q)` or `unified_eq(p, q)`
                if (p->__eq__(*q)) {
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
        return true;
    }
    return false;
}

int MatrixAdd::compare(const Basic &o) const
{
    SYMENGINE_ASSERT(is_a<MatrixAdd>(o));
    const MatrixAdd &other = down_cast<const MatrixAdd &>(o);
    // Not sure if we should compare `terms_` in a similar manner to that in
    // `__eq__` function
    return unified_compare(terms_, other.terms_);
}

bool MatrixAdd::is_canonical(const vec_basic &terms) const
{
    if (terms.size() < 2) {
        return false;
    }
    size_t num_diag = 0;
    size_t num_dense = 0;
    for (auto term : terms) {
        if (is_a<ZeroMatrix>(*term) || is_a<MatrixAdd>(*term)) {
            return false;
        } else if (is_a<DiagonalMatrix>(*term)) {
            num_diag++;
        } else if (is_a<ImmutableDenseMatrix>(*term)) {
            num_dense++;
        }
    }
    if (num_diag > 1 || num_dense > 1) {
        return false;
    }
    if (num_diag == 1 && num_dense == 1) {
        return false;
    }
    return true;
}

void check_matching_sizes(const vec_basic &vec)
{
    for (size_t i = 0; i < vec.size() - 1; i++) {
        auto first_size = size(down_cast<const MatrixExpr &>(*vec[i]));
        if (first_size.first.is_null()) {
            continue;
        }
        for (size_t j = 1; j < vec.size(); j++) {
            auto second_size = size(down_cast<const MatrixExpr &>(*vec[j]));
            if (second_size.first.is_null()) {
                continue;
            }
            auto rowdiff = sub(first_size.first, second_size.first);
            tribool rowmatch = is_zero(*rowdiff);
            if (is_false(rowmatch)) {
                throw DomainError("Matrix dimension mismatch");
            }
            auto coldiff = sub(first_size.second, second_size.second);
            tribool colmatch = is_zero(*coldiff);
            if (is_false(colmatch)) {
                throw DomainError("Matrix dimension mismatch");
            }
        }
    }
}

RCP<const MatrixExpr> matrix_add(const vec_basic &terms)
{
    if (terms.size() == 0) {
        throw DomainError("Empty sum of matrices");
    }
    if (terms.size() == 1) {
        return rcp_static_cast<const MatrixExpr>(terms[0]);
    }
    // extract nested MatrixAdd
    vec_basic expanded;
    for (auto &term : terms) {
        if (is_a<const MatrixAdd>(*term)) {
            auto container = down_cast<const MatrixAdd &>(*term).get_terms();
            expanded.insert(expanded.end(), container.begin(), container.end());
        } else {
            expanded.push_back(term);
        }
    }
    check_matching_sizes(expanded);
    // Coefficients for kept terms, which together with corresponding terms
    // will be converted to (matrix) multiplication
    vec_basic coef_keep;
    vec_basic keep;
    RCP<const DiagonalMatrix> diag;
    RCP<const ImmutableDenseMatrix> dense;
    RCP<const ZeroMatrix> zero;
    for (auto &term : expanded) {
        if (is_a<ZeroMatrix>(*term)) {
            zero = rcp_static_cast<const ZeroMatrix>(term);
        } else if (is_a<DiagonalMatrix>(*term)) {
            if (diag.is_null()) {
                diag = rcp_static_cast<const DiagonalMatrix>(term);
            } else {
                vec_basic container;
                for (size_t i = 0; i < diag->get_container().size(); i++) {
                    container.push_back(
                        add(diag->get_container()[i],
                            down_cast<const DiagonalMatrix &>(*term)
                                .get_container()[i]));
                }
                diag = make_rcp<const DiagonalMatrix>(container);
            }
        } else if (is_a<ImmutableDenseMatrix>(*term)) {
            if (dense.is_null()) {
                dense = rcp_static_cast<const ImmutableDenseMatrix>(term);
            } else {
                const vec_basic &vec1
                    = down_cast<const ImmutableDenseMatrix &>(*term)
                          .get_values();
                const vec_basic &vec2 = dense->get_values();
                vec_basic sum(vec1.size());
                for (size_t i = 0; i < vec1.size(); i++) {
                    sum[i] = add(vec1[i], vec2[i]);
                }
                dense = make_rcp<const ImmutableDenseMatrix>(
                    dense->nrows(), dense->ncols(), sum);
            }
        } else {
            // For matrix multiplication, we extract its scalar and factors
            RCP<const Basic> coef_term;
            RCP<const Basic> new_term;
            if (is_a<MatrixMul>(*term)) {
                coef_term = down_cast<const MatrixMul &>(*term).get_scalar();
                new_term = matrix_mul(down_cast<const MatrixMul &>(*term).get_factors());
            }
            else {
                coef_term = integer(1);
                new_term = term;
            }
            // Check if the term already exists
            bool not_exist = true;
            for (std::size_t i=0; i<coef_keep.size(); ++i) {
                if (eq(*keep[i], *new_term)) {
                    coef_keep[i] = add(coef_keep[i], coef_term);
                    not_exist = false;
                    break;
                }
            }
            if (not_exist) {
                coef_keep.push_back(coef_term);
                keep.push_back(new_term);
            }
        }
    }
    // Make product of each term and its coefficient
    for (std::size_t i=0; i<coef_keep.size(); ++i) {
        if (is_a_Number(*coef_keep[i])
            && rcp_static_cast<const Number>(coef_keep[i])->is_one())
        //if ((is_a<Integer>(*coef_keep[i])
        //    && down_cast<const Integer>(*coef_keep[i]).is_one())
        //    || (is_a<Rational>(*coef_keep[i])
        //    && down_cast<const Rational>(*coef_keep[i]).is_one()))
            continue;
        if (is_a<Trace>(*keep[i])) {
            keep[i] = mul(coef_keep[i], keep[i]);
        }
        else {
            keep[i] = matrix_mul({coef_keep[i], keep[i]});
        }
    }
    if (!diag.is_null()) {
        if (!dense.is_null()) {
            // Add diagonal with dense matrix
            auto vec = dense->get_values();
            vec_basic sum;
            for (size_t i = 0; i < dense->nrows(); i++) {
                for (size_t j = 0; j < dense->ncols(); j++) {
                    if (i == j) {
                        sum.push_back(add(dense->get(i, j), diag->get(i)));
                    } else {
                        sum.push_back(dense->get(i, j));
                    }
                }
            }
            dense = make_rcp<const ImmutableDenseMatrix>(dense->nrows(),
                                                         dense->ncols(), sum);
        } else {
            keep.push_back(diag);
        }
    }
    if (!dense.is_null()) {
        keep.push_back(dense);
    }
    if (keep.size() == 1) {
        return rcp_static_cast<const MatrixExpr>(keep[0]);
    }
    if (keep.size() == 0 && !zero.is_null()) {
        return zero;
    }
    return make_rcp<const MatrixAdd>(keep);
}

} // namespace SymEngine

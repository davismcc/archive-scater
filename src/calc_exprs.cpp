#include "scater.h"

template <typename T, class V, class M>
Rcpp::RObject calc_exprs_internal (M mat, Rcpp::RObject size_fac,  
        Rcpp::RObject prior_count, Rcpp::RObject log, 
        Rcpp::RObject sum, Rcpp::RObject subset) {

    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();
    if (size_fac.sexp_type()!=REALSXP || LENGTH(size_fac)!=ncells) { 
        throw std::runtime_error("length of 'size_fac' does not equal number of columns");
    }
    Rcpp::NumericVector sf(size_fac);
    Rcpp::IntegerVector subout=process_subset_vector(subset, mat, true);
    const size_t slen=subout.size();
    
    // Checking scalars.
    if (prior_count.sexp_type()!=REALSXP || LENGTH(prior_count)!=1) { 
        throw std::runtime_error("'prior_count' should be a numeric scalar");
    }
    const double prior=Rcpp::NumericVector(prior_count)[0];
    if (log.sexp_type()!=LGLSXP || LENGTH(log)!=1) {
        throw std::runtime_error("log specification should be a logical scalar"); 
    }
    const bool dolog=Rcpp::LogicalVector(log)[0];
    if (sum.sexp_type()!=LGLSXP || LENGTH(sum)!=1) {
        throw std::runtime_error("sum specification should be a sumical scalar"); 
    }
    const bool dosum=Rcpp::LogicalVector(sum)[0];

    V input(ngenes);
    Rcpp::NumericVector output(dosum ? slen : slen*ncells, 0);
    auto oIt=output.begin();

    // Computing normalized expression values for each cell, plus a prior.
    // We may or may not log-transform, and we may or may not sum across genes.
    for (size_t c=0; c<ncells; ++c) {
        mat->get_col(c, input.begin());
        const double& cursize=sf[c];

        for (auto sIt=subout.begin(); sIt!=subout.end(); ++sIt) {
            double tmp=double(input[*sIt])/cursize + prior;
            if (dosum) { 
                (*oIt)+=tmp;
            } else if (dolog) { 
                (*oIt)=std::log(tmp)/M_LN2;
                ++oIt;
            } else {
                (*oIt)=tmp;
                ++oIt;
            }
        }
        
        if (dosum) {
            ++oIt;
        }
    }

    // Cleaning up expected output.
    if (dosum) {
        if (dolog) {
            for (oIt=output.begin(); oIt!=output.end(); ++oIt) {
                (*oIt)=std::log(*oIt)/M_LN2;
            }
        }
    } else {
        output.attr("dim") = Rcpp::IntegerVector::create(slen, ncells);
    }
    return output;
}


SEXP calc_exprs(SEXP counts, SEXP size_fac, SEXP prior_count, SEXP log, SEXP sum, SEXP subset) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(counts);
        return calc_exprs_internal<int, Rcpp::IntegerVector>(mat.get(), size_fac, prior_count, log, sum, subset);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(counts);
        return calc_exprs_internal<double, Rcpp::NumericVector>(mat.get(), size_fac, prior_count, log, sum, subset);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

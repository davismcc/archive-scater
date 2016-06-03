#include "scater.h"

/*********************************************************
 This contains a number of small functions, written to improve 
 speed or memory efficiency over a native R implementation. 
 **********************************************************/

/* A function to get the proportion of reads in a subset. */

template <typename T>
SEXP colsum_subset_internal (const T* ptr, const matrix_info& MAT, SEXP subset) {
    if (!isInteger(subset)) { 
        throw std::runtime_error("subset vector must be an integer vector");
    }
    const int slen=LENGTH(subset);
    const int* sptr=INTEGER(subset);
    for (int s=0; s<slen; ++s) {
        if (sptr[s] < 1 || sptr[s] > MAT.nrow) { 
            throw std::runtime_error("subset indices out of range");
        }
    }

    SEXP output=PROTECT(allocVector(REALSXP, MAT.ncol));
    try {
        double* optr=REAL(output);
        
        // Summing across, using 1-indexed pointers.
        --ptr;
        int s;
        for (size_t c=0; c<MAT.ncol; ++c) {
            optr[c]=0;
            for (s=0; s<slen; ++s) {
                optr[c]+=ptr[sptr[s]];
            }
            ptr+=MAT.nrow;
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }
    UNPROTECT(1);
    return output;
}

SEXP colsum_subset(SEXP matrix, SEXP subset) try {
    matrix_info MAT=check_matrix(matrix);
    if (MAT.is_integer){
        return colsum_subset_internal<int>(MAT.iptr, MAT, subset);
    } else {
        return colsum_subset_internal<double>(MAT.dptr, MAT, subset);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}

/* A function to get the percentage of counts/expression 
   taken up by the top set of genes. */

template <typename T>
SEXP calc_top_features_internal(const T* ptr, const matrix_info& MAT, SEXP top) {
    if (!isInteger(top)) { 
        throw std::runtime_error("top specification must be an integer vector");
    }
    const int ntop=LENGTH(top);
    const int *tptr=INTEGER(top);
    for (size_t t=1; t<ntop; ++t) {
        if (tptr[t] < tptr[t-1]) { 
            throw std::runtime_error("numbers of top genes must be sorted"); 
        }
    }
    if (ntop && (tptr[0] < 1 || tptr[ntop-1] > MAT.nrow)) {
        throw std::runtime_error("number of top genes is out of index range");
    }
    
    SEXP output=PROTECT(allocMatrix(REALSXP, MAT.ncol, ntop));
    try {
        double** optrs=(double**)R_alloc(ntop, sizeof(double*));
        size_t t;
        if (ntop) {
            optrs[0]=REAL(output);
        }
        for (t=1; t<ntop; ++t) {
            optrs[t]=optrs[t-1]+MAT.ncol;
        }

        std::deque<T> values(MAT.nrow);
        size_t r, x;
        size_t target_index;
        double accumulated, total;

        for (size_t c=0; c<MAT.ncol; ++c) {
            total=0;
            for (r=0; r<MAT.nrow; ++r) {
                values[r]=ptr[r];
                total+=ptr[r];
            }
            
            // Sorting in descending order, and computing the accumulated total.
            std::sort(values.begin(), values.end(), std::greater<T>());
            x=0;
            accumulated=0;
            for (t=0; t<ntop; ++t) {
                target_index=size_t(tptr[t] - 1); // 0-based index.
                while (x<=target_index && x<MAT.nrow) {
                    accumulated+=values[x];
                    ++x;
                }
                optrs[t][c]=accumulated/total;
            }

            ptr+=MAT.nrow;
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
}

SEXP calc_top_features (SEXP matrix, SEXP top) try {
    matrix_info MAT=check_matrix(matrix);
    if (MAT.is_integer){
        return calc_top_features_internal<int>(MAT.iptr, MAT, top);
    } else {
        return calc_top_features_internal<double>(MAT.dptr, MAT, top);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}



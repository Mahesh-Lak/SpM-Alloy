#include "dia_matrix.h"

dia_matrix *dia_matrix_new(uint nrows, uint ncols, uint ndiags) {
    dia_matrix *dia = (dia_matrix *) calloc(1, sizeof(dia_matrix));
    dia->nrows = nrows;
    dia->ncols = ncols;
    dia->ndiags = ndiags;

    if (nrows > 0 && ndiags > 0) {
        dia->offsets = (int*) calloc(ndiags, sizeof(int));
        dia->vals = (real*) calloc(nrows * ndiags, sizeof(real));
    } else {
        dia->offsets = NULL;
        dia->vals = NULL;
    }

    return dia;
}

dia_matrix *dia_matrix_copy(dia_matrix const *in) {
    dia_matrix *dia = NULL;
    if (in != NULL) {
        dia = (dia_matrix *) calloc(1, sizeof(dia_matrix));
        dia->nrows = in->nrows;
        dia->ncols = in->ncols;
        dia->ndiags = in->ndiags;

        size_t size = dia->ndiags * sizeof(uint);
        dia->offsets = (int*) malloc(size);
        memcpy(dia->offsets, in->offsets, size);

        size = dia->ndiags * dia->nrows * sizeof(real);
        dia->vals = (real*) malloc(size);
        memcpy(dia->vals, in->vals, size);
    }

    return dia;
}

void dia_matrix_free(dia_matrix *dia) {
    if (dia != NULL) {
        if (dia->offsets != NULL) {
            free(dia->offsets);
        }
        if (dia->vals != NULL) {
            free(dia->vals);
        }
        free(dia);
    }
}

void dia_matrix_init(dia_matrix *dia, uint nnz, uint *rows, uint *cols, real *vals) {
    // Convert from COO to CSR...
    csr_matrix *csr = csr_matrix_new(dia->nrows, dia->ncols, nnz);
    csr_matrix_init(csr, rows, cols, vals);

    // Convert from CSR to DIA (inspector)...
    csr_dia_convert(csr, dia);
    csr_matrix_free(csr);
}

void dia_matrix_out(dia_matrix *dia, FILE *out) {
    uint i;
    fprintf(out, "dia_matrix: offsets: [ ");
    for (i = 0; i < dia->ndiags; i++) {
        fprintf(out, "%d ", dia->offsets[i]);
    }
    fprintf(out, "], vals: [ ");

    for (uint d = 0; d < dia->ndiags; d++) {
        fprintf(out, "[ ");
        for (uint i = 0; i < dia->nrows; i++) {
            fprintf(out, "%g ", dia->vals[offset2(i, d, dia->ndiags)]);
        }
        fprintf(out, "] ");
    }
    fprintf(out, "]\n");
}

void csr_dia_convert(csr_matrix const *csr, dia_matrix *dia) {
    uint i, j;
    uint N = csr->nrows;
    uint *col = csr->cols;
    uint *index = csr->rows;
    real *A = csr->vals;
    real *Aprime = NULL;

    /* Generate Inspector */
//    Dset = {[k'] | ∃j,k' = cols(j) − i ∧ rows(i) <= j < rows(i+1)}}
//    ND = count(Dset)
//    c = order(Dset)
//    Aprime = calloc(N ∗ ND, sizeof(real))
//    R_A→Aprime = {[j] → [i, d] | 0 <= d < ND ∧ ∃k', k' = cols(j) − i ∧ c(d) = k'}

#define order(i,j) k=col[(j)]-(i)+(N-1);\
                   if(c[k]==INT_MAX)c[k]=ND++

    // order
    uint d, k, ND = 0;
    uint maxdiag = N+N-1;
    uint *c = calloc(maxdiag, sizeof(int));
    for (uint i = 0; i < maxdiag; i++) c[i] = INT_MAX;

    for (i = 0; i < N; i++) {
        for (j = index[i]; j < index[i+1]; j++) {
            order(i,j);
        }
    }

    // count
//    ND=0;
//    for (i = 0; i < N; i++) {
//        for (j = index[i]; j < index[i+1]; j++) {
//            count(i,j);
//        }
//    }

#define copy(i,j) k=col[(j)]-(i)+(N-1);\
                  Aprime[offset2((i),c[k],ND)]=A[j]

    // copy
    Aprime = (real*) calloc(N * ND, sizeof(real));

    for (i = 0; i < N; i++) {
        for (j = index[i]; j < index[i+1]; j++) {
            copy(i,j);
        }
    }

#define invert(i,j) k=col[(j)]-(i);\
                    offsets[c[k+(N-1)]]=k

    // invert
    int *offsets = calloc(ND, sizeof(int));

    for (i = 0; i < N; i++) {
        for (j = index[i]; j < index[i+1]; j++) {
            invert(i,j);
        }
    }

    dia->nrows = csr->nrows;
    dia->ncols = csr->ncols;
    dia->ndiags = ND;
    dia->offsets = offsets;
    dia->vals = Aprime;
}

vector *dia_matvec_mul(dia_matrix *dia, vector *vec) {
    //assert(dia->ncols == vec->nvals);
    vector *out = vector_new(dia->nrows);

    real *A = dia->vals;
    real *x = vec->vals;
    real *y = out->vals;

    uint i, d, j, k;
    uint N = dia->nrows;
    uint ND = dia->ndiags;
    int *offsets = dia->offsets;

    for (i = 0; i < N; i++) {
        for (d = 0; d < ND; d++) {
            k = ND * i + d;
            //k = N * d + i;
            //if (A[k] != 0.0) {
                j = (i + offsets[d]) % N;   // Modding works to remove guard but not desirable!
                y[i] += A[k] * x[j];
                //fprintf(stderr, "i=%d,d=%d,offset=%d,j=%d,k=%d,A[%d][%d]=%g,x[%d]=%g,y[%d]=%g\n",
                //        i, d, offsets[d], j, k, i, d, A[k], j, x[j], i, y[i]);
            //}
        }
    }


    return out;
}

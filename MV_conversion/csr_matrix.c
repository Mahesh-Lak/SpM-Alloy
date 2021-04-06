#include "csr_matrix.h"

csr_matrix *csr_matrix_new(uint nrows, uint ncols, uint nnz) {
    csr_matrix *csr = (csr_matrix *) calloc(1, sizeof(csr_matrix));
    csr->nrows = nrows;
    csr->ncols = ncols;
    csr->nnz = nnz;

    if (nrows > 0 && nnz > 0) {
        csr->rows = (uint *) calloc(nrows + 1, sizeof(int));
        csr->cols = (uint *) calloc(nnz, sizeof(int));
        csr->vals = (real *) calloc(nnz, sizeof(real));
    } else {
        csr->rows = csr->cols = NULL;
        csr->vals = NULL;
    }

    return csr;
}

csr_matrix *csr_matrix_copy(csr_matrix const *in) {
    csr_matrix *csr = NULL;
    if (in != NULL) {
        csr = (csr_matrix *) calloc(1, sizeof(csr_matrix));
        csr->nrows = in->nrows;
        csr->ncols = in->ncols;
        csr->nnz = in->nnz;

        size_t size = (csr->nrows + 1) * sizeof(int);
        csr->rows = (uint *) malloc(size);
        memcpy(csr->rows, in->rows, size);

        size = csr->nnz * sizeof(int);
        csr->cols = (uint *) malloc(size);
        memcpy(csr->cols, in->cols, size);

        size = csr->nnz * sizeof(real);
        csr->vals = (real *) malloc(size);
        memcpy(csr->vals, in->vals, size);
    }

    return csr;
}

coo_matrix *csr_coo_convert(csr_matrix const *csr) {
    // known constants...
    uint N_R = csr->nrows;
    uint N_C = csr->ncols;
    uint NNZ = csr->nnz;

    coo_matrix *coo = calloc(1, sizeof(coo_matrix));
    coo->nrows = N_R;
    coo->ncols = N_C;
    coo->nnz = NNZ;

    coo->rows = calloc(NNZ, sizeof(uint));
    coo->cols = calloc(NNZ, sizeof(uint));
    coo->vals = calloc(NNZ, sizeof(real));

    // Merged order, invert, and copy...
    uint i, j;
    for (i = 0; i < N_R; i++) {
        for (j = csr->rows[i]; j < csr->rows[i + 1]; j++) {
            coo->rows[j] = i;
            coo->cols[j] = csr->cols[j];
            coo->vals[j] = csr->vals[j];
        }
    }

    return coo;
}

csr_matrix *csr_ell_convert(csr_matrix const *csr) {
    csr_matrix *ell = calloc(1, sizeof(csr_matrix));
    ell->nrows = csr->nrows;
    ell->ncols = csr->ncols;
    ell->nnz = csr->nnz;
    itype N = csr->nrows;
    itype M = csr->ncols;

    // 1) Find K (max # of nonzeros per row)
    itype K = 0;

    // 2) copy...
    itype K_est = csr->nnz / 2;
    itype* col_prime = calloc(K_est * N, sizeof(itype));
    real* A_prime = calloc(K_est * N, sizeof(real));

    for (itype i = 0; i < N; i++) {
        K = MAX(K, csr->rows[i+1] - csr->rows[i]);
        for (itype j = csr->rows[i]; j < csr->rows[i+1]; j++) {
            itype k = j - csr->rows[i];
            A_prime[offset2(k,i,N)] = csr->vals[j];
            col_prime[offset2(k,i,N)] = csr->cols[j];
        }
    }

    // 3) Resize arrays...
    col_prime = realloc(col_prime, K * N * sizeof(itype));\
    A_prime = realloc(A_prime, K * N * sizeof(real));

    ell->cols = col_prime;
    ell->vals = A_prime;
    ell->ncols = K;

    return ell;
}

void csr_matrix_init(csr_matrix *csr, uint *rows, uint *cols, real *vals) {
    // Matrix is already initialized, vals and cols are the same as COO input
    for (int i = 0; i < csr->nnz; ++i) {
        csr->cols[i] = cols[i];
        csr->vals[i] = vals[i];
        csr->rows[rows[i]] += 1;
    }
    csr->rows[csr->nrows] = csr->nnz;
    for (int i = csr->nrows - 1; i >= 0; --i) {
        csr->rows[i] = csr->rows[i + 1] - csr->rows[i];
    }
}

void csr_matrix_free(csr_matrix *csr) {
    if (csr != NULL && csr->nnz > 0) {
        free(csr->cols);
        free(csr->rows);
        free(csr->vals);
        free(csr);
    }
}

void csr_matrix_out(csr_matrix *csr, FILE *out) {
    uint i;
    fprintf(out, "csr_matrix: rows: [ ");
    for (i = 0; i <= csr->nrows; i++) {
        fprintf(out, "%d ", csr->rows[i]);
    }
    fprintf(out, "], cols: [ ");
    for (i = 0; i < csr->nnz; i++) {
        fprintf(out, "%d ", csr->cols[i]);
    }
    fprintf(out, "], vals: [ ");
    for (i = 0; i < csr->nnz; i++) {
        fprintf(out, "%g ", csr->vals[i]);
    }
    fprintf(out, "]\n");
}

// CSR-SpMV
vector *csr_matvec_mul(csr_matrix *csr, vector *vec) {
    //assert(csr->ncols == vec->nvals);
    vector *out = vector_new(csr->nrows);

    uint N = csr->nrows;
    real *A = csr->vals;
    real *x = vec->vals;
    real *y = out->vals;

    for (uint i = 0; i < N; i++) {
        for (uint j = csr->rows[i]; j < csr->rows[i+1]; j++) {
            uint k = csr->cols[j];
            y[i] += A[j] * x[k];
            //fprintf(stderr, "row=%d,col=%d,ndx=%d,A[%d]=%g,x[%d]=%g,y[%d]=%g\n", i, k, j, j, A[j], k, x[k], i, y[i]);
        }
    }

    return out;
}

// ELL-SpMV
vector *ell_matvec_mul(csr_matrix *ell, vector *vec) {
    vector *out = vector_new(csr->nrows);

    uint N = ell->nrows;
    uint M = ell->ncols;
    uint *cols = ell->cols;
    real *A = ell->vals;
    real *x = vec->vals;
    real *y = out->vals;

    for (uint i = 0; i < N; i++) {
        for (itype j = 0; j < M; j++) {
            uint k = cols[offset2(j,i,N)];
            y[i] += A[offset2(j,i,N)] * x[k];
        }
    }

    return out;
}
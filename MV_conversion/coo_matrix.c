#include "coo_matrix.h"

coo_matrix *coo_matrix_new(uint nrows, uint ncols, uint nnz) {
    coo_matrix *coo = (coo_matrix *) calloc(1, sizeof(coo_matrix));
    coo->nrows = nrows;
    coo->ncols = ncols;
    coo->nnz = nnz;

    if (nrows > 0 && nnz > 0) {
        coo->rows = (uint *) calloc(nnz, sizeof(int));
        coo->cols = (uint *) calloc(nnz, sizeof(int));
        coo->vals = (real *) calloc(nnz, sizeof(real));
    } else {
        coo->rows = coo->cols = NULL;
        coo->vals = NULL;
    }

    return coo;
}

coo_matrix *coo_matrix_copy(coo_matrix const *in) {
    coo_matrix *coo = NULL;
    if (in != NULL) {
        coo = (coo_matrix *) calloc(1, sizeof(coo_matrix));
        coo->nrows = in->nrows;
        coo->ncols = in->ncols;
        coo->nnz = in->nnz;

        size_t size = (coo->nrows + 1) * sizeof(int);
        coo->rows = (uint *) malloc(size);
        memcpy(coo->rows, in->rows, size);

        size = coo->nnz * sizeof(int);
        coo->cols = (uint *) malloc(size);
        memcpy(coo->cols, in->cols, size);

        size = coo->nnz * sizeof(real);
        coo->vals = (real *) malloc(size);
        memcpy(coo->vals, in->vals, size);
    }

    return coo;
}

void coo_matrix_init(coo_matrix *coo, uint *rows, uint *cols, real *vals) {
    // Matrix is already initialized, vals and cols are the same as COO input
    for (uint i = 0; i < coo->nnz; ++i) {
        coo->rows[i] = rows[i];
        coo->cols[i] = cols[i];
        coo->vals[i] = vals[i];
    }
}

void coo_matrix_free(coo_matrix *coo) {
    if (coo != NULL && coo->nnz > 0) {
        free(coo->cols);
        free(coo->rows);
        free(coo->vals);
        free(coo);
    }
}

void coo_matrix_out(coo_matrix *coo, FILE *out) {
    uint i;
    fprintf(out, "coo_matrix: rows: [ ");
    for (i = 0; i < coo->nnz; i++) {
        fprintf(out, "%d ", coo->rows[i]);
    }
    fprintf(out, "], cols: [ ");
    for (i = 0; i < coo->nnz; i++) {
        fprintf(out, "%d ", coo->cols[i]);
    }
    fprintf(out, "], vals: [ ");
    for (i = 0; i < coo->nnz; i++) {
        fprintf(out, "%g ", coo->vals[i]);
    }
    fprintf(out, "]\n");
}

vector *coo_matvec_mul(coo_matrix *coo, vector *vec) {
    vector *out = vector_new(coo->nrows);
    real *x = vec->vals;
    real *y = out->vals;

    for (uint i = 0; i < coo->nnz; ++i) {
        y[coo->rows[i]] += coo->vals[i] * x[coo->cols[i]];
    }

    return out;
}

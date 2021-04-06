#ifndef COO_MATRIX_H
#define COO_MATRIX_H

#include "matrix.h"
#include "vector.h"

typedef struct _coo_matrix {
    uint nrows;
    uint ncols;
    uint nnz;
    uint *rows;
    uint *cols;
    real *vals;
} coo_matrix;

coo_matrix *coo_matrix_new(uint nrows, uint ncols, uint nnz);
coo_matrix *coo_matrix_copy(coo_matrix const *in);
void coo_matrix_init(coo_matrix *coo, uint *rows, uint *cols, real *vals);
void coo_matrix_free(coo_matrix *coo);
void coo_matrix_out(coo_matrix *coo, FILE *out);
vector *coo_matvec_mul(coo_matrix *coo, vector *vec);

#endif // COO_MATRIX_H

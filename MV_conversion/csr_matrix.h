#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H

#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"
#include "vector.h"
#include "coo_matrix.h"

typedef struct _csr_matrix {
    uint nrows;
    uint ncols;
    uint nnz;
    uint *rows;
    uint *cols;
    real *vals;
} csr_matrix;

csr_matrix *csr_matrix_new(uint nrows, uint ncols, uint nnz);
csr_matrix *csr_matrix_copy(csr_matrix const *in);
coo_matrix *csr_coo_convert(csr_matrix const *csr);
csr_matrix *csr_ell_convert(csr_matrix const *csr);
void csr_matrix_init(csr_matrix *csr, uint *rows, uint *cols, real *vals);
void csr_matrix_free(csr_matrix *csr);
void csr_matrix_out(csr_matrix *csr, FILE *out);
vector *csr_matvec_mul(csr_matrix *csr, vector *vec);
vector *ell_matvec_mul(csr_matrix *csr, vector *vec);

#endif // CSR_MATRIX_H

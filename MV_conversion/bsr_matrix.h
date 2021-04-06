#ifndef BSR_MATRIX_H
#define BSR_MATRIX_H

#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include "csr_matrix.h"
#include "matrix.h"
#include "vector.h"

typedef struct _bsr_matrix {
    uint nrows;
    uint ncols;
    uint nb;        // # NZ blocks
    uint nd;        // Block dimensionality
    uint bsizer;    // Block size R
    uint bsizec;    // Block size C
    uint *rows;     // block-row
    uint *cols;     // block-col
    real *vals;
} bsr_matrix;

bsr_matrix *bsr_matrix_new(uint nrows, uint ncols, uint bsizer, uint bsizec);
bsr_matrix *bsr_matrix_copy(bsr_matrix const *in);
void bsr_matrix_free(bsr_matrix *bsr);
void bsr_matrix_init(bsr_matrix *bsr, uint nnz, uint *rows, uint *cols, real *vals);
void bsr_matrix_out(bsr_matrix *bsr, FILE *out);
void csr_bsr_convert(csr_matrix const *csr, bsr_matrix *bsr);
vector *bsr_matvec_mul(bsr_matrix *bsr, vector *vec);

#endif // BSR_MATRIX_H

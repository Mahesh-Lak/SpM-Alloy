//
// DIA is useful if the fill ratio defined below is greater than one:
//
// DIA Fill Ratio = (Ndiags * Nrows) / NNZ
//

#ifndef DIA_MATRIX_H
#define DIA_MATRIX_H

#include <limits.h>
#include "matrix.h"
#include "csr_matrix.h"
#include "vector.h"

typedef struct _dia_matrix {
    uint nrows;
    uint ncols;
    uint ndiags;
    int *offsets;
    real *vals;
} dia_matrix;

dia_matrix *dia_matrix_new(uint nrows, uint ncols, uint ndiags);
dia_matrix *dia_matrix_copy(dia_matrix const *in);
void dia_matrix_free(dia_matrix *dia);
void dia_matrix_init(dia_matrix *dia, uint nnz, uint *rows, uint *cols, real *vals);
void dia_matrix_out(dia_matrix *dia, FILE *out);
void csr_dia_convert(csr_matrix const *csr, dia_matrix *dia);
vector *dia_matvec_mul(dia_matrix *dia, vector *vec);

#endif // DIA_MATRIX_H

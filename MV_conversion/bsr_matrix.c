#include "bsr_matrix.h"

bsr_matrix *bsr_matrix_new(uint nrows, uint ncols, uint bsizer, uint bsizec) {
    bsr_matrix *bsr = (bsr_matrix *) calloc(1, sizeof(bsr_matrix));
    bsr->nrows = nrows;
    bsr->ncols = ncols;
    bsr->bsizer = bsizer;
    bsr->bsizec = bsizec;
    bsr->nd = nrows / bsizer;

    // Find a suitable block size if none defined...
    uint size;
    for (size = bsr->bsizer; bsr->nrows % size != 0; ++size);
    bsr->bsizer = size;

    for (size = bsr->bsizec; bsr->ncols % size != 0; ++size);
    bsr->bsizec = size;

    // These will be allocated later...
    bsr->nb = 0;
    bsr->rows = NULL;
    bsr->cols = NULL;
    bsr->vals = NULL;

    return bsr;
}

bsr_matrix *bsr_matrix_copy(bsr_matrix const *in) {
    bsr_matrix *out = (bsr_matrix *) calloc(1, sizeof(bsr_matrix));
    out->nrows = in->nrows;
    out->ncols = in->ncols;
    out->bsizer = in->bsizer;
    out->bsizec = in->bsizec;
    out->nb = in->nb;

    if (in->rows != NULL && in->cols != NULL) {
        out->rows = (uint *) calloc(in->nb + 1, sizeof(uint));
        out->cols = (uint *) calloc(in->nb, sizeof(uint));
        out->vals = (real *) calloc(in->nb, sizeof(real));

        uint bsize = in->bsizer * in->bsizec;
        for (uint i = 0; i < in->nb; ++i) {
            out->rows[i] = in->rows[i];
            out->cols[i] = in->cols[i];
            for (uint j = 0; j < bsize; j++) {
                uint k = offset2(i, j, bsize);
                out->vals[k] = in->vals[k];
            }
        }
        out->rows[out->nb] = in->nb;
    }

    return out;
}

/// This version is using the I_{bcsr} iteration space, a la the PLDI'15 paper.
/// \param csr Input CSR matrix
/// \param bsr Output BCSR matrix
void csr_bsr_convert(csr_matrix const *csr, bsr_matrix *bsr) {
//void csr_bsr_convert_dev(csr_matrix const *csr, bsr_matrix *bsr) {
    uint N = csr->nrows;
    uint R = bsr->bsizer;
    uint C = bsr->bsizec;
    real *A = csr->vals;
    uint *index = csr->rows;
    uint *col = csr->cols;
    uint bb, ii, kk, ri, ck, i, j, k;

    // 1) P-C relationship, statement from one block READS data from another, we can shift and fuse.
    // 2) If iteration space relies on data produced by a node, then ALL of the data needs to be produced first.

//#define makesetd(ii,kk) bset[(ii)*(N/C)+(kk)]|=1

    // makeset
    uint nblocks = (N/R)*(N/C);      // max # blocks

#define countd(ii,kk) bb=(ii)*(N/C)+(kk);\
                      if(!marked[bb]){NB+=1;marked[bb]=1;}
//if(!marked[bb]){NB+=bset[bb];marked[bb]=1;}

    // count
    uint NB = 0;

    // Need to mark blocks so as not to count the same one twice.
    uint *marked = calloc(nblocks, sizeof(uint));

    for (ii = 0; ii < N/R; ii++) {
        for (ri = 0; ri < R; ri++) {
            for (j = index[ii*R+ri]; j < index[ii*R+ri+1]; j++) {
                countd(ii,(col[(j)]/C));
            }
        }
    }

    free(marked);
    bsr->nb = NB;
    bsr->nd = N/R;

// b+=bset[bb];
#define orderd(ii,kk) bb=(ii)*(N/C)+(kk);\
                      if(!ordered[bb]){\
                        b+=1;\
                        ordered[bb]=b;\
                      }

    // order
    uint *ordered = calloc(nblocks, sizeof(uint));
    uint b = 0;

    for (ii = 0; ii < N/R; ii++) {
        for (ri = 0; ri < R; ri++) {
            for (j = index[ii*R+ri]; j < index[ii*R+ri+1]; j++) {
                orderd(ii, (col[(j)]/C));
            }
        }
    }

#define copyd(i,j) kk=col[(j)]/C;\
                   bb=(ii)*(N/C)+(kk);\
                   b=ordered[bb]-1;\
                   ck=col[(j)]-kk*C;\
                   k=offset3(b,(ri),(ck),R,C);\
                   A_prime[k]=A[(j)]

    // copy
    bsr->vals = calloc(NB*R*C, sizeof(real));
    real *A_prime = bsr->vals;

    for (ii = 0; ii < N/R; ii++) {
        for (ri = 0; ri < R; ri++) {
            for (j = index[ii*R+ri]; j < index[ii*R+ri+1]; j++) {
                copyd((ii*R+ri), j);
            }
        }
    }

//    b_index[(ii)+1]=ordered[bb];\
//    ordered[bb]=0;
#define invertd(ii,kk) bb=(ii)*(N/C)+(kk);\
                       if(ordered[bb]){\
                         b=ordered[bb]-1;\
                         b_row[b]=(ii);\
                         b_col[b]=(kk);\
                       }

    // invert
    uint *b_col = calloc(NB, sizeof(uint));
    uint *b_row = calloc(NB, sizeof(uint));
    bsr->cols = b_col;
    uint *b_index = calloc(N/R+1, sizeof(uint));
    bsr->rows = b_index; //block_row;

    for (ii = 0; ii < N/R; ii++) {
        for (ri = 0; ri < R; ri++) {
            for (j = index[ii*R+ri]; j < index[ii*R+ri+1]; j++) {
                invertd(ii, (col[(j)]/C));
            }
        }
    }

    for (b = 0; b < NB; b++) {
        b_index[b_row[b] + 1] = b + 1;
    }

    // Cleanup...
    //free(bset);
    free(ordered);
}

void bsr_matrix_free(bsr_matrix *bsr) {
    if (bsr != NULL) {
        if (bsr->rows != NULL) {
            free(bsr->rows);
        }
        if (bsr->cols != NULL) {
            free(bsr->cols);
        }
        if (bsr->vals != NULL) {
            free(bsr->vals);
        }
        free(bsr);
    }
}

void bsr_matrix_init(bsr_matrix *bsr, uint nnz, uint *rows, uint *cols, real *vals) {
    // Convert from COO to CSR...
    csr_matrix *csr = csr_matrix_new(bsr->nrows, bsr->ncols, nnz);
    csr_matrix_init(csr, rows, cols, vals);

    // Convert from CSR to BSR (inspector)...
    csr_bsr_convert(csr, bsr);
    csr_matrix_free(csr);
}

void bsr_matrix_out(bsr_matrix *bsr, FILE *out) {
    fprintf(out, "bsr_matrix: rows = [ ");
    for (uint i = 0; i <= bsr->nd; i++) {
        fprintf(out, "%d ", bsr->rows[i]);
    }
    fprintf(out, "], ");

    fprintf(out, "cols = [ ");
    for (uint i = 0; i < bsr->nb; i++) {
        fprintf(out, "%d ", bsr->cols[i]);
    }
    fprintf(out, "], ");

    fprintf(out, "vals = [ ");
    for (uint i = 0; i < bsr->nb; i++) {
        uint bsize = bsr->bsizer * bsr->bsizec;
        fprintf(out, "[ ");
        for (uint j = 0; j < bsize; j++) {
            fprintf(out, "%g ", bsr->vals[offset2(i, j, bsize)]);
        }
        fprintf(out, "] ");
    }
    fprintf(out, "]\n");
}

void bsr_matrix_write(bsr_matrix *bsr, FILE *out) {
    uint datasize = bsr->nb * bsr->bsizer * bsr->bsizec;
    fprintf(out,"%d  %d  %d  %d  %d  %d  %d\n", bsr->nrows, bsr->ncols, bsr->bsizer, bsr->bsizec, bsr->nd, bsr->nb, datasize);

    for (uint i = 0; i <= bsr->nb; i++) {
        fprintf(out, "%d\n", bsr->rows[i]);
    }

    for (uint i = 0; i < bsr->nb; i++) {
        fprintf(out, "%d\n", bsr->cols[i]);
    }

    for (uint i = 0; i < datasize; i++) {
        fprintf(out, "%.6lf\n", bsr->vals[i]);
    }
}

vector *bsr_matvec_mul(bsr_matrix *bsr, vector *vec) {
    //assert(bsr->ncols == vec->nvals);
    vector *out = vector_new(bsr->nrows);

    uint N = bsr->nrows;
    uint R = bsr->bsizer;
    uint C = bsr->bsizec;

    uint *bcol = bsr->cols;
    uint *brow = bsr->rows;

    real *a = bsr->vals;
    real *x = vec->vals;
    real *y = out->vals;

    // PLDI'15 Executor:
    for (uint ii = 0; ii < N/R; ii++) {
        for (uint jj = brow[ii]; jj < brow[ii+1]; jj++) {  // brow <=> bsr->rows <=> offset_index <=> bsr_matrix.indptr
            uint kk = bcol[jj];                            // bcol <=> bsr->cols <=> explicit_index <=> bsr_matrix.indices
            for (uint ri = 0; ri < R; ri++) {
                for (uint ck = 0; ck < C; ck++) {
                    uint i = ii * R + ri;
                    uint k = kk * C + ck;
                    uint m = offset3(jj, ri, ck, R, C);
                    //fprintf(stderr, "ii=%d,kk=%d,jj=%d,ri=%d,ck=%d,i=%d,k=%d,m=%d,v=%g\n", ii, kk, jj, ri, ck, i, j, k, bsr->vals[k]);
                    y[i] += a[m] * x[k];
                }
            }
        }
    }

    // IMPACT'18 Executor:
//    for (uint b = 0; b < bsr->nb; b++) {
//        uint ii = bsr->rows[b];
//        uint kk = bsr->cols[b];
//        for (uint ri = 0; ri < R; ri++) {
//            for (uint ck = 0; ck < C; ck++) {
//                uint i = ii * R + ri;
//                uint j = kk * C + ck;
//                uint k = offset3(ii, ri, ck, R, C);
//                //fprintf(stderr, "ii=%d,kk=%d,ri=%d,ck=%d,i=%d,j=%d,k=%d,v=%g\n", ii, kk, ri, ck, i, j, k, bsr->vals[k]);
//                y[i] += bsr->vals[k] * x[j];
//            }
//        }
//    }

    return out;
}

/*
 * Specification of matrix-vector multiplication for dense
 * matrices and CSR sparse format, showing refinement
 * using CSR abstraction function
 * 
 * Authors: Tristan Dyer, Alper Altunutas, John Baugh
 * Date: August 16, 2019
 * Alloy Analyzer 5.0.0.201804081720
 *
 * For a detailed description, see:
 * 
 *  Bounded Verification of Sparse Matrix Computations
 *  Third International Workshop on Software Correctness for HPC Applications
 *  (submitted)
 *   
 *    => Section V.B. Matrix-Vector Multiplication
 */

open structure
open behavior

-- a sum of products
sig SumProd {
  vals: Int->lone Value->Value
}

-- dense matrix-vector multiplication
pred MVM [m: Matrix, x: seq Value, b: seq SumProd] {
  m.rows = #b
  m.cols = #x
  all i: range[m.rows] |
    b[i].vals = { j: Int, p, q: Value-Zero |
      j in range[m.cols] and p = m.vals[i][j] and q = x[j] }
}

-- CSR matrix-vector multiplication
pred MVM [c: CSR, x: seq Value, b: seq SumProd] {
  c.rows = #b
  c.cols = #x
  all i: range[c.rows] |
    b[i].vals = { j: Int, p, q: Value-Zero |
      some k: range[c.IA[i], c.IA[i.add[1]]] |
        j = c.JA[k] and p = c.A[k] and q = x[c.JA[k]] }
}

/*
 * Below we perform two checks:
 * 
 *   1. alphaValid checks that the abstraction function maps
 *      valid CSR representations to valid dense representations
 *
 *   2. mvmRefines checks that the CSR MVM method is a
 *      valid refinement of the dense MVM method
 */  

assert alphaValid {
  all c: CSR, m: Matrix |
    I[c] and alpha[c, m] => I[m]
}

assert mvmRefines {
  all m: Matrix, c: CSR, x: seq Value, b: seq SumProd |
    I[c] and alpha[c, m] and MVM[c, x, b]
      => MVM[m, x, b]
}

check alphaValid for 1 Matrix, 1 CSR, 0 ELL, 3 Value, 0 SumProd
check mvmRefines for 1 Matrix, 1 CSR, 0 ELL, 3 Value, 6 SumProd, 6 seq

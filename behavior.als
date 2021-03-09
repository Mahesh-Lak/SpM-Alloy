module behavior

/*
 * Specification of the ELL and CSR abstraction functions, dense and 
 * ELL update operations, and correctness checks
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
 *    => Section IV. Matrix Behavior
 */

open structure

-- ELL abstraction function
pred alpha [e: ELL, m: Matrix] {
  m.rows = e.rows
  m.cols = e.cols
  m.vals = {
    i: range[e.rows], j: range[e.cols], v: Value |
      let k = e.jcoef[i].j |
        some k => v = e.coef[i][k] else v = Zero
  }
}

-- CSR abstraction function
pred alpha [c: CSR, m: Matrix] {
  m.rows = c.rows
  m.cols = c.cols
  m.vals = {
    i: range[c.rows], j: range[c.cols], v: Value |
      let k = { k: range[c.IA[i], c.IA[add[i, 1]]] | c.JA[k] = j } |
        one k => v = c.A[k] else v = Zero
  }
}

-- update single value of dense matrix
pred update [m, m': Matrix, i, j: Int, v: Value] {
  i->j in indices[m]
  let u = m.vals[i][j] |
    m'.vals = m.vals - i->j->u + i->j->v
  m'.sameDimensions[m]
}

-- update single value of ELL matrix
pred update [e, e': ELL, i, j: Int, v: Value] {
  i->j in indices[e]
  v = Zero => toZero[e, e', i, j] 
    else toNonZero[e, e', i, j, v]
  e'.sameDimensions[e]
}

-- set single value of ELL to be zero
pred toZero [e, e': ELL, i, j: Int] {
  let k = e.jcoef[i].j {
    e'.jcoef = e.jcoef.setAt[i, k, -1]
    e'.coef = e.coef.setAt[i, k, Zero]
  }
}

-- set single value of ELL matrix to be nonzero
pred toNonZero [e, e': ELL, i, j: Int, v: Value] {
  some e.jcoef[i].j => {
   let k = e.jcoef[i].j |
    e'.coef = e.coef.setAt[i, k, v] and
    e'.jcoef = e.jcoef
  } else {
    one k: Int |
      e.jcoef[i][k] = -1 and
      e'.coef = e.coef.setAt[i, k, v] and
      e'.jcoef = e.jcoef.setAt[i, k, j]
  }
}

-- set of valid indices for abstract matrix
fun indices [m: Matrix]: Int->Int {
  range[m.rows]->range[m.cols]
}

-- set of valid indices for ELL matrix
fun indices [e: ELL]: Int->Int {
  range[e.rows]->range[e.cols]
}

-- two abstract matrices have the same dimensions
pred sameDimensions [m, m': Matrix] {
  m.rows = m'.rows
  m.cols = m'.cols
}

-- two ELL matrices have the same dimensions and maxnz
pred sameDimensions [e, e': ELL] {
  e.rows = e'.rows
  e.cols = e'.cols
  e.maxnz = e'.maxnz
}

-- two dense matrices are equivalent
pred equivalent [m, m': Matrix] {
  m.rows = m'.rows
  m.cols = m'.cols
  m.vals = m'.vals
}

-- generate duplicate set with single value changed
fun setAt [s: Int->Int->univ, m, n: Int, v: univ]: Int->Int->univ {
  let c = s[m][n] | s - m->n->c + m->n->v
}

/*
 * Below we perform four checks:
 * 
 *   1. alphaValid checks that the abstraction function maps
 *      valid ELL representations to valid dense representations
 *
 *   2. preservesRep checks that the ELL update method preserves
 *      the ELL concrete invariant
 *
 *   3. updateRefines checks that the ELL update method is a
 *      valid refinement of the dense update method
 *
 *   4. updateDeterministic checks that the dense update method
 *      is deterministic in order to help ensure our spec is correct
 */      

assert alphaValid {
  all e: ELL, m: Matrix |
    I[e] and alpha[e, m] => I[m]
}

assert preservesRep {
  all e, e': ELL, i, j: Int, v: Value |
    I[e] and update[e, e', i, j, v] => I[e']
}

assert updateRefines {
  all e, e': ELL, m, m': Matrix, i, j: Int, v: Value |
    I[e] and alpha[e, m] and alpha[e', m'] and
      update[e, e', i, j, v] => update[m, m', i, j, v]
}

assert updateDeterministic {
  all m, m', m'': Matrix, i, j: Int, v: Value |
    I[m] and update[m, m', i, j, v] and update[m, m'', i, j, v]
      => equivalent[m', m'']
}

check alphaValid for 1 Matrix, 0 CSR, 1 ELL, 5 Value, 4 Int
check preservesRep for 0 Matrix, 0 CSR, 2 ELL, 5 Value, 4 Int
check updateRefines for 2 Matrix, 0 CSR, 2 ELL, 5 Value, 4 Int
check updateDeterministic for 2 Matrix, 0 CSR, 0 ELL, 5 Value, 4 Int

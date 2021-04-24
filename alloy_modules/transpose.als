/*
 * Specification of the CSR transpose algorithm found
 * in the SPARSKIT library
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
 *    => Section V.D. CSR Transpose
 */

open structure
open behavior

-- Dense matrix transpose
pred transpose [m, m': Matrix] {
  m'.rows = m.cols
  m'.cols = m.rows
  m'.vals = { j, i: Int, v: Value | i->j->v in m.vals } //??use this to represent relation?
}

-- CSR matrix transpose
pred transpose [c, c': CSR] {
  c'.rows = c.cols
  c'.cols = c.rows
  some iao, iao', iao'', iao''': seq Int {   //seq Int: Int -> Int
    countcols[c, iao]
    setptrs[iao, iao']
    copyvals[c, c', iao', iao'']
    shift[iao'', iao''']
    c'.IA = iao'''
    #c'.JA = c'.IA.end
    #c'.A = c'.IA.end
  }
}

-- Phase 1 of the transpose algorithm:
-- Count the number of non-zero values in each column of initial matrix
pred countcols[c: CSR, iao: seq Int] {
  #iao = c.cols.add[1]
  iao[0] = 0
  all j: range[c.cols] |
    iao[j.add[1]] = #c.JA.j //c.JA.j: all elements in JA whose value is j
}

-- Phase 2 of the transpose algorithm:
-- Assemble output IA array using number of non-zeros in each column
pred setptrs [iao, iao': seq Int] {
  #iao' = #iao
  iao'[0] = 0
  all i: iao.butlast.inds |   //butlast: (seq Int) except the last element; inds:index set
    let i' = i.add[1] |
      iao'[i'] = iao'[i].add[iao[i']]
}
//? This is endding index not starting index


-- Phase 3 of the transpose algorithm:
-- Copy values into output matrix.  Has side effect of right-shifting output IA
pred copyvals [c, c': CSR, iao_in, iao_out: seq Int] {
  #iao_out = #iao_in
  let loopvars = { i: range[c.rows], k: range[c.IA[i], c.IA[i.add[1]]] } | //loopvars: i -> k
    some iter: Int->Int->Int, iao: Int->Int->Int {
      table[loopvars, iter]
      iao[0] = iao_in
      #iao.Int.Int = add[#iter, 1]
      all i: range[c.rows] {
        all k: range[c.IA[i], c.IA[i.add[1]]] {
          let t = iter[i][k],
              t' = t.add[1],
              j = c.JA[k],
              nxt = iao[t][j] {
            iao[t'] = iao[t] ++ j->nxt.add[1]
            c'.JA[nxt] = i
            c'.A[nxt] = c.A[k]
          }
        }
      }
      iao_out = iao[c.IA.end]
    }
}

-- Phase 4 of the transpose algorithm:
-- Left-shift the output IA to undo side effect from step 3
pred shift [iao, iao': seq Int] {
  #iao' = #iao
  iao'[0] = 0
  all i: range[1, #iao] |
    iao'[i] = iao[i.sub[1]]
}

-- establish iteration table by adding time column to end of relation r
pred table [r: Int->Int, iter: Int->Int->Int] {
  #iter = #r                -- same size as r
  iter.Int = r              -- make left side = r
  c3[iter] = range[#r]      -- make right side = time stamps
  all i, i': r.Int,
      j, j': Int.r,
      t, t': range[#r] {
    i->j->t in iter and i'->j'->t' in iter and t < t' =>
      i <= i' and (i = i' => j < j')
  }
}

-- third column of a ternary relation
fun c3 [r: Int->Int->Int]: Int {
  Int.(Int.r)
}

/*
 * Below we perform four checks:
 * 
 *   1. alphaValid checks that the abstraction function maps
 *      valid CSR representations to valid dense representations
 *
 *   2. preservesRep checks that the CSR transpose method preserves
 *      the CSR concrete invariant
 *
 *   3. transposeRefines checks that the CSR transpose method is a
 *      valid refinement of the dense transpose method
 *
 *   4. transposeDeterministic checks that the dense transpose method
 *      is deterministic in order to help ensure our spec is correct
 */ 

assert alphaValid {
  all c: CSR, m: Matrix |
    I[c] and alpha[c, m] => I[m]
}

assert preservesRep {
  all c, c': CSR |
    I[c] and transpose[c, c'] => I[c']
}

assert transposeRefines {
  all c, c': CSR, m, m': Matrix |
    I[c] and transpose[c, c'] and alpha[c, m] and alpha[c', m']
      => transpose[m, m']
}

assert transposeDeterministic {
  all m, m', m'': Matrix |
    I[m] and transpose[m, m'] and transpose[m, m'']
      => equivalent[m', m'']
}

check alphaValid for 1 Matrix, 0 ELL, 1 CSR, 3 Value
check preservesRep for 0 Matrix, 0 ELL, 1 CSR, 3 Value, 7 seq
check transposeRefines for 2 Matrix, 0 ELL, 2 CSR, 3 Value, 7 seq
check transposeDeterministic for 3 Matrix, 0 ELL, 0 CSR, 3 Value

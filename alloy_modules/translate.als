/*
 * Specification of the ELL to CSR translation algorithm
 * found in the SPARSKIT library
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
 *    => Section V.C. ELL to CSR translation
 */

open structure
open behavior

-- translate ELL to CSR
pred ellcsr [e: ELL, c: CSR] {
  c.rows = e.rows
  c.cols = e.cols
  c.IA[0] = 0 //c.IA.0=0
  #c.IA = c.rows.add[1] // maximum size of IA
				   // IA is the relation inside CSR 
				   // `add` is the built-in function in int
  some iter: Int->Int->Int, kpos: Int->Int {
    table[range[e.rows]->range[e.maxnz], iter]
    kpos[0] = 0
    #kpos = add[#iter, 1] //maximum size of kpos
    #c.A = kpos.end
    #c.JA = kpos.end
    all i: range[e.rows] {
      all k: range[e.maxnz] |
        let t = iter[i][k], t' = t.add[1] |   //A let statement acts as a macro replacing the right-side of the assignment wherever the left-side of the assignment appears
          e.coef[i][k] != Zero => {
            c.A[kpos[t]] = e.coef[i][k]   //c, e are all known input
            c.JA[kpos[t]] = e.jcoef[i][k]
            kpos[t'] = kpos[t].add[1]
          } else 
            kpos[t'] = kpos[t]   // If not equal, we skip this index without
                                       // increase A to next time stamp
      c.IA[i.add[1]] = kpos[end[iter, i].add[1]]
    }
  }
}

-- establish iteration table by adding time column to end of relation r
pred table [r: Int->Int, iter: Int->Int->Int] {
  #iter = #r                -- same size as r
  iter.Int = r              -- make left side = r // the last element is int, join the last element and Int obtain left, i.e., Int -> Int
  c3[iter] = range[#r]      -- make right side = time stamps
  all i, i': r.Int,   // row
      j, j': Int.r,   // column
      t, t': range[#r] {  // t is time stamp
    i->j->t in iter and i'->j'->t' in iter and t < t' =>
      i <= i' and (i = i' => j < j')
  }
}

-- third column of a ternary relation
fun c3 [r: Int->Int->Int]: Int {
  Int.(Int.r)
}

-- the last iteration of the inner loop for outer loop i'
fun end [iter: Int->Int->Int, i': Int]: Int {
  let times = {t: Int | some i, k: Int | i->k->t in iter and i <= i' } |
    some times => max[times] else 0 // the last iteartion for one row end
}



/*
 * check for ellcsr
 * Below we perform three checks:
 * 
 *   1. alphaValid checks that the abstraction functions map
 *      valid sparse representations to valid dense representations
 *
 *   2. preservesRep checks that the translation algorithm results
 *      in a valid CSR matrix
 *
 *   3. translateValid checks that the resulting CSR matrix
 *      represents the same dense matrix as the original ELL
 */  

assert alphaValid {
  all e: ELL, m: Matrix |
    I[e] and alpha[e, m] => I[e]
  all c: CSR, m: Matrix |
    I[c] and alpha[c, m] => I[m]
}

assert preservesRep {
  all e: ELL, c: CSR |
    I[e] and ellcsr[e, c] => I[c]
}

assert translateValid {
  all e: ELL, c: CSR, m: Matrix |
    I[e] and ellcsr[e, c] => (alpha[e, m] <=> alpha[c, m])
}

check alphaValid for 1 Matrix, 1 ELL, 1 CSR, 3 Value
check preservesRep for 0 Matrix, 1 ELL, 1 CSR, 3 Value
check translateValid for 1 Matrix, 1 ELL, 1 CSR, 3 Value

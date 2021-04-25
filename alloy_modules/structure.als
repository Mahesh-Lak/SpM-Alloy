module structure

/*
 * Specification of dense matrices and ELL and CSR sparse
 * format representations
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
 *    => Section III. Matrix Structure
 */

sig Value {}
one sig Zero extends Value {}

sig Matrix {
  rows, cols: Int,
  vals: Int->Int->lone Value
}

sig ELL {
  rows, cols, maxnz: Int,
  coef: Int->Int->lone Value,
  jcoef: Int->Int->lone Int
}

sig CSR {
  rows, cols: Int,
  A: Int->lone Value,
  IA, JA: Int->lone Int
}

-- arrays of the ELL format are 2D, arrays of CSR format are 1D
fact {
  all e: ELL | Array2D[e.coef] and Array2D[e.jcoef]
  all c: CSR | Array1D[c.A] and Array1D[c.IA] and Array1D[c.JA]
}

-- dense matrix abstract invariant
pred I [m: Matrix] {
  m.rows >= 0
  m.cols >= 0
  //m.rows != 0 <=> m.cols != 0
  m.rows = 0 <=> m.cols = 0 
  m.vals.univ = range[m.rows]->range[m.cols]
}


-- ELL concrete invariant
pred I [e: ELL] {
  e.rows >= 0
  e.cols >= 0
  //e.rows != 0 <=> e.cols != 0
  e.rows = 0 <=> e.cols = 0
  e.maxnz >= 0
  e.maxnz <= e.cols
  -- if all values are 0, this matrix is not necessary  
  e.coef.univ = range[e.rows]->range[e.maxnz]   -- set i->j for coef (Int -> Int) ?? When use .univ, when use coef[Int]
  e.jcoef.univ = range[e.rows]->range[e.maxnz]  -- set i->j for jcoef (Int -> Int)

  -- column indices are valid (include -1 placeholder)
  all j: e.jcoef[Int][Int] |
    j in range[-1, e.cols]

  -- column indices can appear once per row (-1 excepted)
  all i: range[e.rows] |   
    all j: range[e.cols] |   
      #e.jcoef[i].j <= 1    // coef[i] ith row. coef[i].j ith row join j value. this value will <=1

  -- enforce placeholders are at same locations within coef and jcoef
  all i: range[e.rows] |
    all k: range[e.maxnz] |
      e.coef[i][k] = Zero <=> e.jcoef[i][k] = -1

  -- there exists some rows in coef where all values are non-zero
  e.rows>0 => {some i: range[e.rows] |
    Zero not in (i.(e.coef))[Int]}

}

-- CSR concrete invariant
pred I [c: CSR] {
  c.rows >= 0
  c.cols >= 0
  //c.rows != 0 <=> c.cols != 0
  c.rows = 0 <=> c.cols = 0
  c.IA[0] = 0
  end[c.IA] = #c.A        -- last value of IA is length of A
  #c.IA = c.rows.add[1]  -- length of IA is rows + 1
  #c.A = #c.JA           -- A and JA are same length
  Zero not in c.A[Int]   -- A cannot contain 0

  -- values of IA are monotonically increasing
  all i: range[c.rows] |
    c.IA[i] <= c.IA[i.add[1]]

  -- column indices are valid
  all j: c.JA[Int] |
    j in range[c.cols]

  -- column indices appear once per row
  all i: range[c.rows] |
    let a = c.IA[i], b = c.IA[i.add[1]] |
      c.JA.slice[a, b].valsUnique


}

-- relation is a 1D array (indices range from 0 to n-1)
pred Array1D [a: Int->univ] { // used in CSR, received Int -> Value or Int -> Int 
                                          // could be replaced by seq Int??
  a.univ = range[#a]
}

-- relation is a 2D array (indices range from 0 to rows-1, 0 cols-1)
pred Array2D [a: Int->Int->univ] {
  a.univ.univ = range[#a.univ.univ]
  let s = #a[0] | all i: a.univ.univ | a[i].univ = range[s] //#a[0] is the number of columns
}

-- values in 1D array are unique
pred valsUnique [s: Int->univ] {
  #s = #s[Int]
}

-- last value in 1D array
fun end [s: Int->univ]: univ {
  s[max[s.univ]]
}

-- subset of 1D array in range [m, n-1]
fun slice [s: Int->univ, m, n: Int]: Int->univ { //s is an existing relation, return a set of (Int->Univ)
  0 <= m and m < n and n <= #s => { i: Int, v: univ |
      (i.add[m])->v in s and i in range[n.sub[m]] } else {none->none}
}

-- the set [0, n-1]
fun range [n: Int]: set Int {
  n>0 => { i: Int | 0 <= i and i < n } else none
}

-- the set [m, n-1]
fun range [m, n: Int]: set Int {
  n>0 => { i: Int | m <= i and i < n } else none
}

-- generate a 5x5 dense matrix
pred showDense {
  some m: Matrix {
    m.rows = 5
    m.cols = 5
    I[m]
  }
}

run showDense for 1 Matrix, 0 ELL, 0 CSR, 5 Value

-- generate a 5x5 ELL matrix
pred showELL {
  some e: ELL {
    e.rows = 5
    e.cols = 5
    I[e]
  }
}

run showELL for 0 Matrix, 1 ELL, 0 CSR, 5 Value

-- generate a 5x5 CSR matrix

pred showCSR {
  some c: CSR {
    c.rows = 5
    c.cols = 5
    I[c]
  }
}

run showCSR for 0 Matrix, 0 ELL, 1 CSR, 5 Value

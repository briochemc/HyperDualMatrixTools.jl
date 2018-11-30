# HyperDualMatrixTools.jl

<p>
  <a href="https://doi.org/10.5281/zenodo.1734670">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1734670.svg" alt="DOI">
  </a>
  <a href="https://briochemc.github.io/HyperDualMatrixTools.jl/stable">
    <img src=https://img.shields.io/badge/docs-stable-blue.svg>
  </a>
  <a href="https://ci.appveyor.com/project/briochemc/hyperdualmatrixtools-jl">
    <img src=https://ci.appveyor.com/api/projects/status/udbwakr621jbyvj1?svg=true>
  </a>
  <a href="https://travis-ci.com/briochemc/HyperDualMatrixTools.jl">
    <img alt="Build Status" src="https://travis-ci.com/briochemc/HyperDualMatrixTools.jl.svg?branch=master">
  </a>
  <a href='https://coveralls.io/github/briochemc/HyperDualMatrixTools.jl?branch=master'>
    <img src='https://coveralls.io/repos/github/briochemc/HyperDualMatrixTools.jl/badge.svg?branch=master' alt='Coverage Status' />
  </a>
  <a href="https://codecov.io/gh/briochemc/HyperDualMatrixTools.jl">
    <img src="https://codecov.io/gh/briochemc/HyperDualMatrixTools.jl/branch/master/graph/badge.svg" />
  </a>
  <a href="https://github.com/briochemc/HyperDualMatrixTools.jl/blob/master/LICENSE">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-yellow.svg">
  </a>
</p>

This package provides an overloaded `factorize` and `\` that work with hyperdual-valued arrays.

It uses the hyper dual type defined by the [HyperDualNumbers.jl](https://github.com/JuliaDiff/HyperDualNumbers.jl) package.
The idea is that for a hyperdual-valued matrix

<a href="https://www.codecogs.com/eqnedit.php?latex=\fn_phv&space;M&space;=&space;A&space;&plus;&space;\varepsilon_1&space;B&space;&plus;&space;\varepsilon_2&space;C&space;&plus;&space;\varepsilon_1&space;\varepsilon_2&space;D" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\fn_phv&space;M&space;=&space;A&space;&plus;&space;\varepsilon_1&space;B&space;&plus;&space;\varepsilon_2&space;C&space;&plus;&space;\varepsilon_1&space;\varepsilon_2&space;D" title="M = A + \varepsilon_1 B + \varepsilon_2 C + \varepsilon_1 \varepsilon_2 D" /></a>,

its inverse is given by

<a href="https://www.codecogs.com/eqnedit.php?latex=\fn_phv&space;M^{-1}&space;=&space;(I&space;-&space;\varepsilon_1&space;A^{-1}&space;B&space;-&space;\varepsilon_2&space;A^{-1}&space;C&space;-&space;\varepsilon_1&space;\varepsilon_2&space;A^{-1}&space;(D&space;-&space;C&space;A^{-1}&space;B&space;-&space;B&space;A^{-1}&space;C))&space;A^{-1}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\fn_phv&space;M^{-1}&space;=&space;(I&space;-&space;\varepsilon_1&space;A^{-1}&space;B&space;-&space;\varepsilon_2&space;A^{-1}&space;C&space;-&space;\varepsilon_1&space;\varepsilon_2&space;A^{-1}&space;(D&space;-&space;C&space;A^{-1}&space;B&space;-&space;B&space;A^{-1}&space;C))&space;A^{-1}" title="M^{-1} = (I - \varepsilon_1 A^{-1} B - \varepsilon_2 A^{-1} C - \varepsilon_1 \varepsilon_2 A^{-1} (D - C A^{-1} B - B A^{-1} C)) A^{-1}" /></a>.

Therefore, only the inverse of <a href="https://www.codecogs.com/eqnedit.php?latex=\fn_phv&space;A" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\fn_phv&space;A" title="A" /></a> is required to evaluate the inverse of <a href="https://www.codecogs.com/eqnedit.php?latex=\fn_phv&space;M" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\fn_phv&space;M" title="M" /></a>.
This package makes available a `HyperDualFactors` type which containts the factors of <a href="https://www.codecogs.com/eqnedit.php?latex=\fn_phv&space;A" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\fn_phv&space;A" title="A" /></a> and the non-real parts of <a href="https://www.codecogs.com/eqnedit.php?latex=\fn_phv&space;M" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\fn_phv&space;M" title="M" /></a>, and overloads `factorize` to create an instance of `HyperDualFactors`, which can then be called with `\` to efficiently solve hyperdual-valued linear systems of the type <a href="https://www.codecogs.com/eqnedit.php?latex=\fn_phv&space;M&space;x&space;=&space;b" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\fn_phv&space;M&space;x&space;=&space;b" title="M x = b" /></a>. 

This package should be useful for evaluation of second derivatives of functions that use `\` (e.g., with iterative solvers).

## Usage

- Create your hyperdual-valued matrix `M`:
    ```julia
    julia> M = A + ε₁ * B + ε₂ * C + ε₁ε₂ * D
    ```

- Factorize `M`:
    ```julia
    julia> Mf = factorize(M)
    ```

- Apply `\` to solve systems of the type `M * x = b`
    ```julia
    julia> x = Mf \ b
    ```





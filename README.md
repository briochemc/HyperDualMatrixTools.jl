# HyperDualMatrixTools.jl

This module overloads `factorize` and `\` to work with hyperdual-valued arrays.

It uses the HyperDualNumbers.jl package.

The idea is that for a hyperdual-valued matrix

<a href="https://www.codecogs.com/eqnedit.php?latex=\fn_phv&space;M&space;=&space;A&space;&plus;&space;\varepsilon_1&space;B&space;&plus;&space;\varepsilon_2&space;C&space;&plus;&space;\varepsilon_1&space;\varepsilon_2&space;D" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\fn_phv&space;M&space;=&space;A&space;&plus;&space;\varepsilon_1&space;B&space;&plus;&space;\varepsilon_2&space;C&space;&plus;&space;\varepsilon_1&space;\varepsilon_2&space;D" title="M = A + \varepsilon_1 B + \varepsilon_2 C + \varepsilon_1 \varepsilon_2 D" /></a>,

its inverse is given by

<a href="https://www.codecogs.com/eqnedit.php?latex=\fn_phv&space;M^{-1}&space;=&space;(I&space;-&space;\varepsilon_1&space;A^{-1}&space;B&space;-&space;\varepsilon_2&space;A^{-1}&space;C&space;-&space;\varepsilon_1&space;\varepsilon_2&space;A^{-1}&space;(D&space;-&space;C&space;A^{-1}&space;B&space;-&space;B&space;A^{-1}&space;C))&space;A^{-1}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\fn_phv&space;M^{-1}&space;=&space;(I&space;-&space;\varepsilon_1&space;A^{-1}&space;B&space;-&space;\varepsilon_2&space;A^{-1}&space;C&space;-&space;\varepsilon_1&space;\varepsilon_2&space;A^{-1}&space;(D&space;-&space;C&space;A^{-1}&space;B&space;-&space;B&space;A^{-1}&space;C))&space;A^{-1}" title="M^{-1} = (I - \varepsilon_1 A^{-1} B - \varepsilon_2 A^{-1} C - \varepsilon_1 \varepsilon_2 A^{-1} (D - C A^{-1} B - B A^{-1} C)) A^{-1}" /></a>

## Usage

- Create your hyperdual-valued matrix `M`:
    ```julia
    julia> M = A + ε₁ * B + ε₂ * C + ε₁ε₂ * D
    ```

- Factorize `M`:
    ```julia
    julia> Mf = factorize(M)
    ```

- Apply `\` to solve systems of the type `M * x = y`
    ```julia
    julia> x = M \ y
    ```





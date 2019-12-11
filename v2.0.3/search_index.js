var documenterSearchIndex = {"docs":
[{"location":"#HyperDualMatrixTools.jl-Documentation-1","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"","category":"section"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"This package provides an overloaded factorize and \\ that work with hyperdual-valued arrays.","category":"page"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"It is essentially base on the hyper dual type defined by the HyperDualNumbers.jl package.","category":"page"},{"location":"#Motivation-1","page":"HyperDualMatrixTools.jl Documentation","title":"Motivation","text":"","category":"section"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"The idea is that for a hyperdual-valued matrix M = A + varepsilon_1 B + varepsilon_2 C + varepsilon_1 varepsilon_2 D, its inverse is given by M^-1 = (I - varepsilon_1 A^-1 B - varepsilon_2 A^-1 C - varepsilon_1varepsilon_2 A^-1 (D - B A^-1 C - C A^-1 B)) A^-1. Therefore, only the inverse of A is required to evaluate the inverse of M. This package should be useful for evaluation of second derivatives of functions that use \\ (e.g., with iterative solvers).","category":"page"},{"location":"#How-it-works-1","page":"HyperDualMatrixTools.jl Documentation","title":"How it works","text":"","category":"section"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"HyperDualMatrixTools.jl makes available a HyperDualFactors type which contains the factors of A (i.e., the output of factorize, e.g., L and U, or Q and R) and the non-real parts of M (i.e., B, C, and D). HyperDualMatrixTools.jl overloads factorize so that for a hyperdual-valued matrix M, factorize(M) creates an instance of HyperDualFactors. Finally, HyperDualMatrixTools.jl also overloads \\ to efficiently solve hyperdual-valued linear systems of the type M x = y by using the default \\ with the factors of A only.","category":"page"},{"location":"#Usage-1","page":"HyperDualMatrixTools.jl Documentation","title":"Usage","text":"","category":"section"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"DocTestSetup = quote\n    using HyperDualMatrixTools\nend","category":"page"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"Create your hyperdual-valued matrix M:","category":"page"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"n = 4\nA, B, C, D = rand(n, n), randn(n, n), rand(n, n), randn(n, n)\nM = A + ε₁ * B + ε₂ * C + ε₁ε₂ * D\ntypeof(M)\n\n# output\n\nArray{HyperDualNumbers.Hyper{Float64},2}","category":"page"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"(The ε₁, ε₂, and ε₁ε₂ constants are provided by HyperDualMatrixTools.jl for convenience.)","category":"page"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"Factorize M:","category":"page"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"Mf = factorize(M)\ntypeof(Mf)\n\n# output\n\nHyperDualFactors","category":"page"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"Apply \\ to solve systems of the type M * x = y","category":"page"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"y = rand(n, 4) * [1.0, ε₁, ε₂, ε₁ε₂]\nx = Mf \\ y\nM * x ≈ y\n\n# output\n\ntrue","category":"page"},{"location":"#Functions-1","page":"HyperDualMatrixTools.jl Documentation","title":"Functions","text":"","category":"section"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"factorize","category":"page"},{"location":"#LinearAlgebra.factorize","page":"HyperDualMatrixTools.jl Documentation","title":"LinearAlgebra.factorize","text":"factorize(M::SparseMatrixCSC{<:Hyper,<:Int})\n\nInvokes factorize on just the real part of M and stores it along with the dual parts into a HyperDualFactors object.\n\n\n\n\n\nfactorize(M::Array{<:Hyper,2})\n\nInvokes factorize on just the real part of M and stores it along with the dual parts into a HyperDualFactors object.\n\n\n\n\n\n","category":"function"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"\\","category":"page"},{"location":"#Base.:\\","page":"HyperDualMatrixTools.jl Documentation","title":"Base.:\\","text":"\\(M::HyperDualFactors, y::AbstractVecOrMat{Float64})\n\nBacksubstitution for HyperDualFactors. See HyperDualFactors for details.\n\n\n\n\n\n\\(M::HyperDualFactors, y::AbstractVecOrMat{Hyper256})\n\nBacksubstitution for HyperDualFactors. See HyperDualFactors for details.\n\n\n\n\n\n\\(Af::Factorization{Float64}, y::AbstractVecOrMat{Hyper256})\n\nBacksubstitution for HyperDual-valued RHS.\n\n\n\n\n\n\\(M::AbstractArray{<:Hyper,2}, y::AbstractVecOrMat)\n\nBackslash (factorization and backsubstitution) for Dual-valued matrix M.\n\n\n\n\n\n","category":"function"},{"location":"#New-types-1","page":"HyperDualMatrixTools.jl Documentation","title":"New types","text":"","category":"section"},{"location":"#","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.jl Documentation","text":"HyperDualFactors","category":"page"},{"location":"#HyperDualMatrixTools.HyperDualFactors","page":"HyperDualMatrixTools.jl Documentation","title":"HyperDualMatrixTools.HyperDualFactors","text":"HyperDualFactors\n\nContainer type to work efficiently with backslash on hyperdual-valued sparse matrices.\n\nfactorize(M) will create an instance containing\n\nAf = factorize(realpart.(M)) — the factors of the real part\nB = ε₁part.(M) — the varepsilon_1 part\nC = ε₂part.(M) — the varepsilon_2 part\nD = ε₁ε₂part.(M) — the varepsilon_1varepsilon_2 part\n\nfor a hyperdual-valued matrix M.\n\nThis is because only the factors of the real part are needed when solving a linear system of the type M x = b for a hyperdual-valued matrix M = A + varepsilon_1 B + varepsilon_2 C + varepsilon_1 varepsilon_2 D. In fact, the inverse of M is given by M^-1 = (I - varepsilon_1 A^-1 B - varepsilon_2 A^-1 C - varepsilon_1varepsilon_2 A^-1 (D - B A^-1 C - C A^-1 B)) A^-1.\n\n\n\n\n\n","category":"type"}]
}

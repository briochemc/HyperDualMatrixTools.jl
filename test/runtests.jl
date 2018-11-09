using Test

using HyperDualMatrixTools
using HyperDualNumbers, LinearAlgebra, SparseArrays, SuiteSparse

n = 4

A = rand(n, n)
B = rand(n, n)
C = rand(n, n)
D = rand(n, n)

hyperx = rand(n, n) * [1.0, ε₁, ε₂, ε₁ε₂]
realx = rand(n)

hyperM = A + ε₁ * B + ε₂ * C + ε₁ε₂ * D
realM = rand(n, n)

realMrealx = realM * realx
realMhyperx = realM * hyperx
hyperMrealx = hyperM * realx
hyperMhyperx = hyperM * hyperx

# Solving `M * x = y` via `x = M \ y`
realMfac = factorize(realM)
hyperMfac = factorize(hyperM)

realx2 = realMfac \ realMrealx
realx3 = hyperMfac \ hyperMrealx
hyperx2 = realMfac \ realMhyperx
hyperx3 = hyperMfac \ hyperMhyperx


spA = sparse(rand(n, n))
spB = sparse(rand(n, n))
spC = sparse(rand(n, n))
spD = sparse(rand(n, n))

hyperspM = spA + ε₁ * spB + ε₂ * spC + ε₁ε₂ * spD
realspM = sparse(rand(n, n))

realspMrealx = realspM * realx
realspMhyperx = realspM * hyperx
hyperspMrealx = hyperspM * realx
hyperspMhyperx = hyperspM * hyperx

# Solving `spM * x = y` via `x = spM \ y`
realspMfac = factorize(realspM)
hyperspMfac = factorize(hyperspM)

realspx2 = realspMfac \ realspMrealx
realspx3 = hyperspMfac \ hyperspMrealx
hyperspx2 = realspMfac \ realspMhyperx
hyperspx3 = hyperspMfac \ hyperspMhyperx

@test realx2 ≈ realx
@test realx3 ≈ realx
@test hyperx2 ≈ hyperx
@test hyperx3 ≈ hyperx

@test realspx2 ≈ realx
@test realspx3 ≈ realx
@test hyperspx2 ≈ hyperx
@test hyperspx3 ≈ hyperx

module HyperDualMatrixTools

using HyperDualNumbers, LinearAlgebra, SparseArrays, SuiteSparse

const ε₁ = hyper(0.0, 1.0, 0.0, 0.0)
const ε₂ = hyper(0.0, 0.0, 1.0, 0.0)
const ε₁ε₂ = ε₁ * ε₂

import LinearAlgebra.factorize
import Base.\

"""
    HyperDualFactors

Container type to work efficiently with backslash on hyperdual-valued sparse matrices.

`factorize(M)` will create an instance containing
- `Af = factorize(real.(M))` — the factors of the real part
- `B = eps1.(M)` — the ``\\varepsilon_1`` part
- `C = eps2.(M)` — the ``\\varepsilon_2`` part
- `D = eps1eps2.(M)` — the ``\\varepsilon_1\\varepsilon_2`` part
for a hyperdual-valued matrix `M`.

This is because only the factors of the real part are needed when solving a linear system of the type ``M x = b`` for a hyperdual-valued matrix ``M = A + B \\varepsilon_1 + C \\varepsilon_2 + D \\varepsilon_1 \\varepsilon_2``.
In fact, the inverse of ``M`` is given by
``M^{-1} = (I - \\varepsilon_1 A^{-1} B - \\varepsilon_2 A^{-1} C + \\varepsilon_1\\varepsilon_2 (-A^{-1} D + A^{-1} B A^{-1} C + A^{-1} C A^{-1} B)) A^{-1}``.
"""
mutable struct HyperDualFactors
    Af # the factors of the real part
    B  # the ε₁ part
    C  # the ε₂ part
    D  # the ε₁ε₂ part
end


"""
    factorize(M::Union{SparseMatrixCSC{Hyper{Float64},Int64}, Array{Hyper{Float64},2}})

Efficient factorization of hyperdual-valued sparse matrices.
See `HyperDualFactors` for details.
"""
function factorize(M::Union{SparseMatrixCSC{Hyper{Float64},Int64}, Array{Hyper{Float64},2}})
    return HyperDualFactors(factorize(real.(M)), eps1.(M), eps2.(M), eps1eps2.(M))
end

"""
    \\(M::HyperDualFactors, a::AbstractVecOrMat{Float64})

Backsubstitution for `HyperDualFactors`.
See `HyperDualFactors` for details.
"""
function \(M::HyperDualFactors, y::AbstractVecOrMat{Float64})
    A, B, C, D = M.Af, M.B, M.C, M.D
    A⁻¹y = A \ y
    DA⁻¹y = D * A⁻¹y
    A⁻¹BA⁻¹y = A \ (B * A⁻¹y)
    A⁻¹CA⁻¹y = A \ (C * A⁻¹y)
    return A⁻¹y - ε₁ * A⁻¹BA⁻¹y - ε₂ * A⁻¹CA⁻¹y - ε₁ε₂ * (A \ DA⁻¹y) +
        ε₁ε₂ * (A \ (B * A⁻¹CA⁻¹y) + A \ (C * A⁻¹BA⁻¹y))
end

"""
    \\(M::HyperDualFactors, y::AbstractVecOrMat{HyperDual{Float64}})

Backsubstitution for `HyperDualFactors`.
See `HyperDualFactors` for details.
"""
function \(M::HyperDualFactors, y::AbstractVecOrMat{Hyper{Float64}})
    a, b, c, d = real.(y), eps1.(y), eps2.(y), eps1eps2.(y)
    A, B, C, D = M.Af, M.B, M.C, M.D
    A⁻¹a = A \ a
    DA⁻¹a = D * A⁻¹a
    A⁻¹BA⁻¹a = A \ (B * A⁻¹a)
    A⁻¹CA⁻¹a = A \ (C * A⁻¹a)
    A⁻¹b = A \ b
    A⁻¹CA⁻¹b = A \ (C * A⁻¹b)
    A⁻¹c = A \ c
    A⁻¹BA⁻¹c = A \ (B * A⁻¹c)
    A⁻¹d = A \ d
    return A⁻¹a - ε₁ * A⁻¹BA⁻¹a - ε₂ * A⁻¹CA⁻¹a - ε₁ε₂ * (A \ DA⁻¹a) +
        ε₁ε₂ * (A \ (B * A⁻¹CA⁻¹a) + A \ (C * A⁻¹BA⁻¹a)) +
        ε₁ * A⁻¹b + ε₂ * A⁻¹c + ε₁ε₂ * A⁻¹d - ε₁ε₂ * (A⁻¹CA⁻¹b + A⁻¹BA⁻¹c)
end

"""
    \\(Af::SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64}, y::AbstractVecOrMat{HyperDual{Float64}})

Backsubstitution for HyperDual-valued RHS.
"""
function \(Af::Union{SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64}, LU{Float64,Array{Float64,2}}}, y::AbstractVecOrMat{Hyper{Float64}})
    return (Af \ real.(y)) + ε₁ * (Af \ eps1.(y)) + ε₂ * (Af \ eps2.(y)) + ε₁ε₂ * (Af \ eps1eps2.(y))
end

function Base.isapprox(x::AbstractVecOrMat{Hyper{Float64}}, y::AbstractVecOrMat{Hyper{Float64}})
    bigx = [real.(x) eps1.(x) eps2.(x) eps1eps2.(x)]
    bigy = [real.(y) eps1.(y) eps2.(y) eps1eps2.(y)]
    return isapprox(bigx, bigy)
end
Base.isapprox(x::AbstractVecOrMat, y::AbstractVecOrMat{Hyper{Float64}}) = isapprox(hyper.(x), y)
Base.isapprox(x::AbstractVecOrMat{Hyper{Float64}}, y::AbstractVecOrMat) = isapprox(x, hyper.(y))
function Base.isapprox(x::Hyper{Float64}, y::Hyper{Float64})
    bigx = [real(x) eps1(x) eps2(x) eps1eps2(x)]
    bigy = [real(y) eps1(y) eps2(y) eps1eps2(y)]
    return isapprox(bigx, bigy)
end
Base.isapprox(x::Float64, y::Hyper{Float64}) = isapprox(hyper(x), y)
Base.isapprox(x::Hyper{Float64}, y::Float64) = isapprox(x, hyper(y))

export ε₁, ε₂, ε₁ε₂, HyperDualFactors

end # module

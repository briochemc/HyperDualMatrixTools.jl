using Test

using HyperDualMatrixTools
using HyperDualNumbers, LinearAlgebra, SparseArrays, SuiteSparse

@testset "Testing HyperDualMatrixTools" begin
    # Chose a size for matrices
    n = 10

    # Create a real-valued random vector
    y = randn(n)
    # Create a hyperdual-valued random vector
    x = randn(n, 4) * [1.0, ε₁, ε₂, ε₁ε₂]

    @testset "Full matrices" begin
        # Create a real-valued random matrix
        A = randn(n, n)
        # Create a hyperdual-valued random matrix
        B = randn(n, n)
        C = randn(n, n)
        D = randn(n, n)
        M = A + ε₁ * B + ε₂ * C + ε₁ε₂ * D

        # Check that `\` works for full matrices
        Af = factorize(A)
        Mf = factorize(M)
        @test Af \ (A * x) ≈ x
        @test Mf \ (M * x) ≈ x
        @test Af \ (A * y) ≈ y
        @test Mf \ (M * y) ≈ y
        @test A * (Af \ x) ≈ x
        @test M * (Mf \ x) ≈ x
        @test A * (Af \ y) ≈ y
        @test M * (Mf \ y) ≈ y
    end

    @testset "Sparse matrices" begin
        # Create a real-valued sparse matrix
        A = sparse(randn(n, n))
        # Create a hyperdual-valued sparse matrix
        B = sparse(randn(n, n))
        C = sparse(randn(n, n))
        D = sparse(randn(n, n))
        M = A + ε₁ * B + ε₂ * C + ε₁ε₂ * D

        # Check that `\` works for sparse matrices
        Af = factorize(A)
        Mf = factorize(M)
        @test Af \ (A * x) ≈ x
        @test Mf \ (M * x) ≈ x
        @test Af \ (A * y) ≈ y
        @test Mf \ (M * y) ≈ y
        @test A * (Af \ x) ≈ x
        @test M * (Mf \ x) ≈ x
        @test A * (Af \ y) ≈ y
        @test M * (Mf \ y) ≈ y
    end
end

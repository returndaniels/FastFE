@doc raw"""
Módulo principal do `FastFE` -- um pacote de alto desempenho para cálculos utilizando o método dos elementos finitos (FEM).

Este pacote fornece ferramentas para resolver problemas usando o método dos elementos finitos, com foco em eficiência e análise de erros. Ele inclui implementações de rotinas FEM tanto vetorizadas quanto seriais.

# Visão Geral

- **Método dos Elementos Finitos (FEM):**
  - Cálculo eficiente de matrizes de rigidez e vetores de forças externas.
  - Suporte para abordagens seriais e vetorizadas de FEM.

- **Análise de Erros:**
  - Ferramentas para analisar erros de aproximação e taxas de convergência.
  - Suporte para métodos de análise de erros vetorizados e seriais.

- **Discretização de Equações Diferenciais:**
  - Utilitários para discretização espacial e temporal de equações diferenciais modernas (MDE).

# Exportações

$(EXPORTS)
"""
module FastFE

# Import standard libraries
import LinearAlgebra
import SparseArrays
import Plots

# Import external dependencies
using DocStringExtensions: EXPORTS
using Test: @testset, @test

# Import internal modules
include("./VectorizedFiniteElements.jl")
include("./SerialFiniteElements.jl")
include("./MDEDiscretization.jl")
include("./VectorizedErrorAnalysis.jl")
include("./SerialErrorAnalysis.jl")

# Expose relevant symbols from internal modules
export K_serial, F_serial!, G_serial!, erro_serial, C0_options,
    K_vectorized, F_vectorized!, G_vectorized!,
    calculate_errors_serial, erro_vectorized, calculate_errors_vectorized,
    beta, alpha, gamma, a, b, T, npg, ne, m, h, tau, N, P, W, 
    φ1, φ2, dφ1, dφ2, φ1P, φ2P, dφ1P, dφ2P, Wφ1P, Wφ2P, WφP, 
    u, u0, du0, g, f, monta_LG, monta_EQ

# Submodules and Aliases
using .VectorizedFiniteElements: K_vectorized, F_vectorized!, G_vectorized!
using .SerialFiniteElements: K_serial, F_serial!, G_serial!, erro_serial, C0_options
using .SerialErrorAnalysis: calculate_errors_serial
using .VectorizedErrorAnalysis: erro_vectorized, calculate_errors_vectorized
using .MDEDiscretization: 
    beta, alpha, gamma, a, b, T, npg, ne, m, h, tau, N, P, W, φ1, φ2, dφ1, dφ2, 
    φ1P, φ2P, dφ1P, dφ2P, Wφ1P, Wφ2P, WφP, u, u0, du0, g, f, monta_LG, monta_EQ

# Version number of FastFE
# const FASTFE_VERSION = let
#     project = joinpath(dirname(dirname(pathof(FastFE))), "Project.toml")
#     Base.include_dependency(project)
#     toml = read(project, String)
#     m = match(r"(*ANYCRLF)^version\s*=\s\"(.*)\"$"m, toml)
#     VersionNumber(m[1])
# end

# Documentation of key exports
@doc raw"""
    K_vectorized(ne, m, h, npg, alpha, beta, gamma, EQoLG)

Calcula a matriz global de rigidez utilizando o método vetorizado.
"""
K_vectorized

@doc raw"""
    F_vectorized!(F_ext_vectorized, X, f_eval, values, f, ne, m, h, npg, EQoLG)

Calcula o vetor de forças externas utilizando quadratura de Gauss vetorizada.
"""
F_vectorized!

@doc raw"""
    G_vectorized!(G_ext_vectorized, g_eval, values, C, ne, m, h, npg, EQoLG)

Calcula o vetor estendido de forças associado à função `g` usando vetorização.
"""
G_vectorized!

@doc raw"""
    K_serial()

Calcula a matriz de rigidez global K a partir das matrizes locais, utilizando pontos de Gauss e pesos para a integração numérica.
Essa função é implementada de forma sequencial.
"""
K_serial

@doc raw"""
    F_serial!()

Calcula a quadratura gaussiana do vetor de forças externas `F_ext` associada à equação diferencial.
"""
F_serial!

@doc raw"""
    G_serial!

Calcula a quadratura gaussiana associada ao vetor `G_ext` considerando os coeficientes `C` da solução aproximada.
"""
G_serial!

@doc raw"""
    erro_serial()

Calcula o erro integral L2 entre a solução exata `u` e a solução aproximada representada pelos coeficientes `C`.
"""
erro_serial

@doc raw"""
    C0_options()

Calcula o vetor `C0` de coeficientes iniciais, dado o tipo de inicialização especificado.
"""
C0_options

@doc raw"""
    calculate_errors_serial()

Calculates and plots the error convergence for different discretizations using a serial approach.
"""
calculate_errors_serial

@doc raw"""
    erro_vectorized()

Calculates the error for a given approximation using a vectorized approach.
"""
erro_vectorized

@doc raw"""
    calculate_errors()

Calculates and plots the error convergence for different discretizations.
"""
calculate_errors

# Precompilation setup
# import PrecompileTools
# PrecompileTools.@compile_workload begin
#     include("precompile/make.jl")
# end

end # module

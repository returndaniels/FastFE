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
    calculate_errors_serial, erro_vectorized, calculate_errors,
    beta, alpha, gamma, a, b, T, npg, ne, m, h, tau, N, P, W, 
    φ1, φ2, dφ1, dφ2, φ1P, φ2P, dφ1P, dφ2P, Wφ1P, Wφ2P, WφP, 
    u, u0, du0, g, f, monta_LG, monta_EQ

# Submodules and Aliases
using .VectorizedFiniteElements: K_vectorized, F_vectorized!, G_vectorized!
using .SerialFiniteElements: K_serial, F_serial!, G_serial!, erro_serial, C0_options
using .SerialErrorAnalysis: calculate_errors_serial
using .VectorizedErrorAnalysis: erro_vectorized, calculate_errors
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


# Precompilation setup
# import PrecompileTools
# PrecompileTools.@compile_workload begin
#     include("precompile/make.jl")
# end

end # module

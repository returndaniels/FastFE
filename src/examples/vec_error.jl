include("../FastFE.jl")

using ..FastFE
using Plots

tam = 5
option = 1

erros, hs = calculate_errors_vectorized(tam, u, u0, f, FastFE.alpha, FastFE.beta, FastFE.gamma, FastFE.a, FastFE.b, npg, option)

plot(hs, hs.^2, label="h^2", xlabel="h", title="Ordem de Convergência do Erro (Escala Log)", legend=:topright)
plot!(hs, erros, yscale=:log10, xscale=:log10, label="Erro")
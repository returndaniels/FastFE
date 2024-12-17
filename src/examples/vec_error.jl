include("../FastFE.jl")

using ..FastFE

tam = 6
option = 1
EQoLG = monta_EQ(ne)[monta_LG(ne)]
K = K_vectorized(ne, m, h, npg, FastFE.alpha, beta, gamma, EQoLG)
C0 = C0_options(1, u0, a, ne, m, h, npg, EQoLG)

erros, hs = calculate_errors(tam, u, u0, f, EQoLG, K, C0, alpha, beta, gamma, a, b, npg, option)

plot(hs, hs.^2, label="h^2", xlabel="h", title="Ordem de Convergência do Erro (Escala Log)", legend=:topright)
plot!(hs, erros, yscale=:log10, xscale=:log10, label="Erro")
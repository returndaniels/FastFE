include("./FastFE.jl")

using .FastFE

using LinearAlgebra
using FastGaussQuadrature
using SparseArrays
using DataFrames

using BenchmarkTools # Usar @btime ou @benchmark na função que deseja avaliar

option = 1 # Opção de como calcula C0
# Define as discretizações
EQoLG = monta_EQ(ne)[monta_LG(ne)]
# Aproximando U0
C0 = C0_options(option, u0, a, ne, m, h, npg, EQoLG)
############################ Variáveis das funções ############################
x = h*(P .+ 1)/2 .+ a
F_ext_serial = zeros(Float64, ne)
G_ext_serial = zeros(Float64, ne);

output_K = @benchmark K_serial(ne, m, h, npg, alpha, beta, gamma, EQoLG)
display(output_K)

output_F = @benchmark F_serial!(F_ext_serial, x, x -> f(x, 0.0), ne, m, h, npg, EQoLG)
display(output_F)

output_G = @benchmark G_serial!(G_ext_serial, C0, ne, m, h, npg, EQoLG)
display(output_G)

output_sys = @benchmark begin
    # Monta estrutura local global
    EQoLG = monta_EQ(ne)[monta_LG(ne)]

    # Estruturas externas usadas pelas funções
    x = h*(P .+ 1)/2 .+ a
    F_ext_serial = zeros(Float64, ne)
    G_ext_serial = zeros(Float64, ne)
    
    # Montando o sistema linear
    M = K_serial(ne, m, h, npg, 0., 1., 0., EQoLG)
    K = K_serial(ne, m, h, npg, alpha, beta, gamma, EQoLG)
    MK = M/tau - K/2
    # Decomposição LU para resolução de múltiplos sistemas
    LU_dec = lu(M/tau + K/2)
    
    # Aproximação de U0
    C0 = C0_options(1, u0, a, ne, m, h, npg, EQoLG)
    # Aproximação de U1
    C1_tiu = LU_dec\(F_serial!(F_ext_serial, x, x -> f(x, tau*0.5), ne, m, h, npg, EQoLG) + MK*C0 -
             G_serial!(G_ext_serial, C0, ne, m, h, npg, EQoLG))
    C1 = LU_dec\(F_serial!(F_ext_serial, x, x -> f(x, tau*0.5), ne, m, h, npg, EQoLG) + MK*C0 -
         G_serial!(G_ext_serial, (C0 + C1_tiu)/2, ne, m, h, npg, EQoLG))

    for n in 0:N
        Cn = LU_dec\(F_serial!(F_ext_serial, x, x -> f(x, tau*(n-0.5)), ne, m, h, npg, EQoLG) + MK*C1 -
             G_serial!(G_ext_serial, (3*C1 - C0)/2, ne, m, h, npg, EQoLG))
        C0 = C1
        C1 = Cn
    end
end
display(output_sys)

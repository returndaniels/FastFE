include("../FastFE.jl")

using ..FastFE
using LinearAlgebra
using Plots

# Intervalo espacial da solução exata
malha_x = a:0.01:b
# Intervalos da solução aproximada
intervalo_x = a:h:b
intervalo_t = 0:tau:T

# Monta estrutura local global
EQoLG = monta_EQ(ne)[monta_LG(ne)]

# Estruturas externas usadas pelas funções
x = h*(P .+ 1)/2 .+ a
F_ext_serial = zeros(Float64, ne)
G_ext_serial = zeros(Float64, ne)

# Montando o sistema linear
M = K_serial(ne, m, h, npg, 0., 1., 0., EQoLG)
K = K_serial(ne, m, h, npg, FastFE.alpha, beta, gamma, EQoLG)
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

let C0 = C0, C1 = C1 # Resolvendo problemas com troca de contexto
    anim = @animate for n in 0:N
        plt = plot(size=(800, 400), title="Soluções Exata e Aproximada em t_($n) = $(intervalo_t[n+1])",
                       xlabel="x", legend=:topright, ylim=(0, 0.110))
        # Solução exata
        plot!(plt, malha_x, u.(malha_x, intervalo_t[n+1]), lw=3, color=:blue, label="u(x,t_$(n))")
        if(n == 0)
            # Solução aproximada
            plot!(plt, intervalo_x, [0; C0; 0], lw=3, color=:red, linestyle=:dash, markershape=:circle, label="U^$(n)(x)")
        elseif(n == 1)
            # Solução aproximada
            plot!(plt, intervalo_x, [0; C1; 0], lw=3, color=:red, linestyle=:dash, markershape=:circle, label="U^$(n)(x)")
        else
           Cn = LU_dec\(F_serial!(F_ext_serial, x, x -> f(x, tau*(n-0.5)), ne, m, h, npg, EQoLG) + MK*C1 -
                G_serial!(G_ext_serial, (3*C1 - C0)/2, ne, m, h, npg, EQoLG))
            # Solução aproximada
            plot!(plt, intervalo_x, [0; Cn; 0], lw=3, color=:red, linestyle=:dash, markershape=:circle, label="U^$(n)(x)")
            C0 = C1
            C1 = Cn
        end
    end
    gif(anim, "solucao_exata_vs_aproximada.gif", fps=1)
end
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
X = ((h/2)*(P .+ 1) .+ a)' .+ range(a, step=h, stop=b-h)
f_eval = Matrix{Float64}(undef, ne, npg)
g_eval = Matrix{Float64}(undef, ne, npg)
values = Matrix{Float64}(undef, ne, 2)
F_ext_vectorized = zeros(Float64, ne)
G_ext_vectorized = zeros(Float64, ne)

# Pré-aloca mmeória para o lado direito do sistema linear
B = zeros(Float64, ne-1)
# Montando o sistema linear
M = K_vectorized(ne, m, h, npg, 0., 1., 0., EQoLG)[1:ne-1, 1:ne-1]
K = K_vectorized(ne, m, h, npg, FastFE.alpha, FastFE.beta, FastFE.gamma, EQoLG)[1:ne-1, 1:ne-1]
MK = M/tau - K/2
# Decomposição LU para resolução de múltiplos sistemas
LU_dec = lu(M/tau + K/2)
# Aproximação de U0
C0_ext = [C0_options(1, u0, a, ne, m, h, npg, EQoLG); 0]
C0 = @view C0_ext[1:ne-1]
MKC0 = MK*C0
# Aproximação de U1
F_vectorized!(F_ext_vectorized, X, f_eval, values, x -> f(x, tau*0.5), ne, m, h, npg, EQoLG)
F = @view F_ext_vectorized[1:ne-1]
G_vectorized!(G_ext_vectorized, g_eval, values, C0_ext, ne, m, h, npg, EQoLG)
G = @view G_ext_vectorized[1:ne-1]
C1_ext_tiu = zeros(Float64, ne)
B .= F
B .+= MKC0
B .-= G
C1_ext_tiu[1:ne-1] .= LU_dec\B
# Calculando C1
G_vectorized!(G_ext_vectorized, g_eval, values, (C0_ext + C1_ext_tiu)/2, ne, m, h, npg, EQoLG)
G = @view G_ext_vectorized[1:ne-1]
C1_ext = zeros(Float64, ne)
B .= F
B .+= MKC0
B .-= G
C1_ext[1:ne-1] = LU_dec\B
# Pré-aloca memória para operações repetidas
Cn_ext = zeros(Float64, ne)
G_C_ext = zeros(Float64, ne)
MKC1 = zeros(Float64, ne-1)
B = zeros(Float64, ne-1)
let C0 = C0, C1_ext = C1_ext # Resolvendo problemas com troca de contexto
    anim = @views @animate for n in 0:N
        plt = plot(size=(800, 400), title="Soluções Exata e Aproximada em t_($n) = $(intervalo_t[n+1])",
                    xlabel="x", legend=:topright, ylim=(0, 0.110))
        # Solução exata
        plot!(plt, malha_x, u.(malha_x, intervalo_t[n+1]), lw=3, color=:blue, label="u(x,t_$(n))")
        if(n == 0)
            # Solução aproximada
            plot!(plt, intervalo_x, [0; C0; 0], lw=3, color=:red, linestyle=:dash, markershape=:circle, label="U^$(n)(x)")
        elseif(n == 1)
            # Solução aproximada
            C1_ext[end] = 0
            plot!(plt, intervalo_x, [0; C1_ext], lw=3, color=:red, linestyle=:dash, markershape=:circle, label="U^$(n)(x)")
        else
            ############ Cálculo sem consumo extra de memória ############
            F_vectorized!(F_ext_vectorized, X, f_eval, values, x -> f(x, tau*(n-0.5)), ne, m, h, npg, EQoLG)
            F = F_ext_vectorized[1:ne-1]
            MKC1 .= MK*(C1_ext[1:ne-1])
            G_C_ext[1:ne-1] .= 3*(C1_ext[1:ne-1])
            G_C_ext[1:ne-1] .-= C0_ext[1:ne-1]
            G_C_ext[1:ne-1] ./= 2
            G_vectorized!(G_ext_vectorized, g_eval, values, G_C_ext, ne, m, h, npg, EQoLG)
            G = G_ext_vectorized[1:ne-1]
            B .= F
            B .+= MKC1
            B .-= G
            ##############################################################
            # Calcula o próximo coeficiente e atualiza
            Cn_ext[1:ne-1] .= LU_dec\B
            C0_ext[1:ne-1] .= C1_ext[1:ne-1]
            C1_ext[1:ne-1] .= Cn_ext[1:ne-1]
            C1_ext[end] = 0
            plot!(plt, intervalo_x, [0; C1_ext], lw=3, color=:red, linestyle=:dash, markershape=:circle, label="U^$(n)(x)")
        end
    end
    gif(anim, "solucao_exata_vs_aproximada_(vetorizada).gif", fps=1)
end
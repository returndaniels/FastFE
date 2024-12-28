include("./FastFE.jl")

using .FastFE

using LinearAlgebra
using FastGaussQuadrature
using SparseArrays
using DataFrames

using BenchmarkTools # Usar @btime ou @benchmark na função que deseja avaliar

samples = 1000 # Quantidade de amostras para o benchmark de funções
evals = 1 # Quantidade de avaliações por amostras para o benchmark de funções
samples_sys = 10 # Quantidade de amostras para o benchmark do sistema linear
evals_sys = 1 # Quantidade de avaliações por amostras para o benchmark do sistema linear
# Recebe uma função para ser testada no benchmark
function benchmark(func::Function, args::Any, samples::Int64; evals=1)
    func(args...) # Compila a função antes de começar para ignorar tempo de compilação
    resultado = @benchmark ($func($args...)) samples=samples evals=evals
    return resultado
end

option = 1 # Opção de como calcula C0
# Monta estrutura local global
EQoLG = monta_EQ(ne)[monta_LG(ne)]
# Calcula o C0 inicial
C0 = C0_options(option, u0, a, ne, m, h, npg, EQoLG)

# Estruturas externas usadas pelas funções
X = ((h/2)*(P .+ 1) .+ a)' .+ range(a, step=h, stop=b-h)
f_eval = Matrix{Float64}(undef, ne, npg)
g_eval = Matrix{Float64}(undef, ne, npg)
u_eval = Matrix{Float64}(undef, ne, npg)
values = Matrix{Float64}(undef, ne, 2)
F_ext_vectorized = zeros(Float64, ne)
G_ext_vectorized = zeros(Float64, ne);

output_K = benchmark(K_vectorized, (ne, m, h, npg, alpha, beta, gamma, EQoLG), samples, evals=evals)
display(output_K)

output_F = benchmark(F_vectorized!, (F_ext_vectorized, X, f_eval, values, x -> f(x, 0.0), ne, m, h, npg, EQoLG), samples, evals=evals)
display(output_F)

output_G = benchmark(G_vectorized!, (G_ext_vectorized, g_eval, values, C0, ne, m, h, npg, EQoLG), samples, evals=evals)
display(output_G)

output_sys = @benchmark begin
    # Monta estrutura local global
    EQoLG = monta_EQ(ne)[monta_LG(ne)]
    
    # Estruturas externas usadas pelas funções
    X = ((h/2)*(P .+ 1) .+ a)' .+ range(a, step=h, stop=b-h)
    f_eval = Matrix{Float64}(undef, ne, npg)
    g_eval = Matrix{Float64}(undef, ne, npg)
    u_eval = Matrix{Float64}(undef, ne, npg)
    values = Matrix{Float64}(undef, ne, 2)
    F_ext_vectorized = zeros(Float64, ne)
    G_ext_vectorized = zeros(Float64, ne)
    
    # Montando o sistema linear
    M = K_vectorized(ne, m, h, npg, 0., 1., 0., EQoLG)
    K = K_vectorized(ne, m, h, npg, alpha, beta, gamma, EQoLG)
    MK = M/tau - K/2
    # Decomposição LU para resolução de múltiplos sistemas
    LU_dec = lu(M/tau + K/2)
    
    # Aproximação de U0
    C0 = C0_options(1, u0, a, ne, m, h, npg, EQoLG)
    # Aproximação de U1
    C1_tiu = LU_dec\(F_vectorized!(F_ext_vectorized, X, f_eval, values, x -> f(x, tau*0.5), ne, m, h, npg, EQoLG) + MK*C0 -
             G_vectorized!(G_ext_vectorized, g_eval, values, C0, ne, m, h, npg, EQoLG))
    C1 = LU_dec\(F_vectorized!(F_ext_vectorized, X, f_eval, values, x -> f(x, tau*0.5), ne, m, h, npg, EQoLG) + MK*C0 -
         G_vectorized!(G_ext_vectorized, g_eval, values, (C0 + C1_tiu)/2, ne, m, h, npg, EQoLG))
    
    # Pré-aloca memória para operações repetidas
    Cn = zeros(Float64, ne-1)
    G_C = zeros(Float64, ne-1)
    MK_C1 = zeros(Float64, ne-1)
    @views @simd for n in 2:N
        # Cálculo sem consumo extra de memória
        G_C .= 3*C1
        G_C .-= C0
        G_C ./= 2
        MK_C1 .= MK*C1
        ######################################
        Cn .= LU_dec\(F_vectorized!(F_ext_vectorized, X, f_eval, values, x -> f(x, tau*(n-0.5)), ne, m, h, npg, EQoLG) + MK*C1 -
            G_vectorized!(G_ext_vectorized, g_eval, values, G_C, ne, m, h, npg, EQoLG))
        # Solução aproximada
        C0 .= C1
        C1 .= Cn
    end
end samples=samples_sys evals=evals_sys
display(output_sys)

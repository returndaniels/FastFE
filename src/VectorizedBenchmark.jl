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
C0_ext = [C0_options(option, u0, a, ne, m, h, npg, EQoLG); 0]

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

output_G = benchmark(G_vectorized!, (G_ext_vectorized, g_eval, values, C0_ext, ne, m, h, npg, EQoLG), samples, evals=evals)
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

    # Pré-aloca mmeória para o lado direito do sistema linear
    B = zeros(Float64, ne-1)

    # Montando o sistema linear
    M = K_vectorized(ne, m, h, npg, 0., 1., 0., EQoLG)[1:ne-1, 1:ne-1]
    K = K_vectorized(ne, m, h, npg, alpha, beta, gamma, EQoLG)[1:ne-1, 1:ne-1]
    
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
    @views @simd for n in 2:N
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
    end
end samples=samples_sys evals=evals_sys
display(output_sys)

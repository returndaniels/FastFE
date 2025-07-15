module SerialErrorAnalysis

include("./SerialFiniteElements.jl")
include("./MDEDiscretization.jl")

using .SerialFiniteElements
using .MDEDiscretization
using SparseArrays
using LinearAlgebra
using Plots

export calculate_errors_serial

@doc raw"""
    calculate_errors_serial(tam::Int64, u::Function, u0::Function, f::Function, EQoLG_func::Function, 
                            K::SparseMatrixCSC{Float64, Int64}, C0_options::Vector{Float64}, α::Float64, β::Float64, γ::Float64, 
                            a::Float64, b::Float64, npg::Int64, option::Int64)

Calcula e plota a convergência do erro para diferentes discretizações usando uma abordagem serial.
"""
function calculate_errors_serial(tam::Int64, u::Function, u0::Function, f::Function, 
                                 α::Float64, β::Float64, γ::Float64,
                                 a::Float64, b::Float64, npg::Int64, option::Int64)
    erros = zeros(tam - 1)
    hs = zeros(tam - 1)
    for i in 2:tam
        # Define a discretização
        ne_er = 2^i
        N_er = ne_er
        m_er = ne_er - 1
        h_er = (b - a) / ne_er
        tau_er = h_er
        hs[i - 1] = h_er
        # Monta estrutura local global
        EQoLG = monta_EQ(ne_er)[monta_LG(ne_er)]
        # Estruturas externas usadas pelas funções
        x = h_er*(P .+ 1)/2 .+ a
        F_ext_serial = zeros(Float64, ne_er)
        G_ext_serial = zeros(Float64, ne_er)
        # Pré-aloca mmeória para o lado direito do sistema linear
        B = zeros(Float64, ne_er-1)
        # Montando o sistema linear
        M = K_serial(ne_er, m_er, h_er, npg, 0., 1., 0., EQoLG)[1:ne_er-1, 1:ne_er-1]
        K = K_serial(ne_er, m_er, h_er, npg, α, β, γ, EQoLG)[1:ne_er-1, 1:ne_er-1]
        MK = M/tau_er - K/2
        # Decomposição LU para resolução de múltiplos sistemas
        LU_dec = lu(M/tau_er + K/2)
        # Aproximação de U0
        C0_ext = [C0_options(option, u0, a, ne_er, m_er, h_er, npg, EQoLG); 0]
        C0 = @view C0_ext[1:ne_er-1]
        MKC0 = MK*C0
        # Aproximação de U1
        F_serial!(F_ext_serial, x, x -> f(x, tau_er*0.5), ne_er, m_er, h_er, npg, EQoLG)
        F = @view F_ext_serial[1:ne_er-1]
        G_serial!(G_ext_serial, C0_ext, ne_er, m_er, h_er, npg, EQoLG)
        G = @view G_ext_serial[1:ne_er-1]
        C1_ext_tiu = zeros(Float64, ne_er)
        B .= F
        B .+= MKC0
        B .-= G
        C1_ext_tiu[1:ne_er-1] .= LU_dec\B
        # Calcula C1
        G_serial!(G_ext_serial, (C0_ext + C1_ext_tiu)/2, ne_er, m_er, h_er, npg, EQoLG)
        G = @view G_ext_serial[1:ne_er-1]
        C1_ext = zeros(Float64, ne_er)
        B .= F
        B .+= MKC0
        B .-= G
        C1_ext[1:ne_er-1] = LU_dec\B
        # calcula os erros iniciais
        erro = zeros(N_er + 1)
        erro[1] = erro_serial(u0, x, ne_er, m_er, h_er, npg, C0_ext, EQoLG)
        erro[2] = erro_serial(x -> u(x, tau_er), x, ne_er, m_er, h_er, npg, C1_ext, EQoLG)
        # Pré-aloca memória para operações repetidas
        Cn_ext = zeros(Float64, ne_er)
        G_C_ext = zeros(Float64, ne_er)
        MKC1 = zeros(Float64, ne_er-1)
        B = zeros(Float64, ne_er-1)
        @views @simd for n in 2:N_er
            ############ Cálculo sem consumo extra de memória ############
            F_serial!(F_ext_serial, x, x -> f(x, tau_er*(n-0.5)), ne_er, m_er, h_er, npg, EQoLG)
            F = F_ext_serial[1:ne_er-1]
            MKC1 .= MK*(C1_ext[1:ne_er-1])
            G_C_ext[1:ne_er-1] .= 3*(C1_ext[1:ne_er-1])
            G_C_ext[1:ne_er-1] .-= C0_ext[1:ne_er-1]
            G_C_ext[1:ne_er-1] ./= 2
            G_serial!(G_ext_serial, G_C_ext, ne_er, m_er, h_er, npg, EQoLG)
            G = G_ext_serial[1:ne_er-1]
            B .= F
            B .+= MKC1
            B .-= G
            ##############################################################
            # Calcula o próximo coeficiente e atualiza
            Cn_ext[1:ne_er-1] .= LU_dec\B
            C0_ext[1:ne_er-1] .= C1_ext[1:ne_er-1]
            C1_ext[1:ne_er-1] .= Cn_ext[1:ne_er-1]
            erro[n + 1] = erro_serial(x -> u(x, tau_er * n), x, ne_er, m_er, h_er, npg, C1_ext, EQoLG)
        end
        erros[i-1] = maximum(erro)
    end
    return erros, hs
end

end # module

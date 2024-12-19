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
function calculate_errors_serial(tam::Int64, u::Function, u0::Function, f::Function, EQoLG_func::Function, 
                                 K::SparseMatrixCSC{Float64, Int64}, C0_options::Vector{Float64}, α::Float64, β::Float64, γ::Float64, 
                                 a::Float64, b::Float64, npg::Int64, option::Int64)
    erros = zeros(tam - 1)
    hs = zeros(tam - 1)

    for i in 2:tam
        # Define discretizations
        ne_er = 2^i
        N_er = ne_er
        m_er = ne_er - 1
        h_er = (b - a) / ne_er
        tau_er = h_er
        hs[i - 1] = h_er
        EQoLG = EQoLG_func(ne_er)

        # Initialize variables
        X = ((h_er / 2) * (P .+ 1) .+ a)' .+ range(a, step=h_er, stop=b - h_er)
        f_eval = Matrix{Float64}(undef, ne_er, npg)
        g_eval = Matrix{Float64}(undef, ne_er, npg)
        u_eval = Matrix{Float64}(undef, ne_er, npg)
        values = Matrix{Float64}(undef, ne_er, 2)
        F_ext_serial = zeros(Float64, ne_er)
        G_ext_serial = zeros(Float64, ne_er)
        x = h_er * (P .+ 1) / 2 .+ a

        # Compute linear system
        M = K_serial(ne_er, m_er, h_er, npg, 0., 1., 0., EQoLG)
        MK = M / tau_er - K / 2
        LU_dec = lu(M / tau_er + K / 2)

        # Approximate initial conditions
        C0 = C0_options(option, u0, a, ne_er, m_er, h_er, npg, EQoLG)
        C1_tiu = LU_dec \ (F_serial!(F_ext_serial, X, f_eval, values, x -> f(x, tau_er * 0.5), ne_er, m_er, h_er, npg, EQoLG) +
                            MK * C0 - G_serial!(G_ext_serial, g_eval, values, C0, ne_er, m_er, h_er, npg, EQoLG))
        C1 = LU_dec \ (F_serial!(F_ext_serial, X, f_eval, values, x -> f(x, tau_er * 0.5), ne_er, m_er, h_er, npg, EQoLG) +
                        MK * C0 - G_serial!(G_ext_serial, g_eval, values, (C0 + C1_tiu) / 2, ne_er, m_er, h_er, npg, EQoLG))

        # Calculate error
        erro = zeros(N_er + 1)
        erro[1] = erro_serial(u0, X, u_eval, ne_er, m_er, h_er, npg, C0, EQoLG)
        erro[2] = erro_serial(x -> u(x, tau_er), X, u_eval, ne_er, m_er, h_er, npg, C1, EQoLG)
        Cn_menos_2, Cn_menos_1 = C0, C1

        for n in 2:N_er
            temp = LU_dec \ (F_serial!(F_ext_serial, X, f_eval, values, x -> f(x, tau_er * (n - 0.5)), ne_er, m_er, h_er, npg, EQoLG) +
                              MK * Cn_menos_1 - G_serial!(G_ext_serial, g_eval, values, (3 * Cn_menos_1 - Cn_menos_2) / 2, ne_er, m_er, h_er, npg, EQoLG))
            Cn_menos_2, Cn_menos_1 = Cn_menos_1, temp
            erro[n + 1] = erro_serial(x -> u(x, tau_er * n), X, u_eval, ne_er, m_er, h_er, npg, Cn_menos_1, EQoLG)
        end

        erros[i - 1] = maximum(erro)
    end

    # Plot results
    # plot(hs, hs .^ 2, label="h^2", xlabel="h", title="Ordem de Convergência do Erro (Escala Log)", legend=:topright)
    # plot!(hs, erros, yscale=:log10, xscale=:log10, label="Erro")

    return erros
end

end # module

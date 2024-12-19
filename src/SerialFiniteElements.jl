module SerialFiniteElements

include("./MDEDiscretization.jl")

using LinearAlgebra
using SparseArrays
using .MDEDiscretization

# Funções de elementos finitos

@doc raw"""
    ($name)$(ne.value, m.value, h.value, npg.value, alpha.value, beta.value, gamma.value, EQoLG.value) -> SparseMatrixCSC

Calcula a matriz de rigidez global `K` a partir das matrizes locais de elementos finitos, utilizando integração numérica por pontos de Gauss.

# Argumentos
- `$ne.value`: Número de elementos finitos no domínio.
- `$m.value`: Número de pontos no intervalo.
- `$h.value`: Tamanho do passo no domínio espacial.
- `$npg.value`: Número de pontos de Gauss utilizados na integração.
- `$alpha.value`: Coeficiente do termo proporcional à derivada segunda na equação diferencial.
- `$beta.value`: Coeficiente do termo proporcional à derivada primeira.
- `$gamma.value`: Coeficiente do termo proporcional à função.
- `$EQoLG.value`: Matriz de conectividade que mapeia os nós locais para os globais.

# Retorno
Retorna uma matriz esparsa (`SparseMatrixCSC`) representando a matriz de rigidez global.
"""
function K_serial(ne::Int64, m::Int64, h::Float64, npg::Int64,
    alpha::Float64, beta::Float64, gamma::Float64, EQoLG::Matrix{Int64})::SparseMatrixCSC{Float64, Int64}
    @views begin
        K = spzeros(Float64, ne, ne)
        K_local = zeros(Float64, 2, 2)
        K_local[1, 1] = dot(W, 2*alpha*(dφ1P.*dφ1P)/h + gamma*(φ1P.*dφ1P) + beta*h*(φ1P.*φ1P)/2)
        K_local[1, 2] = dot(W, 2*alpha*(dφ1P.*dφ2P)/h + gamma*(φ1P.*dφ2P) + beta*h*(φ1P.*φ2P)/2)
        K_local[2, 1] = dot(W, 2*alpha*(dφ2P.*dφ1P)/h + gamma*(φ2P.*dφ1P) + beta*h*(φ2P.*φ1P)/2)
        K_local[2, 2] = dot(W, 2*alpha*(dφ2P.*dφ2P)/h + gamma*(φ2P.*dφ2P) + beta*h*(φ2P.*φ2P)/2)
        I, J, V = zeros(Int64, 4*ne), zeros(Int64, 4*ne), zeros(Float64, 4*ne)
        @simd for idx in 0:3
            i = div(idx, 2) + 1
            j = (idx % 2) + 1
            offset = idx*ne
            @simd for e in 1:ne
                pos = EQoLG[e, :]
                I[e + offset] = pos[i]
                J[e + offset] = pos[j]
                V[e + offset] = K_local[i, j]
            end
        end
        K = sparse(I, J, V, ne, ne)
        return K[1:(ne-1),1:(ne-1)]
    end
end

@doc raw"""
    ($name)$(F_ext_serial.value, x.value, f.value, ne.value, m.value, h.value, npg.value, EQoLG.value) -> Vector{Float64}

Calcula a contribuição das forças externas utilizando integração numérica por quadratura de Gauss.

# Argumentos
- `$F_ext_serial.value`: Vetor de forças externas a ser preenchido.
- `$x.value`: Pontos de Gauss no domínio.
- `$f.value`: Função que representa o termo fonte na equação diferencial.
- `$ne.value`: Número de elementos finitos.
- `$m.value`: Número de pontos no intervalo.
- `$h.value`: Tamanho do passo no domínio espacial.
- `$npg.value`: Número de pontos de Gauss.
- `$EQoLG.value`: Matriz de conectividade.

# Retorno
Retorna o vetor de forças externas calculado.
"""
function F_serial!(F_ext_serial::Vector{Float64}, x::Vector{Float64},
    f::Function, ne::Int64, m::Int64, h::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}
    @views begin
        fill!(F_ext_serial, 0.0)
        @simd for e in 1:ne
            offset = (e-1)*h
            i = EQoLG[e, 1]
            j = EQoLG[e, 2]
            @simd for g in 1:npg
                f_eval = h*f(x[g] + offset)/2
                F_ext_serial[i] += Wφ1P[g]*f_eval
                F_ext_serial[j] += Wφ2P[g]*f_eval
            end
        end
        return F_ext_serial[1:(ne-1)]
    end
end

@doc raw"""
    ($name)$(G_ext_serial.value, C.value, ne.value, m.value, h.value, npg.value, EQoLG.value) -> Vector{Float64}

Calcula a contribuição do termo não-linear associado à solução aproximada utilizando quadratura de Gauss.

# Argumentos
- `$G_ext_serial.value`: Vetor para armazenar o resultado.
- `$C.value`: Vetor de coeficientes da solução aproximada.
- `$ne.value`: Número de elementos finitos.
- `$m.value`: Número de pontos no intervalo.
- `$h.value`: Tamanho do passo no domínio espacial.
- `$npg.value`: Número de pontos de Gauss.
- `$EQoLG.value`: Matriz de conectividade.

# Retorno
Retorna o vetor `G_ext_serial` com os valores calculados.
"""
function G_serial!(G_ext_serial::Vector{Float64}, C::Vector{Float64},
    ne::Int64, m::Int64, h::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}
    @views begin
        fill!(G_ext_serial, 0.0)
        C_ext = [C;0]
        @simd for e in 1:ne
            i = EQoLG[e, 1]
            j = EQoLG[e, 2]
            c1 = C_ext[i]
            c2 = C_ext[j]
            @simd for p in 1:npg
                g_eval = h*g(c1*φ1P[p] + c2*φ2P[p])/2
                G_ext_serial[i] += Wφ1P[p]*g_eval
                G_ext_serial[j] += Wφ2P[p]*g_eval
            end
        end
        return G_ext_serial[1:(ne-1)]
    end
end

@doc raw"""
    ($name)$(u.value, x.value, ne.value, m.value, h.value, npg.value, C.value, EQoLG.value) -> Float64

Calcula o erro `L_2` entre a solução exata `u` e a solução aproximada representada pelos coeficientes `C`.

# Argumentos
- `$u.value: Função representando a solução exata.
- `$x.value`: Pontos de Gauss no domínio.
- `$ne.value`: Número de elementos finitos.
- `$m.value`: Número de pontos no intervalo.
- `$h.value`: Tamanho do passo no domínio espacial.
- `$npg.value`: Número de pontos de Gauss.
- `$C.value`: Vetor de coeficientes da solução aproximada.
- `$EQoLG.value`: Matriz de conectividade.

# Retorno
Retorna o valor do erro `L_2`.
"""
function erro_serial(u::Function, x::Vector{Float64}, ne::Int64, m::Int64, h::Float64, npg::Int64,
    C::Vector{Float64}, EQoLG::Matrix{Int64})::Float64
    @views begin
        C_ext = [C;0]
        sum_er = 0
        @simd for e in 1:ne
            C1 = C_ext[EQoLG[e, 1]]
            C2 = C_ext[EQoLG[e, 2]]
            offset = (e-1)*h
            @simd for g in 1:npg
                sum_er += W[g]*((u(x[g] + offset) - (C1*φ1P[g] + C2*φ2P[g])).^2)
            end
        end
        return (h*sum_er/2)^0.5
    end
end

@doc raw"""
    ($name)$(op.value, u0.value, a.value, ne.value, m.value, h.value, alpha.value, beta.value, gamma.value, npg.value, EQoLG.value) -> Vector{Float64}

Inicializa o vetor de coeficientes `C_0` com base na escolha de método.

# Argumentos
- `op::Int64`: Tipo de inicialização:
    - `1`: Interpolação da condição inicial.
    - `2`: Projeção `L_2`.
    - `3`: Projeção `H`.
    - `4`: Projeção com operador de rigidez.
- `$u0.value`: Função da condição inicial.
- `$a.value`: Limite inferior do domínio.
- `@ne.value`: Número de elementos finitos.
- `@m.value`: Número de pontos no intervalo.
- `$h.value`: Tamanho do passo no domínio.
- `$alpha.value`: Coeficiente da equação diferencial.
- `$beta.value`: Coeficiente da equação diferencial.
- `$gamma.value`: Coeficiente da equação diferencial.
- `$npg.value`: Número de pontos de Gauss.
- `$EQoLG.value`: Matriz de conectividade.

# Retorno
Retorna o vetor de coeficientes `C_0`.
"""
function C0_options(op::Int64, u0::Function, a::Float64, ne::Int64, m::Int64, h::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}
    # Interpolante de u0
    if(op == 1)
        @views begin
            C0 = zeros(Float64, ne-1)
            @simd for i in 1:(ne-1)
                C0[i] = u0(a + i*h)
            end
            return C0
        end
    # Projeção L2 de u0
    elseif(op == 2)
        M = K_QG(ne, m, h, 0, 1, 0, npg)
        return M\F_QG(u0, a, ne, m, h, npg, EQoLG)
    # Projeção H de u0
    elseif(op == 3)
        # Gerando o vetor de termos independentes (notação local)
        vet_local = zeros(2)
        # Preenchendo o vetor global
        vet = zeros(ne)
        for e in 1:ne
            vet_local[1] = (u0(a + (e-1)*h) - u0(a + e*h))/h
            vet_local[2] = (u0(a + e*h) - u0(a + (e-1)*h))/h
            for a in 1:2
                i = EQoLG[a, e]
                vet[i] += vet_local[a]
            end 
        end
        D = K_QG(ne, m, h, 1, 0, 0, npg)
        return D\vet[1:(ne-1)]
    # Operador k(u, v) como projeção de u0
    elseif(op == 4)
        phi1(x) = (1-x)/2
        phi2(x) = (1+x)/2
        # Gerando o vetor de termos independentes (notação local)
        vet_local = zeros(2)
        # Quadratura gaussiana
        X, W = gausslegendre(npg)
        # Preenchendo o vetor global
        vet = zeros(ne)
        for e in 1:ne
            vet_local[1] = alpha*(u0(a+(e-1)*h)-u0(a+e*h))/h + gamma*W'*(du0.(h*(X .+ 1)/2 .+ (a+(e-1)*h)).*phi1.(X))*h/2 +
                           beta*W'*(u0.(h*(X .+ 1)/2 .+ (a+(e-1)*h)).*(1 .- X)/2)*h/2
            vet_local[2] = alpha*(u0(a+e*h)-u0(a+(e-1)*h))/h + gamma*W'*(du0.(h*(X .+ 1)/2 .+ (a+(e-1)*h)).*phi2.(X))*h/2 +
                           beta*W'*(u0.(h*(X .+ 1)/2 .+ (a+(e-1)*h)).*(1 .+ X)/2)*h/2
            for a in 1:2
                i = EQoLG[a, e]
                vet[i] += vet_local[a]
            end
        end
        K = K_QG(ne, m, h, alpha, beta, gamma, npg)
        return K\vet[1:(ne-1)]
    else
        error("Opção inválida!")
    end
end

export K_serial, F_serial!, G_serial!, erro_serial, C0_options

end # módulo FiniteElements

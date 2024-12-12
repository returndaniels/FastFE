module VectorizedFiniteElements

using SparseArrays

"""
    K_vectorized(ne::Int64, m::Int64, h::Float64, npg::Int64,
                 alpha::Float64, beta::Float64, gamma::Float64, EQoLG::Matrix{Int64})::SparseMatrixCSC{Float64, Int64}

Calcula a matriz global de rigidez \( K \) para um sistema de elementos finitos. A matriz é uma soma ponderada das contribuições de cada elemento, utilizando pontos de Gauss e pesos para a integração numérica. Esta versão da função faz uso de vetorização para otimizar o processo de cálculo e melhorar a eficiência computacional.

# Argumentos
- `ne::Int64`: Número de elementos finitos no domínio.
- `m::Int64`: Número de pontos no intervalo.
- `h::Float64`: Tamanho do passo no espaço.
- `npg::Int64`: Número de pontos de Gauss utilizados na integração.
- `alpha::Float64`: Coeficiente \(\alpha\) na equação diferencial.
- `beta::Float64`: Coeficiente \(\beta\) na equação diferencial.
- `gamma::Float64`: Coeficiente \(\gamma\) na equação diferencial.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

# Retorno
Retorna uma matriz esparsa do tipo `SparseMatrixCSC{Float64, Int64}`, que representa a matriz de rigidez global \( K \).
"""
function K_vectorized(ne::Int64, m::Int64, h::Float64, npg::Int64,
    alpha::Float64, beta::Float64, gamma::Float64, EQoLG::Matrix{Int64})::SparseMatrixCSC{Float64, Int64}
    @views begin
        i = EQoLG[:, 1]
        j = EQoLG[:, 2]
        lin_idx = vcat(i, i, j, j)
        col_idx = vcat(i, j, i, j)
        k_values = vcat(fill(dot(W, 2*alpha*(dφ1P.*dφ1P)/h + gamma*(φ1P.*dφ1P) + beta*h*(φ1P.*φ1P)/2), ne),
            fill(dot(W, 2*alpha*(dφ1P.*dφ2P)/h + gamma*(φ1P.*dφ2P) + beta*h*(φ1P.*φ2P)/2), ne),
            fill(dot(W, 2*alpha*(dφ2P.*dφ1P)/h + gamma*(φ2P.*dφ1P) + beta*h*(φ2P.*φ1P)/2), ne),
            fill(dot(W, 2*alpha*(dφ2P.*dφ2P)/h + gamma*(φ2P.*dφ2P) + beta*h*(φ2P.*φ2P)/2), ne))
        K = sparse(lin_idx, col_idx, k_values, ne, ne)
        return K[1:(ne-1),1:(ne-1)]
    end
end

"""
    F_vectorized!(F_ext_vectorized::Vector{Float64}, X::Matrix{Float64}, f_eval::Matrix{Float64}, values::Matrix{Float64},
                  f::Function, ne::Int64, m::Int64, h::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}

Calcula o vetor de forças externas \( F \) por meio de uma quadratura de Gauss, utilizando vetorização para otimizar o processo de cálculo. A função é responsável por computar a contribuição de cada ponto de Gauss para o vetor de forças externas, com base na função fornecida \( f \).

# Argumentos
- `F_ext_vectorized::Vector{Float64}`: Vetor de forças externas que será preenchido pela função.
- `X::Matrix{Float64}`: Matriz contendo as coordenadas dos pontos de Gauss no espaço.
- `f_eval::Matrix{Float64}`: Matriz contendo os valores da função \( f \) nos pontos de Gauss.
- `values::Matrix{Float64}`: Matriz de valores utilizados na avaliação da função \( f \) e no cálculo do vetor de forças.
- `f::Function`: Função a ser integrada numericamente sobre os elementos finitos.
- `ne::Int64`: Número de elementos finitos no domínio.
- `m::Int64`: Número de pontos no intervalo.
- `h::Float64`: Tamanho do passo no espaço.
- `npg::Int64`: Número de pontos de Gauss.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

# Retorno
Retorna o vetor `F_ext_vectorized` preenchido com as contribuições das forças externas, calculadas pela quadratura de Gauss.
"""
function F_vectorized!(F_ext_vectorized::Vector{Float64}, X::Matrix{Float64}, f_eval::Matrix{Float64}, values::Matrix{Float64},
    f::Function, ne::Int64, m::Int64, h::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}
    @views begin
        fill!(F_ext_vectorized, 0.0)
        @. f_eval = f(X)
        mul!(values, f_eval, WφP)
        @simd for i in 1:2
        F_ext_vectorized[EQoLG[:, i]] .+= values[:, i]
        end
            return (h/2)*F_ext_vectorized[1:ne-1]
    end
end

"""
    G_vectorized!(G_ext_vectorized::Vector{Float64}, g_eval::Matrix{Float64}, values::Matrix{Float64},
                  C::Vector{Float64}, ne::Int64, m::Int64, h::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}

Calcula o vetor de forças estendido \( G \) associado à função \( g \), utilizando vetorização para melhorar a eficiência dos cálculos. A função realiza a integração numérica utilizando pontos de Gauss e pesos, com o vetor de coeficientes \( C \) representando a solução aproximada.

# Argumentos
- `G_ext_vectorized::Vector{Float64}`: Vetor de forças estendidas que será preenchido pela função.
- `g_eval::Matrix{Float64}`: Matriz contendo os valores da função \( g \) nos pontos de Gauss.
- `values::Matrix{Float64}`: Matriz de valores utilizados na avaliação da função \( g \) e no cálculo do vetor de forças.
- `C::Vector{Float64}`: Vetor de coeficientes da solução aproximada \( C \).
- `ne::Int64`: Número de elementos finitos no domínio.
- `m::Int64`: Número de pontos no intervalo.
- `h::Float64`: Tamanho do passo no espaço.
- `npg::Int64`: Número de pontos de Gauss.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

# Retorno
Retorna o vetor `G_ext_vectorized` preenchido com as contribuições das forças associadas à função \( g \), calculadas pela quadratura de Gauss.
"""
function G_vectorized!(G_ext_vectorized::Vector{Float64}, g_eval::Matrix{Float64}, values::Matrix{Float64},
    C::Vector{Float64}, ne::Int64, m::Int64, h::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}
    @views begin
        fill!(G_ext_vectorized, 0.0)
        C_ext = [C;0]
        @. g_eval = g((φ1P')*C_ext[EQoLG[:, 1]] + (φ2P')*C_ext[EQoLG[:, 2]])
        mul!(values, g_eval, WφP)
        @simd for i in 1:2
            G_ext_vectorized[EQoLG[:, i]] .+= values[:, i]
        end
        return (h/2)*G_ext_vectorized[1:(ne-1)]
    end
end

# Exportando as funções para que sejam acessíveis externamente
export K_vectorized, F_vectorized!, G_vectorized!

end

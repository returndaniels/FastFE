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
    # (função já definida anteriormente)
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
    # (função já definida anteriormente)
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
    # (função já definida anteriormente)
end

# Exportando as funções para que sejam acessíveis externamente
export K_vectorized, F_vectorized!, G_vectorized!

end

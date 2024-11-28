module FiniteElements

using LinearAlgebra
using SparseArrays

# Funções de elementos finitos

"""
    K_serial(ne::Int64, m::Int64, h::Float64, npg::Int64,
             alpha::Float64, beta::Float64, gamma::Float64, EQoLG::Matrix{Int64})::SparseMatrixCSC{Float64, Int64}

Calcula a matriz de rigidez global K a partir das matrizes locais, utilizando pontos de Gauss e pesos para a integração numérica.
Essa função é implementada de forma sequencial.

# Argumentos
- `ne::Int64`: Número de elementos finitos no domínio de `u`.
- `m::Int64`: Número de pontos no intervalo `(a, b)`.
- `h::Float64`: Tamanho do passo no espaço.
- `npg::Int64`: Número de pontos de Gauss.
- `alpha::Float64`: Coeficiente alpha na equação diferencial.
- `beta::Float64`: Coeficiente beta na equação diferencial.
- `gamma::Float64`: Coeficiente gamma na equação diferencial.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

# Retorno
Retorna uma matriz esparsa de tipo `SparseMatrixCSC{Float64, Int64}` com a matriz de rigidez global.
"""
function K_serial(ne::Int64, m::Int64, h::Float64, npg::Int64,
                  alpha::Float64, beta::Float64, gamma::Float64, EQoLG::Matrix{Int64})::SparseMatrixCSC{Float64, Int64}
end

"""
    F_serial!(F_ext_serial::Vector{Float64}, x::Vector{Float64},
              f::Function, ne::Int64, m::Int64, h::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}

Calcula a quadratura gaussiana do vetor de forças externas `F_ext` associada à equação diferencial.

# Argumentos
- `F_ext_serial::Vector{Float64}`: Vetor de forças externas estendido a ser preenchido.
- `x::Vector{Float64}`: Vetor de pontos de Gauss no domínio de `f`.
- `f::Function`: Função que representa a equação diferencial.
- `ne::Int64`: Número de elementos finitos no intervalo.
- `m::Int64`: Número de pontos no intervalo `(a, b)`.
- `h::Float64`: Tamanho do passo no espaço.
- `npg::Int64`: Número de pontos de Gauss.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

# Retorno
Retorna o vetor `F_ext_serial` com os valores calculados pela quadratura gaussiana.
"""

function F_serial!(F_ext_serial::Vector{Float64}, x::Vector{Float64},
                   f::Function, ne::Int64, m::Int64, h::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}
end

"""
    G_serial!(G_ext_serial::Vector{Float64}, C::Vector{Float64},
              ne::Int64, m::Int64, h::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}

Calcula a quadratura gaussiana associada ao vetor `G_ext` considerando os coeficientes `C` da solução aproximada.

# Argumentos
- `G_ext_serial::Vector{Float64}`: Vetor de forças estendidas a ser preenchido.
- `C::Vector{Float64}`: Vetor de coeficientes da solução aproximada `u`.
- `ne::Int64`: Número de elementos finitos no intervalo.
- `m::Int64`: Número de pontos no intervalo `(a, b)`.
- `h::Float64`: Passo no domínio espacial.
- `npg::Int64`: Número de pontos de Gauss.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

# Retorno
Retorna o vetor `G_ext_serial` preenchido com os valores calculados pela quadratura gaussiana.
"""
function G_serial!(G_ext_serial::Vector{Float64}, C::Vector{Float64},
                   ne::Int64, m::Int64, h::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}
end

"""
    erro_QG(u::Function, a::Float64, ne::Int64, m::Int64, h::Float64, 
            pg::Int64, C::Vector{Float64}, EQoLG::Matrix{Int64})::Float64

Calcula o erro integral L2 entre a solução exata `u` e a solução aproximada representada pelos coeficientes `C`.

# Argumentos
- `u::Function`: Função solução exata.
- `a::Float64`: Limite inferior do domínio de `u`.
- `ne::Int64`: Número de elementos finitos no domínio de `u`.
- `m::Int64`: Número de pontos no intervalo `(a, b)`.
- `h::Float64`: Tamanho do passo no espaço.
- `pg::Int64`: Número de pontos de Gauss.
- `C::Vector{Float64}`: Vetor com os coeficientes da solução aproximada.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

# Retorno
Retorna o erro L2 entre a solução exata `u` e a solução aproximada.
"""
function erro_QG(u, a, ne, m, h, pg, C, EQoLG)
end

"""
    C0_options(op::Int64, u0::Function, a::Float64, ne::Int64, 
                m::Int64, h::Float64, alpha::Float64, beta::Float64, 
                gamma::Float64, npg::Int64, EQoLG::Matrix{Int64})::Vector{Float64}

Calcula o vetor `C0` de coeficientes iniciais, dado o tipo de inicialização especificado.

# Argumentos
- `op::Int64`: Opção de inicialização do vetor `C0`.
    - `1`: Interpolação da condição inicial `u0`.
    - `2`: Projeção L2 de `u0`.
    - `3`: Projeção H de `u0`.
    - `4`: Operador `k(u, v)` como projeção de `u0`.
- `u0::Function`: Função que representa a condição inicial no tempo.
- `a::Float64`: Limite inferior do domínio de `u`.
- `ne::Int64`: Número de elementos finitos no domínio espacial.
- `m::Int64`: Número de pontos no intervalo `(a, b)`.
- `h::Float64`: Tamanho do passo no espaço.
- `alpha::Float64`: Coeficiente alpha na equação diferencial.
- `beta::Float64`: Coeficiente beta na equação diferencial.
- `gamma::Float64`: Coeficiente gamma na equação diferencial.
- `npg::Int64`: Número de pontos de Gauss.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

# Retorno
Retorna o vetor `C0` com os coeficientes iniciais, de acordo com a opção selecionada.
"""
function C0_options(op, u0, a, ne, m, h, alpha, beta, gamma, npg, EQoLG)
end

export K_serial, F_serial!, G_serial!, erro_QG, C0_options

end # módulo FiniteElements


# Discretization of a Modern Differential Equations
module VectorizedSolution

# Importa o pacote necessário para obter os pontos de Gauss e pesos
using Quadrature

# Parâmetros da EDO
const beta = 1.0           # Beta da EDO
const alpha = 1.0          # Alpha da EDO
const gamma = 1.0          # Gamma da EDO

# Intervalos espaciais e temporais
const a = 0.0              # Limite inferior do intervalo espacial
const b = 1.0              # Limite superior do intervalo espacial
const T = 1.0              # Limite superior do intervalo temporal
const npg = 5              # Número de pontos de Gauss

# Discretização no espaço
const ne = 2^6             # Número de intervalos no espaço
const m = ne - 1           # Número de pontos no espaço (internos)
const h = (b - a) / ne     # Passo no espaço

# Discretização no tempo
const tau = h              # Passo no tempo
const N = trunc(Int, T / tau)  # Número de intervalos no tempo

# Obtém os pontos de Gauss e seus pesos
const P, W = gausslegendre(npg)

# Funções de forma e suas derivadas no elemento padrão
const φ1(ξ) = (1 - ξ) / 2
const φ2(ξ) = (1 + ξ) / 2
const dφ1(ξ) = -0.5
const dφ2(ξ) = 0.5

# Avaliações de funções de forma e suas derivadas nos pontos de Gauss
const φ1P = φ1.(P)
const φ2P = φ2.(P)
const dφ1P = dφ1.(P)
const dφ2P = dφ2.(P)

# Avaliações de funções de forma multiplicadas pelos pesos nos pontos de Gauss
const Wφ1P = W.*φ1.(P)
const Wφ2P = W.*φ2.(P)
##################### Funções #####################

"""
    u(x::Float64, t::Float64) -> Float64

Função solução da EDO que avalia a solução u no ponto (x, t).

# Argumentos
- `x::Float64`: Ponto no domínio espacial de u.
- `t::Float64`: Ponto no domínio temporal de u.

# Retorno
- Valor de u avaliado no ponto (x, t).
"""
function u(x::Float64, t::Float64) # Função solução da EDO
    return sin(pi*x)*exp(-t)/(pi^2)
end

"""
    u0(x::Float64) -> Float64

Função inicial que avalia u0 no ponto x.

# Argumentos
- `x::Float64`: Ponto no domínio espacial de u.

# Retorno
- Valor de u0 avaliado no ponto x.
"""
function u0(x::Float64)
    return sin(pi*x)/(pi^2)
end

"""
    du0(x::Float64) -> Float64

Função que calcula a derivada de u0 no ponto x.

# Argumentos
- `x::Float64`: Ponto no domínio espacial de u.

# Retorno
- Valor da derivada de u0 no ponto x.
"""
function du0(x::Float64)
    return cos(pi*x)/pi
end

"""
    g(s::Float64) -> Float64

Função no domínio de g.

# Argumentos
- `s::Float64`: Ponto no domínio de g.

# Retorno
- Valor de g no ponto s.
"""
function g(s::Float64)
    return s^3 - s
end

"""
    f(x::Float64, t::Float64) -> Float64

Função f que avalia uma expressão no ponto (x, t).

# Argumentos
- `x::Float64`: Ponto no domínio espacial de f.
- `t::Float64`: Ponto no domínio temporal de f.

# Retorno
- Valor de f avaliado no ponto (x, t).
"""
function f(x::Float64, t::Float64)
    func = ((alpha*(pi^2) + beta - 1)*sin(pi*x) + gamma*pi*cos(pi*x))*exp(-t)/(pi^2)
    return func + g(u(x, t))
end

"""
    monta_LG(ne::Int64) -> Matrix{Int64}

Função que monta a matriz de conectividade do elemento finito.

# Argumentos
- `ne::Int64`: Número de elementos finitos.

# Retorno
- Matriz de conectividade.
"""
function monta_LG(ne::Int64)::Matrix{Int64}
    return hcat(1:ne, 2:ne+1) # Cocatena duas linhas horizontais [1, 2, ..., ne] e [2, 3, ..., ne+1]
end

"""
    monta_EQ(ne::Int64) -> Vector{Int64}

Função que monta o vetor de numeração das funções phi.

# Argumentos
- `ne::Int64`: Número de elementos finitos.

# Retorno
- Vetor com a numeração das funções phi.
"""
function monta_EQ(ne::Int64)::Vector{Int64}
    return vcat(ne, 1:(ne-1), ne) # Concatena o número de elementos na esquerda e na direita do vetor [1, 2,..., ne-1]
end

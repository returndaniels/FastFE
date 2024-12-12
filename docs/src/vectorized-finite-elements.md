# Módulo `VectorizedFiniteElements`

Este módulo contém funções otimizadas para o cálculo de matrizes de rigidez e vetores de forças externas utilizando o método dos elementos finitos (FEM), com o uso de vetorização para melhorar a eficiência computacional.

## Funções

### `K_vectorized`

Calcula a matriz global de rigidez \( K \) para um sistema de elementos finitos utilizando vetorização.

#### Argumentos

- `ne::Int64`: Número de elementos finitos no domínio.
- `m::Int64`: Número de pontos no intervalo.
- `h::Float64`: Tamanho do passo no espaço.
- `npg::Int64`: Número de pontos de Gauss utilizados na integração.
- `alpha::Float64`: Coeficiente \(\alpha\) na equação diferencial.
- `beta::Float64`: Coeficiente \(\beta\) na equação diferencial.
- `gamma::Float64`: Coeficiente \(\gamma\) na equação diferencial.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

#### Retorno

Retorna uma matriz esparsa do tipo `SparseMatrixCSC{Float64, Int64}`, que representa a matriz de rigidez global \( K \).

#### Exemplo de uso

```julia
using SparseArrays
using VectorizedFiniteElements

# Definir parâmetros
ne = 10
m = 5
h = 0.1
npg = 3
alpha = 1.0
beta = 2.0
gamma = 3.0
EQoLG = Matrix{Int64}(rand(10, 2))

# Calcular a matriz de rigidez
K = K_vectorized(ne, m, h, npg, alpha, beta, gamma, EQoLG)
println(K)
```

### `F_vectorized!`

Calcula o vetor de forças externas \( F \) por meio de uma quadratura de Gauss, utilizando vetorização.

#### Argumentos

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

#### Retorno

Retorna o vetor `F_ext_vectorized` preenchido com as contribuições das forças externas, calculadas pela quadratura de Gauss.

#### Exemplo de uso

```julia
using LinearAlgebra
using VectorizedFiniteElements

# Definir parâmetros
ne = 10
m = 5
h = 0.1
npg = 3
f(x) = sin(x)  # Função f a ser integrada
X = rand(ne, m)
f_eval = zeros(ne, m)
values = zeros(ne, m)
EQoLG = Matrix{Int64}(rand(ne, 2))
F_ext = zeros(ne)

# Calcular o vetor de forças externas
F = F_vectorized!(F_ext, X, f_eval, values, f, ne, m, h, npg, EQoLG)
println(F)
```

### `G_vectorized!`

Calcula o vetor de forças estendido \( G \) associado à função \( g \), utilizando vetorização.

#### Argumentos

- `G_ext_vectorized::Vector{Float64}`: Vetor de forças estendidas que será preenchido pela função.
- `g_eval::Matrix{Float64}`: Matriz contendo os valores da função \( g \) nos pontos de Gauss.
- `values::Matrix{Float64}`: Matriz de valores utilizados na avaliação da função \( g \) e no cálculo do vetor de forças.
- `C::Vector{Float64}`: Vetor de coeficientes da solução aproximada \( C \).
- `ne::Int64`: Número de elementos finitos no domínio.
- `m::Int64`: Número de pontos no intervalo.
- `h::Float64`: Tamanho do passo no espaço.
- `npg::Int64`: Número de pontos de Gauss.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

#### Retorno

Retorna o vetor `G_ext_vectorized` preenchido com as contribuições das forças associadas à função \( g \), calculadas pela quadratura de Gauss.

#### Exemplo de uso

```julia
using LinearAlgebra
using VectorizedFiniteElements

# Definir parâmetros
ne = 10
m = 5
h = 0.1
npg = 3
g(x) = cos(x)  # Função g a ser integrada
C = rand(ne)
g_eval = zeros(ne, m)
values = zeros(ne, m)
EQoLG = Matrix{Int64}(rand(ne, 2))
G_ext = zeros(ne)

# Calcular o vetor de forças estendidas
G = G_vectorized!(G_ext, g_eval, values, C, ne, m, h, npg, EQoLG)
println(G)
```

## Exportação das funções

As funções `K_vectorized`, `F_vectorized!` e `G_vectorized!` são exportadas pelo módulo para que possam ser usadas externamente:

```julia
export K_vectorized, F_vectorized!, G_vectorized!
```

## Requisitos

Este módulo depende de:

- `SparseArrays` para o cálculo de matrizes esparsas.
- `LinearAlgebra` para operações matriciais.

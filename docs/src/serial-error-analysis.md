# SerialErrorAnalysis.jl

Este módulo fornece uma abordagem serial para o cálculo da análise de erro em discretizações de elementos finitos. Ele inclui funções para calcular os erros nas aproximações das soluções e para plotar a convergência do erro para diferentes discretizações.

## Funções

### `calculate_errors_serial(tam::Int64, u::Function, u0::Function, f::Function, EQoLG_func::Function,

                             K_func::Function, C0_options::Function, α::Float64, β::Float64, γ::Float64,
                             a::Float64, b::Float64, npg::Int64, option::Int64)`

Calcula e plota a convergência do erro para diferentes discretizações utilizando uma abordagem serial.

#### Argumentos:

- `tam::Int64`: O número máximo de elementos na discretização.
- `u::Function`: A função que representa a solução exata.
- `u0::Function`: A função que representa a condição inicial.
- `f::Function`: A função do termo fonte.
- `EQoLG_func::Function`: Uma função para gerar a matriz de conectividade dos elementos.
- `K_func::Function`: Uma função para gerar a matriz de rigidez.
- `C0_options::Function`: Uma função para gerar o vetor inicial de coeficientes.
- `α::Float64`: Um parâmetro para a matriz de rigidez.
- `β::Float64`: Um parâmetro para a matriz de rigidez.
- `γ::Float64`: Um parâmetro para a matriz de rigidez.
- `a::Float64`: A fronteira esquerda do domínio.
- `b::Float64`: A fronteira direita do domínio.
- `npg::Int64`: O número de pontos de Gauss por elemento.
- `option::Int64`: Uma flag para selecionar diferentes opções de coeficientes.

#### Retorna:

- `Vector{Float64}`: Um vetor com os valores de erro para cada discretização.

#### Exemplo de uso:

```julia
# Defina as funções necessárias
u_exact(x) = sin(pi * x)
u0(x) = 0.0
f(x, t) = 0.0
EQoLG_func(ne) = [1 2; 2 3]
K_func(ne, m, h, npg, α, β, γ, EQoLG) = eye(ne)  # Exemplo simples de matriz de rigidez
C0_options(option, u0, a, ne, m, h, npg, EQoLG) = [1.0, 0.5]

# Defina os parâmetros
tam = 10
α = 1.0
β = 0.0
γ = 0.0
a = 0.0
b = 1.0
npg = 2
option = 1

# Calcule os erros
erros = calculate_errors_serial(tam, u_exact, u0, f, EQoLG_func, K_func, C0_options, α, β, γ, a, b, npg, option)
println("Erros: $erros")
```

---

## Dependências

- `SerialFiniteElements`: Utilizado para operações de elementos finitos serializados.
- `MDEDiscretization`: Fornece esquemas de discretização.
- `LinearAlgebra`: Oferece utilitários de álgebra linear.
- `Plots`: Usado para plotar a convergência do erro.

## Plotando a Convergência do Erro

A função `calculate_errors_serial` retorna os erros para diferentes discretizações, que podem ser plotados para visualizar a convergência do erro. O gráfico pode ser personalizado utilizando o pacote `Plots`.

```julia
using Plots

# Exemplo de plotagem da convergência do erro
plot(erros, label="Erro", xlabel="Número de Discretizações", ylabel="Erro", title="Convergência do Erro")
```

---

## Notas

- O módulo utiliza operações serializadas para calcular os erros, otimizando o desempenho para problemas menores ou onde a paralelização não é necessária.
- A convergência do erro é calculada para diferentes refinamentos da malha, e os resultados podem ser usados para avaliar a ordem de convergência.

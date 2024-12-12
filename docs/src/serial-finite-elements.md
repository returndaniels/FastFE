# Documentação do módulo `SerialFiniteElements`

Este módulo fornece funções para resolver problemas de elementos finitos utilizando métodos sequenciais. As funções implementadas incluem o cálculo da matriz de rigidez global, o vetor de forças externas, e a solução para o erro integral L2, entre outras operações de discretização e integração numérica.

## Funções

### `K_serial`

Calcula a matriz de rigidez global `K` a partir das matrizes locais, utilizando pontos de Gauss e pesos para a integração numérica. Esta função é implementada de forma sequencial.

#### Argumentos:

- `ne::Int64`: Número de elementos finitos no domínio de `u`.
- `m::Int64`: Número de pontos no intervalo `(a, b)`.
- `h::Float64`: Tamanho do passo no espaço.
- `npg::Int64`: Número de pontos de Gauss.
- `alpha::Float64`: Coeficiente alpha na equação diferencial.
- `beta::Float64`: Coeficiente beta na equação diferencial.
- `gamma::Float64`: Coeficiente gamma na equação diferencial.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

#### Retorno:

Retorna uma matriz esparsa de tipo `SparseMatrixCSC{Float64, Int64}` com a matriz de rigidez global.

#### Exemplo de Uso:

```julia
# Parâmetros
ne = 10
m = 100
h = 0.1
npg = 3
alpha = 1.0
beta = 0.5
gamma = 0.3
EQoLG = rand(10, 2)  # Exemplo de mapeamento local-global

# Chamada da função
K = K_serial(ne, m, h, npg, alpha, beta, gamma, EQoLG)
println(K)
```

---

### `F_serial!`

Calcula a quadratura gaussiana do vetor de forças externas `F_ext` associada à equação diferencial.

#### Argumentos:

- `F_ext_serial::Vector{Float64}`: Vetor de forças externas estendido a ser preenchido.
- `x::Vector{Float64}`: Vetor de pontos de Gauss no domínio de `f`.
- `f::Function`: Função que representa a equação diferencial.
- `ne::Int64`: Número de elementos finitos no intervalo.
- `m::Int64`: Número de pontos no intervalo `(a, b)`.
- `h::Float64`: Tamanho do passo no espaço.
- `npg::Int64`: Número de pontos de Gauss.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

#### Retorno:

Retorna o vetor `F_ext_serial` com os valores calculados pela quadratura gaussiana.

#### Exemplo de Uso:

```julia
# Parâmetros
F_ext_serial = zeros(Float64, ne)
x = linspace(0, 1, m)
f(x) = x^2  # Exemplo de função f(x)

# Chamada da função
F_ext_serial = F_serial!(F_ext_serial, x, f, ne, m, h, npg, EQoLG)
println(F_ext_serial)
```

---

### `G_serial!`

Calcula a quadratura gaussiana associada ao vetor `G_ext` considerando os coeficientes `C` da solução aproximada.

#### Argumentos:

- `G_ext_serial::Vector{Float64}`: Vetor de forças estendidas a ser preenchido.
- `C::Vector{Float64}`: Vetor de coeficientes da solução aproximada `u`.
- `ne::Int64`: Número de elementos finitos no intervalo.
- `m::Int64`: Número de pontos no intervalo `(a, b)`.
- `h::Float64`: Passo no domínio espacial.
- `npg::Int64`: Número de pontos de Gauss.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

#### Retorno:

Retorna o vetor `G_ext_serial` preenchido com os valores calculados pela quadratura gaussiana.

#### Exemplo de Uso:

```julia
# Parâmetros
G_ext_serial = zeros(Float64, ne)
C = rand(Float64, ne)  # Coeficientes da solução

# Chamada da função
G_ext_serial = G_serial!(G_ext_serial, C, ne, m, h, npg, EQoLG)
println(G_ext_serial)
```

---

### `erro_serial`

Calcula o erro integral L2 entre a solução exata `u` e a solução aproximada representada pelos coeficientes `C`.

#### Argumentos:

- `u::Function`: Função solução exata.
- `x::Vector{Float64}`: Vetor de pontos de Gauss no domínio de `f`.
- `ne::Int64`: Número de elementos finitos no domínio de `u`.
- `m::Int64`: Número de pontos no intervalo `(a, b)`.
- `h::Float64`: Tamanho do passo no espaço.
- `npg::Int64`: Número de pontos de Gauss.
- `C::Vector{Float64}`: Vetor com os coeficientes da solução aproximada.
- `EQoLG::Matrix{Int64}`: Matriz de mapeamento de elementos locais para globais.

#### Retorno:

Retorna o erro L2 entre a solução exata `u` e a solução aproximada.

#### Exemplo de Uso:

```julia
# Função solução exata
u(x) = sin(x)

# Chamada da função
erro = erro_serial(u, x, ne, m, h, npg, C, EQoLG)
println("Erro L2: ", erro)
```

---

### `C0_options`

Calcula o vetor `C0` de coeficientes iniciais, dado o tipo de inicialização especificado.

#### Argumentos:

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

#### Retorno:

Retorna o vetor `C0` com os coeficientes iniciais, de acordo com a opção selecionada.

#### Exemplo de Uso:

```julia
# Função de condição inicial
u0(x) = sin(x)

# Chamada da função com opção de interpolação
C0 = C0_options(1, u0, 0.0, ne, m, h, npg, EQoLG)
println(C0)
```

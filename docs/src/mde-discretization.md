# MDEDiscretization.jl - Discretização de Equações Diferenciais Modernas

O módulo **MDEDiscretization.jl** fornece ferramentas para a discretização de equações diferenciais ordinais (EDO) no contexto de elementos finitos. O objetivo principal deste módulo é facilitar a discretização do problema espacial e temporal, bem como fornecer funções auxiliares para resolver e analisar essas EDOs.

## Funções Exportadas

Este módulo exporta as seguintes variáveis e funções:

### Variáveis de Parâmetros da EDO

- **beta**: Parâmetro beta da equação diferencial.
- **alpha**: Parâmetro alpha da equação diferencial.
- **gamma**: Parâmetro gamma da equação diferencial.

### Variáveis de Intervalos

- **a**: Limite inferior do intervalo espacial.
- **b**: Limite superior do intervalo espacial.
- **T**: Limite superior do intervalo temporal.
- **npg**: Número de pontos de Gauss.
- **ne**: Número de intervalos no espaço.
- **m**: Número de pontos internos no espaço.
- **h**: Passo no espaço.
- **tau**: Passo no tempo.
- **N**: Número de intervalos no tempo.

### Funções de Forma

- **φ1(ξ)**: Função de forma φ1 no elemento finito.
- **φ2(ξ)**: Função de forma φ2 no elemento finito.
- **dφ1(ξ)**: Derivada da função de forma φ1.
- **dφ2(ξ)**: Derivada da função de forma φ2.

### Funções de Avaliação nos Pontos de Gauss

- **φ1P**: Avaliação de φ1 nos pontos de Gauss.
- **φ2P**: Avaliação de φ2 nos pontos de Gauss.
- **dφ1P**: Derivada de φ1 nos pontos de Gauss.
- **dφ2P**: Derivada de φ2 nos pontos de Gauss.
- **Wφ1P**: Avaliação de φ1 multiplicada pelos pesos de Gauss.
- **Wφ2P**: Avaliação de φ2 multiplicada pelos pesos de Gauss.

### Funções Auxiliares

- **u(x, t)**: Função solução da EDO no ponto \( (x, t) \).
- **u0(x)**: Função inicial \( u_0(x) \).
- **du0(x)**: Derivada da função inicial \( u_0(x) \).
- **g(s)**: Função auxiliar \( g(s) \).
- **f(x, t)**: Função \( f(x, t) \), dependente de \( x \) e \( t \).
- **monta_LG(ne)**: Função que monta a matriz de conectividade do elemento finito.
- **monta_EQ(ne)**: Função que monta o vetor de numeração das funções de forma.

## Exemplos de Uso

### 1. Solução da EDO

Para avaliar a solução \( u(x, t) \) da equação diferencial, você pode usar a função `u(x, t)`:

```julia
x = 0.5   # Ponto espacial
t = 0.2   # Ponto temporal
solucao = u(x, t)
println("Solução u($x, $t) = $solucao")
```

### 2. Função Inicial \( u_0(x) \)

Para avaliar a função inicial \( u_0(x) \), você pode usar a função `u0(x)`:

```julia
x = 0.5   # Ponto espacial
u0_value = u0(x)
println("Valor de u0($x) = $u0_value")
```

### 3. Derivada da Função Inicial \( u_0(x) \)

Para calcular a derivada da função inicial \( u_0(x) \), use a função `du0(x)`:

```julia
x = 0.5   # Ponto espacial
du0_value = du0(x)
println("Derivada de u0($x) = $du0_value")
```

### 4. Função \( g(s) \)

Para calcular o valor da função auxiliar \( g(s) \), use a função `g(s)`:

```julia
s = 0.3   # Ponto no domínio de g
g_value = g(s)
println("Valor de g($s) = $g_value")
```

### 5. Função \( f(x, t) \)

Para calcular a função \( f(x, t) \), use a função `f(x, t)`:

```julia
x = 0.5   # Ponto espacial
t = 0.2   # Ponto temporal
f_value = f(x, t)
println("Valor de f($x, $t) = $f_value")
```

### 6. Montagem da Matriz de Conectividade

Para montar a matriz de conectividade do elemento finito, use a função `monta_LG(ne)`:

```julia
ne = 2^6   # Número de elementos finitos
conectividade = monta_LG(ne)
println("Matriz de Conectividade: \n$conectividade")
```

### 7. Montagem do Vetor de Numeração das Funções de Forma

Para montar o vetor de numeração das funções de forma, use a função `monta_EQ(ne)`:

```julia
ne = 2^6   # Número de elementos finitos
numeracao = monta_EQ(ne)
println("Vetor de Numeração: \n$numeracao")
```

## Descrição das Funções

### `u(x::Float64, t::Float64) -> Float64`

Avalia a solução \( u(x, t) \) da equação diferencial no ponto \( (x, t) \).

#### Argumentos

- `x::Float64`: Ponto no domínio espacial de \( u \).
- `t::Float64`: Ponto no domínio temporal de \( u \).

#### Retorno

- Valor de \( u(x, t) \) no ponto \( (x, t) \).

### `u0(x::Float64) -> Float64`

Avalia a função inicial \( u_0(x) \) no ponto \( x \).

#### Argumentos

- `x::Float64`: Ponto no domínio espacial de \( u_0 \).

#### Retorno

- Valor de \( u_0(x) \) no ponto \( x \).

### `du0(x::Float64) -> Float64`

Calcula a derivada de \( u_0(x) \) no ponto \( x \).

#### Argumentos

- `x::Float64`: Ponto no domínio espacial de \( u_0 \).

#### Retorno

- Derivada de \( u_0(x) \) no ponto \( x \).

### `g(s::Float64) -> Float64`

Avalia a função auxiliar \( g(s) \) no ponto \( s \).

#### Argumentos

- `s::Float64`: Ponto no domínio de \( g \).

#### Retorno

- Valor de \( g(s) \) no ponto \( s \).

### `f(x::Float64, t::Float64) -> Float64`

Avalia a função \( f(x, t) \) no ponto \( (x, t) \).

#### Argumentos

- `x::Float64`: Ponto no domínio espacial de \( f \).
- `t::Float64`: Ponto no domínio temporal de \( f \).

#### Retorno

- Valor de \( f(x, t) \) no ponto \( (x, t) \).

### `monta_LG(ne::Int64) -> Matrix{Int64}`

Monta a matriz de conectividade do elemento finito.

#### Argumentos

- `ne::Int64`: Número de elementos finitos.

#### Retorno

- Matriz de conectividade.

### `monta_EQ(ne::Int64) -> Vector{Int64}`

Monta o vetor de numeração das funções de forma.

#### Argumentos

- `ne::Int64`: Número de elementos finitos.

#### Retorno

- Vetor com a numeração das funções de forma.

---

Este módulo pode ser utilizado para discretizar e resolver equações diferenciais em problemas de elementos finitos. Ele oferece funções úteis para avaliação de soluções e montagem das matrizes de conectividade, sendo uma ferramenta valiosa para problemas de simulação numérica e análise de EDOs.

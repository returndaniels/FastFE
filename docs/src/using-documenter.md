# Como utilizar **Documenter.jl** em seu projeto

O pacote **Documenter.jl** é uma ferramenta poderosa para gerar documentação de pacotes na linguagem Julia. Ele permite que você escreva documentação em markdown dentro do seu código e, em seguida, gere documentação estática ou interativa. Vou te guiar nos principais passos para usar o **Documenter.jl** para documentar seu pacote Julia.

### Passos principais:

#### 1. **Adicionar o Documenter.jl ao seu pacote**

Primeiro, você precisa adicionar o pacote `Documenter.jl` ao seu projeto:

```julia
using Pkg
Pkg.add("Documenter")
```

Ou, se você estiver dentro do ambiente do seu pacote, você pode adicionar diretamente:

```julia
using Pkg
Pkg.add("Documenter")
```

#### 2. **Criar o diretório de documentação**

Na estrutura do seu pacote, crie o diretório de documentação (`docs`) na raiz do pacote. A estrutura do diretório ficará assim:

```
MyPackage/
├── src/
│   └── MyPackage.jl
├── docs/
│   ├── make.jl
│   └── ...
└── Project.toml
```

No diretório `docs/`, você precisará de um arquivo `make.jl` que o **Documenter.jl** usará para gerar a documentação.

#### 3. **Escrever a documentação em Markdown**

No diretório `docs/`, crie um arquivo `index.md` que conterá a documentação principal do seu pacote. Este arquivo usará markdown para organizar o conteúdo.

````markdown
# MyPackage.jl Documentation

Bem-vindo à documentação do pacote **MyPackage.jl**.

## Funções

### `erro_vectorized`

```julia
erro_vectorized(u::Function, X::Matrix{Float64}, u_eval::Matrix{Float64}, ne::Int64, m::Int64, h::Float64, npg::Int64, C::Vector{Float64}, EQoLG::Matrix{Int64})::Float64
```
````

Calcula o erro para uma aproximação fornecida usando uma abordagem vetorizada.

...

````

#### 4. **Configurar o `make.jl`**

Crie o arquivo `make.jl` no diretório `docs/` para configurar a geração da documentação.

```julia
# docs/make.jl
using Documenter

# Inclui o pacote que você deseja documentar
using MyPackage

# Configurações do Documenter
makedocs(
    sitename = "MyPackage.jl Documentation",  # Nome do site/documentação
    modules = [MyPackage],                    # Pacote a ser documentado
    pages = [
        "Home" => "index.md",                 # Página principal da documentação
        "API" => "api.md"                     # Pode adicionar outras páginas, se necessário
    ],
    format = Documenter.HTML()                # Você pode escolher o formato HTML, Markdown, ou LaTeX
)
````

#### 5. **Documentação das funções e tipos**

Dentro do código do seu pacote, use os **docstrings** para documentar as funções e tipos. Por exemplo:

```julia
"""
    erro_vectorized(u::Function, X::Matrix{Float64}, u_eval::Matrix{Float64},
                     ne::Int64, m::Int64, h::Float64, npg::Int64, C::Vector{Float64}, EQoLG::Matrix{Int64})::Float64

Calcula o erro para uma aproximação fornecida usando uma abordagem vetorizada.

# Parâmetros
- `u`: Função que representa a solução exata.
- `X`: Matriz de pontos de discretização.
- `u_eval`: Matriz para armazenar os valores calculados da solução.
- `ne`: Número de elementos.
- `m`: Número de pontos de malha.
- `h`: Tamanho do passo.
- `npg`: Número de pontos por elemento.
- `C`: Vetor de coeficientes.
- `EQoLG`: Matriz de equações locais.

# Retorno
Retorna o erro calculado como um número de ponto flutuante.
"""
function erro_vectorized(u::Function, X::Matrix{Float64}, u_eval::Matrix{Float64},
                         ne::Int64, m::Int64, h::Float64, npg::Int64, C::Vector{Float64}, EQoLG::Matrix{Int64})::Float64
    ...
end
```

#### 6. **Gerar a documentação**

Para gerar a documentação, abra o REPL da Julia e execute:

```julia
using Documenter
include("docs/make.jl")
```

Isso vai gerar a documentação na pasta `docs/build/`. Você pode agora abrir a documentação gerada no seu navegador.

#### 7. **Servir a documentação localmente (opcional)**

Se você quiser servir a documentação localmente enquanto faz alterações, use:

```julia
using LiveServer
LiveServer.serve(dir="docs/build")
```

Isso vai abrir um servidor local e você poderá visualizar a documentação no seu navegador.

#### 8. **Publicar a documentação (opcional)**

Você pode publicar a documentação em sites como GitHub Pages. Para isso, você precisa configurar o repositório e garantir que a documentação seja gerada e publicada automaticamente ao fazer push para o repositório. O **Documenter.jl** possui suporte para configurar isso automaticamente com integração ao GitHub.

```julia
makedocs(
    sitename = "MyPackage.jl Documentation",
    modules = [MyPackage],
    pages = ["Home" => "index.md"],
    format = Documenter.HTML(),
    repository = "https://github.com/username/MyPackage.jl",
    deploy = Documenter.GitHubActions()  # Usando GitHub Actions para deploy automático
)
```

Essa configuração irá criar e publicar a documentação automaticamente no GitHub Pages.

### Resumo

Com o **Documenter.jl**, você pode gerar documentação do seu pacote Julia de forma simples e eficiente. O processo envolve:

1. Adicionar o pacote `Documenter.jl` ao seu projeto.
2. Criar um diretório `docs/` contendo arquivos markdown e o script `make.jl`.
3. Escrever a documentação no formato markdown.
4. Rodar o script `make.jl` para gerar a documentação HTML.
5. Publicar a documentação, se desejado.

Esses passos cobrem a criação e publicação básica da documentação do seu pacote Julia.

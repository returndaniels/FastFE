using Documenter

include("../src/FastFE.jl")

# Inclui o pacote que você deseja documentar
using ..FastFE

# Configurações do Documenter
makedocs(
    repo = "github.com/returndaniels/FastFE", # Repositório do GitHub
    sitename = "Documentação FastFE.jl",  # Nome do site/documentação
    # modules = [FastFE],                    # Pacote a ser documentado
    pages = [
        "Introdução" => "index.md",                 # Página principal da documentação
        "SerialFiniteElements.jl" => "serial-finite-elements.md",
        "VectorizedFiniteElements.jl" => "vectorized-finite-elements.md",
        "SerialErrorAnalysis.jl" => "serial-error-analysis.md",
        "VectorizedErrorAnalysis.jl" => "vectorized-error-analysis.md",
        "MDEDiscretization.jl" => "mde-discretization.md",
        "Benchmark.jl" => "benchmark.md",
        "Pacotes em Julia" => "generate-packages.md",
        "Documentação" => "using-documenter.md",
        "Deploy da Documentação" => "deploy.md"
    ],
    format = Documenter.HTML()                # Você pode escolher o formato HTML, Markdown, ou LaTeX
)

# if "deploy" in ARGS
#     include("./env.jl")
#   end
  
#   deploydocs(deps = nothing, make = nothing,
#     repo = "github.com/USER_NAME/PROJECT_NAME.jl.git",
#     target = "build",
#     branch = "main",
#     devbranch = "main",
#   )
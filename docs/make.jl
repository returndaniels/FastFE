using Documenter

# Inclui o pacote que você deseja documentar
using ..FastFE

# Configurações do Documenter
makedocs(
    repo = "github.com/returndaniels/FastFE", # Repositório do GitHub
    sitename = "Documentação FastFE.jl",  # Nome do site/documentação
    modules = [FastFE],                    # Pacote a ser documentado
    pages = [
        "FastFE.jl" => "index.md",                 # Página principal da documentação
        "SerialFiniteElements.jl" => "serial-finite-elements.md",
        "VectorizedFiniteElements.jl" => "vectorized-finite-elements.md",
        "SerialErrorAnalysis.jl" => "serial-error-analysis.md",
        "VectorizedErrorAnalysis.jl" => "vectorized-error-analysis.md",
        "MDEDiscretization.jl" => "mde-discretization.md",
        "Benchmark.jl" => "benchmark.md",
        "Criar e Registarar Pacotes" => "julia-package.md",
        "Configurar pacotes" => "generate-packages.md",
        "Usando o Documenter.jl" => "using-documenter.md"
    ],
    format = Documenter.HTML()                # Você pode escolher o formato HTML, Markdown, ou LaTeX
)
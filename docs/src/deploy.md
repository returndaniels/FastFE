# Deploy de Documentação com Documenter.jl

Este guia apresenta como configurar um pipeline de CI/CD para gerar e publicar automaticamente a documentação do seu pacote Julia usando o [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl). Com isso, você pode manter a documentação atualizada e acessível para os usuários, economizando tempo e reduzindo erros manuais. Você pode verificar informações mais detalhadas em [Documenter.jl: Hosting Documentation](https://documenter.juliadocs.org/stable/man/hosting/).

---

## Por que usar Documenter.jl no CI/CD?

- **Automação**: Garante que a documentação seja gerada sempre que houver alterações relevantes.
- **Acessibilidade**: Facilita a publicação em plataformas como GitHub Pages.
- **Consistência**: Mantém a documentação sincronizada com o código do projeto.

---

## Configurando o Workflow no GitHub Actions

Para gerar e publicar a documentação automaticamente, você pode usar o GitHub Actions. O arquivo YAML abaixo descreve as etapas para configurar o processo:

### Arquivo `.github/workflows/deploy_docs.yml`

```yaml
name: Build and Deploy Documentation

on:
  push:
    branches:
      - main

jobs:
  deploy-docs:
    runs-on: ubuntu-latest

    steps:
      - name: Checkout Code
        uses: actions/checkout@v3

      - name: Set up Julia
        uses: julia-actions/setup-julia@v1
        with:
          version: "1.9"

      - name: Install Dependencies
        run: julia -e 'using Pkg; Pkg.instantiate()'

      - name: Build Documentation
        run: julia --project=docs docs/make.jl

      - name: Deploy Documentation to GitHub Pages
        uses: peaceiris/actions-gh-pages@v3
        with:
          github_token: ${{ secrets.GITHUB_TOKEN }}
          publish_dir: docs/build
```

---

### Explicação das Etapas

1. **`on.push.branches`**:

   - O workflow é acionado sempre que houver um _push_ na branch `main`.

2. **Etapas do Job `deploy-docs`**:
   - **_Checkout_ do código**: Baixa os arquivos do repositório para o ambiente.
   - **Configuração do Julia**: Define a versão do Julia a ser usada no workflow.
   - **Instalação de dependências**: Garante que as dependências sejam baixadas e configuradas com `Pkg.instantiate()`.
   - **Geração da documentação**: Executa o script `docs/make.jl` para criar a documentação.
   - **Publicação no GitHub Pages**: Publica os arquivos gerados no branch `gh-pages`.

---

## Script de Documentação

Crie ou ajuste o arquivo `docs/make.jl` no diretório `docs/`. Um exemplo básico de script:

```julia
using Documenter
using MeuPacote

makedocs(
    sitename = "MeuPacote",
    modules = [MeuPacote],
    format = Documenter.HTML(),
    clean = true
)
```

Esse script realiza as seguintes tarefas:

- Define o nome do site com o parâmetro `sitename`.
- Especifica os módulos para os quais a documentação será gerada.
- Configura a saída no formato HTML.

---

## Publicação no GitHub Pages

1. **Habilite o GitHub Pages no repositório**:

   - Vá em `Settings > Pages`.
   - Configure o branch `gh-pages` como fonte de publicação.

2. **Resultado do Workflow**:
   - A cada execução bem-sucedida, a documentação será publicada automaticamente no endereço:
     `https://<seu-usuario>.github.io/<seu-repositorio>`.

---

## Testando e Validando

1. Faça alterações no repositório e realize um _push_:

   ```bash
   git add .
   git commit -m "Adiciona CI/CD para documentação"
   git push origin main
   ```

2. Acompanhe a execução do workflow na aba `Actions` do repositório no GitHub.
3. Verifique se a documentação foi publicada corretamente no GitHub Pages.

---

## Conclusão

Integrar o Documenter.jl com CI/CD oferece uma maneira eficaz de manter a documentação do seu pacote Julia sempre atualizada. Automatizando esse processo, você:

- Reduz erros manuais na geração e publicação.
- Garante que a documentação esteja sempre alinhada com o código.
- Oferece aos usuários acesso rápido e fácil às informações atualizadas.

Siga as etapas deste guia e otimize o fluxo de trabalho do seu projeto Julia! 🚀

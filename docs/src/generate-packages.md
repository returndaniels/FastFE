# Criar Pacotes em _Julia_

Para transformar seu projeto em Julia em um pacote completo e funcional, siga estas etapas que respeitam as convenções padrão da linguagem. Este guia aborda desde a estruturação inicial até a distribuição e registro do pacote.

---

## Parte 1: Configurando o Pacote

### 1. Estrutura do Projeto

Certifique-se de que o projeto tenha uma estrutura de diretório compatível com pacotes Julia. O layout típico é:

```
MeuPacote/
├── Project.toml
├── src/
│   └── MeuPacote.jl
└── test/
    └── runtests.jl
```

- **`MeuPacote`**: Substitua pelo nome do seu pacote. Use o estilo _CamelCase_.
- **`src/MeuPacote.jl`**: Arquivo principal que define o módulo do pacote.
- **`test/runtests.jl`**: Opcional, mas altamente recomendado para incluir testes automatizados.

### 2. Criar o Arquivo `Project.toml`

O arquivo `Project.toml` contém metadados sobre o pacote, como nome, versão e dependências. Para criá-lo automaticamente:

```julia
using Pkg
Pkg.generate("MeuPacote")
```

Isso cria a estrutura básica do pacote com um arquivo `Project.toml`. Personalize-o para incluir informações relevantes:

```toml
name = "MeuPacote"
uuid = "123e4567-e89b-12d3-a456-426614174000"
version = "0.1.0"

[deps]  # Dependências externas
Example = "7876af07-990d-54b4-ab0e-23690620f79a"

[compat]  # Compatibilidade com versões do Julia
julia = "1.6"
```

### 3. Definir o Módulo Principal

No arquivo `src/MeuPacote.jl`, defina o módulo principal:

```julia
module MeuPacote

export minha_funcao

minha_funcao(x) = x^2

end
```

### 4. Adicionar Dependências

Se o pacote depende de outros pacotes Julia, adicione-os ao ambiente:

```julia
using Pkg
Pkg.add("Example")
```

Isso atualiza automaticamente a seção `[deps]` no `Project.toml`.

### 5. Adicionar Testes

Crie o diretório `test/` e adicione o arquivo `runtests.jl`:

```julia
using Test
using MeuPacote

@testset "Testando minha_funcao" begin
    @test minha_funcao(2) == 4
    @test minha_funcao(0) == 0
end
```

Execute os testes com:

```julia
Pkg.test("MeuPacote")
```

---

## Parte 2: Registrando e Distribuindo o Pacote

### 6. Versionar com Git

Inicie o controle de versão do seu pacote:

```bash
git init
git add .
git commit -m "Primeira versão do MeuPacote"
```

Publique o repositório online, por exemplo, no GitHub ou GitLab.

### 7. Registrar o Pacote

Para compartilhar o pacote com a comunidade Julia:

1. Hospede o código em um repositório público.
2. Adicione o bot [Registrator.jl](https://github.com/JuliaRegistries/Registrator.jl) como colaborador ao repositório.
3. Solicite o registro comentando no GitHub:

   ```
   @JuliaRegistrator register
   ```

4. Certifique-se de que o `Project.toml` está correto e que o pacote segue as regras do [General Registry](https://github.com/JuliaRegistries/General).

### 8. Publicar o Pacote

Depois de registrado, o pacote estará disponível para instalação via `Pkg.add`. Certifique-se de documentar bem o uso e manter atualizações consistentes.

---

## Parte 3: Boas Práticas e Dicas Extras

### Documentação

Considere usar ferramentas como [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) para criar documentação automatizada do pacote. Isso ajuda outros desenvolvedores a entender e usar seu código.

### Testes Abrangentes

Inclua testes para cobrir diferentes cenários de uso. Isso garante maior confiança na funcionalidade e estabilidade do pacote.

### Compatibilidade

Declare explicitamente a compatibilidade com versões do Julia e pacotes no `Project.toml`. Isso evita problemas para outros usuários.

---

Seguindo essas etapas, você transformará seu projeto Julia em um pacote robusto e pronto para distribuição! 🎉

# Como criar seu pacote em _Julia_

Para transformar seu projeto em Julia em um pacote, você deve seguir uma série de etapas que seguem as convenções padrão da linguagem. Isso facilita a instalação, o uso e o compartilhamento do pacote. Aqui está um guia para ajudá-lo:

---

### 1. **Estrutura do Projeto**

Certifique-se de que o projeto tenha uma estrutura de diretório compatível com pacotes Julia. O layout típico é:

```
MyPackage/
├── Project.toml
├── src/
│   └── MyPackage.jl
└── test/
    └── runtests.jl
```

- **`MyPackage`**: Substitua pelo nome do seu pacote. É uma boa prática usar o estilo _CamelCase_.
- **`src/MyPackage.jl`**: Arquivo principal onde você define o módulo do pacote.
- **`test/runtests.jl`**: Opcional, usado para escrever testes automatizados.

---

### 2. **Criar o Arquivo `Project.toml`**

O arquivo `Project.toml` contém metadados do pacote, como nome, versão e dependências. Para criá-lo automaticamente:

1. **Inicialize o ambiente do pacote**:

   ```julia
   using Pkg
   Pkg.generate("MyPackage")
   ```

   Isso criará a estrutura básica do pacote.

2. **Personalize o `Project.toml`**:
   Edite o arquivo para incluir informações relevantes. Exemplo:

   ```toml
   name = "MyPackage"
   uuid = "01234567-89ab-cdef-0123-456789abcdef"  # Gerado automaticamente
   version = "0.1.0"

   [deps]  # Dependências externas
   Example = "7876af07-990d-54b4-ab0e-23690620f79a"

   [compat]  # Compatibilidade com versões do Julia e pacotes
   julia = "1.6"
   ```

---

### 3. **Definir o Módulo Principal**

No arquivo `src/MyPackage.jl`, defina o módulo principal:

```julia
module MyPackage

export say_hello

say_hello() = println("Hello from MyPackage!")

end
```

---

### 4. **Adicionar Dependências**

Se o pacote depende de outros pacotes Julia, adicione-os ao ambiente:

```julia
using Pkg
Pkg.add("Example")
```

Isso atualiza automaticamente a seção `[deps]` no `Project.toml`.

---

### 5. **Adicionar Testes**

Crie testes no arquivo `test/runtests.jl`:

```julia
using Test
using MyPackage

@testset "Basic Tests" begin
    @test say_hello() == nothing
end
```

Execute os testes com:

```julia
Pkg.test("MyPackage")
```

---

### 6. **Registrar e Distribuir o Pacote**

1. Para uso local, você pode ativar e usar o pacote diretamente:

   ```julia
   Pkg.activate(".")
   Pkg.instantiate()
   ```

2. Para compartilhar, registre o pacote no **Julia General Registry**:
   - Hospede o código em um repositório público (por exemplo, GitHub).
   - Siga as instruções do [Julia Package Registry Guide](https://pkgdocs.julialang.org/v1/creating-packages/#Registering-Packages).

---

### 7. **Documentação (Opcional)**

Considere usar ferramentas como [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) para criar documentação.

---

### Exemplo Rápido

1. Crie o pacote:
   ```bash
   julia -e 'using Pkg; Pkg.generate("MyPackage")'
   ```
2. Edite `src/MyPackage.jl` para implementar seu código.
3. Adicione dependências e testes.
4. Ative o ambiente local para desenvolvimento:
   ```julia
   Pkg.activate(".")
   Pkg.instantiate()
   ```
5. Registre o pacote para distribuição pública, se necessário.

Após seguir esses passos, seu projeto estará configurado como um pacote Julia! 🎉

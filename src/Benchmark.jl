using .FiniteElements

# Este módulo contém as funções de benchmark para testar a performance de operações
# serial e vetorizadas em um modelo numérico.

### Constantes e Variáveis ###
# Parâmetros do modelo numérico
a = 0.0
b = 10.0
h = 0.1
tau = 0.01
T = 1.0
alpha = 1.0
beta = 0.5
gamma = 0.5
ne = 1000
npg = 10
m = ne - 1

# Função de benchmark para a comparação de funções serial e vetorizadas.
function run_benchmark()
    # Monta o vetor F com a função monta_EQ e monta_LG
    EQoLG = monta_EQ(ne)[monta_LG(ne)]

    # Vetores para armazenar os resultados
    F_ext_serial = zeros(Float64, ne)
    F_ext_vectorized = zeros(Float64, ne)

    # Avaliações de f para a vetorização
    X = ((h / 2) * (P .+ 1) .+ a)' .+ range(a, step=h, stop=b - h)
    f_eval = Matrix{Float64}(undef, ne, npg)
    values = Matrix{Float64}(undef, ne, 2)

    # Teste de desempenho das funções serial e vetorizadas para o vetor F
    println("Benchmarking F_serial! vs F_vectorized!")
    println(norm(F_serial!(F_ext_serial, x, x -> f(x, 0.5), ne, m, h, npg, EQoLG) -
                 F_vectorized!(F_ext_vectorized, X, f_eval, values, x -> f(x, 0.5), ne, m, h, npg, EQoLG)))

    # Benchmark
    @benchmark F_serial!(F_ext_serial, x, x -> f(x, 0.5), ne, m, h, npg, EQoLG)
    @benchmark F_vectorized!(F_ext_vectorized, X, f_eval, values, x -> f(x, 0.5), ne, m, h, npg, EQoLG)

    # Teste de desempenho para o vetor G
    println("Benchmarking G_serial! vs G_vectorized!")
    C = C0_options(1, u0, a, ne, m, h, alpha, beta, gamma, npg, EQoLG)
    G_ext_serial = zeros(Float64, ne)
    G_ext_vectorized = zeros(Float64, ne)
    g_eval = Matrix{Float64}(undef, ne, npg)

    # Benchmark de G_serial! e G_vectorized!
    println(norm(G_serial!(G_ext_serial, C, ne, m, h, npg, EQoLG) -
                 G_vectorized!(G_ext_vectorized, g_eval, values, C, ne, m, h, npg, EQoLG)))
    
    @benchmark G_serial!(G_ext_serial, C, ne, m, h, npg, EQoLG)
    @benchmark G_vectorized!(G_ext_vectorized, g_eval, values, C, ne, m, h, npg, EQoLG)

    # Teste de desempenho para a matriz K
    println("Benchmarking K_serial vs K_vectorized!")
    println(norm(K_serial(ne, m, h, npg, alpha, beta, gamma, EQoLG) -
                 K_vectorized(ne, m, h, npg, alpha, beta, gamma, EQoLG)))

    @benchmark K_serial(ne, m, h, npg, alpha, beta, gamma, EQoLG)
    @benchmark K_vectorized(ne, m, h, npg, alpha, beta, gamma, EQoLG)
end

### Funções Auxiliares ###

# Função para calcular o erro entre a solução exata e a solução aproximada
function erro_QG(u, a, ne, m, h, npg, C, EQoLG)
    # Implemente a lógica para calcular o erro, pode envolver
    # avaliação da solução exata e a aproximação dada por C
end

# Função de inicialização de parâmetros
function monta_EQ(ne)
    # Implementar a função que gera a estrutura do sistema de equações
end

function monta_LG(ne)
    # Implementar a função que gera a estrutura do sistema local-global
end

# Função principal para inicializar o modelo e executar o benchmark
function main()
    run_benchmark()

    ### Estudo de erro ###
    tam = 6
    erros = zeros(tam - 1)
    hs = zeros(tam - 1)

    option = 1  # Opção de como calcular C0

    # Loop para estudar o erro em diferentes discretizações
    for i in 2:tam
        ne_er = 2^i
        N_er = ne_er
        m_er = ne_er - 1
        h_er = (b - a) / ne_er
        tau_er = h_er
        hs[i - 1] = h_er
        EQoLG = monta_EQ(ne_er)[monta_LG(ne_er)]

        # Inicialização do erro
        erro = zeros(N_er + 1)
        C0 = C0_options(option, u0, a, ne_er, m_er, h_er, alpha, beta, gamma, npg, EQoLG)
        C1_tiu = LU_dec\(F_serial!(F_ext_serial, x, x -> f(x, tau_er * 0.5), ne_er, m_er, h_er, npg, EQoLG) + (M / tau_er - K / 2) * C0 -
                        G_serial!(G_ext_serial, C0, ne_er, m_er, h_er, npg, EQoLG))
        C1 = LU_dec\(F_serial!(F_ext_serial, x, x -> f(x, tau_er * 0.5), ne_er, m_er, h_er, npg, EQoLG) + (M / tau_er - K / 2) * C0 -
                     G_serial!(G_ext_serial, (C0 + C1_tiu) / 2, ne_er, m_er, h_er, npg, EQoLG))

        # Cálculo do erro
        erro[1] = erro_QG(u0, a, ne_er, m_er, h_er, npg, C0, EQoLG)
        erro[2] = erro_QG(x -> u(x, tau_er), a, ne_er, m_er, h_er, npg, C1, EQoLG)
        Cn_menos_2 = C0
        Cn_menos_1 = C1

        for n in 2:N_er
            temp = LU_dec\(F_serial!(F_ext_serial, x, x -> f(x, tau_er * (n - 0.5)), ne_er, m_er, h_er, npg, EQoLG) + (M / tau_er - K / 2) * Cn_menos_1 -
                G_serial!(G_ext_serial, (3 * Cn_menos_1 - Cn_menos_2) / 2, ne_er, m_er, h_er, npg, EQoLG))
            Cn_menos_2 = Cn_menos_1
            Cn_menos_1 = temp
            erro[n + 1] = erro_QG(x -> u(x, tau_er * n), a, ne_er, m_er, h_er, npg, Cn_menos_1, EQoLG)
        end
        erros[i - 1] = maximum(erro)
    end

    plot(hs, hs.^2, label="h^2", xlabel="h", title="Ordem de Convergência do Erro (Escala Log)", legend=:topright)
    plot!(hs, erros, yscale=:log10, xscale=:log10, label="Erro")
end


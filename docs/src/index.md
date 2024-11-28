### 1. Introdução e Formulação do Problema

Dados $\alpha > 0$, $\beta, \gamma \geq 0$ constantes reais e funções $f: [0, 1] \times [0, T] \rightarrow \mathbb{R}$, $u_0: [0, 1] \rightarrow \mathbb{R}$ e $g: \mathbb{R} \to \mathbb{R}$, determine $u: [0, 1] \times [0, T] \to \mathbb{R}$ tal que:

$$
\label{edo}
\begin{cases}
u_t(x, t) -\alpha u_{xx}(x, t) + \gamma u_{x}(x, t) + \beta u(x, t) + g(u(x, t)) = f(x, t), & \forall(x, t) \in (0, 1) \times [0, T] \\
u(0, t) = u(1, t) = 0, & \forall t \in [0, T] \\
u(x, 0) = u_0(x), & \forall x \in (0, 1)
\end{cases}. \tag{1}
$$

Sejam considerados os conjuntos:

$$
H \triangleq \{u \text{ suficientemente suave}: u(0) = u(1) = 0\} \\
V = H.
$$

$H$ é denominado espaço das soluções do problema e $V$ é denominado espaço das funções testes. O primeiro passo que iremos fazer é multiplicar ambos os lados da Equação Diferencial Ordinária (\ref{edo}) por uma função $v(x) \in V$:

$$
v(x)(u_t(x, t) - \alpha u_{xx}(x, t) + \gamma u_x(x, t) + \beta u(x, t) + g(u(x, t))) = v(x)f(x, t) \Rightarrow (\text{Distribui } v(x) ) \\
v(x)u_t(x, t) - \alpha v(x)u_{xx}(x, t) + \gamma v(x)u_x(x, t) + \beta v(x)u(x, t) + v(x)g(u(x, t)) = v(x)f(x, t).
$$

Agora integramos ambos os lados no intervalo $x \in [0, 1]$ e aproveitamos para reduzir a ordem da derivada de $u$:

$$
\int^1_0(v(x)u_t(x, t) - \alpha v(x)u_{xx}(x, t) + \gamma v(x)u_x(x, t) + \beta v(x)u(x, t) + v(x)g(u(x, t)))\,dx = \int^1_0v(x)f(x, t)\,dx \Rightarrow (\text{Linearidade da integral}) \\
\int^1_0v(x)u_t(x, t)\,dx - \alpha\int^1_0v(x)u_{xx}(x, t)\,dx + \gamma\int^1_0v(x)u_x(x, t)\,dx + \beta\int^1_0v(x)u(x, t)\,dx + \int^1_0v(x)g(u(x, t))\,dx = \int^1_0v(x)f(x, t)\,dx \Rightarrow (\text{Integração por partes}) \\
\int^1_0v(x)\frac{\partial u(x, t)}{\partial t}\,dx - \alpha\left(v(x)\frac{\partial u(x, t)}{\partial x}\Bigg|^1_0 -\int^1_0\frac{dv(x)}{\,dx}\frac{\partial u(x, t)}{\partial x}\,dx\right) + \gamma\int^1_0v(x)\frac{\partial u(x, t)}{\partial x}\,dx + \beta\int^1_0v(x)u(x, t)\,dx + \int^1_0v(x)g(u(x, t))\,dx = \int^1_0v(x)f(x, t)\,dx \Rightarrow (v(1) = v(0) = 0) \\
\int^1_0v(x)\frac{\partial u(x, t)}{\partial t}\,dx + \alpha\int^1_0\frac{dv(x)}{\,dx}\frac{\partial u(x, t)}{\partial x}\,dx + \gamma\int^1_0v(x)\frac{\partial u(x, t)}{\partial x}\,dx + \beta\int^1_0v(x)u(x, t)\,dx + \int^1_0v(x)g(u(x, t))\,dx = \int^1_0v(x)f(x, t)\,dx.
$$

Nesse formato, o problema possui apenas derivadas de primeira ordem, relaxando a necessidade da função $u$ ser duas vezes diferenciável em todo o seu domínio. Esse formato é denominado **Formulação Fraca** e será devidamente definido a seguir.

### 2. Formulação Fraca do Problema:

Dados $\alpha > 0$, $\beta, \gamma \geq 0$ constantes reais e funções funções $f: [0, 1] \times [0, T] \rightarrow \mathbb{R}$ e $u_0: [0, 1] \rightarrow \mathbb{R}$, determine $u: [0, 1] \times [0, T] \to \mathbb{R}$ tal que:

$$
\label{ffp}
\int^1_0v(x)\frac{\partial u(x, t)}{\partial t}\,dx + \alpha\int^1_0\frac{dv(x)}{\,dx}\frac{\partial u(x, t)}{\partial x}\,dx + \gamma\int^1_0v(x)\frac{\partial u(x, t)}{\partial x}\,dx + \beta\int^1_0v(x)u(x, t)\,dx + \int^1_0v(x)g(u(x, t))\,dx = \int^1_0v(x)f(x, t)\,dx, \quad \forall v \in V. \tag{2}
$$

Será utilizada a seguinte notação:

$$
\begin{cases}
\kappa(\cdot, \cdot): V \times V \to \mathbb{R} \\
(v, u) \mapsto \alpha\int^1_0\frac{dv(x)}{\,dx}\frac{du(x)}{dx}\,dx + \gamma\int^1_0v(x)\frac{du(x)}{dx}\,dx + \beta\int^1_0v(x)u(x)\,dx \\
\end{cases} \quad
\begin{cases}
(\cdot, \cdot): V \times C_{[0,\ 1]} \\
(v, u) \mapsto \int^1_0v(x)u(x)\,dx.
\end{cases}
$$

Em que $C_D \triangleq \{f: \mathbb{D} \to \mathbb{R}| f \text{ é contínua em todo domínio} D\}$. Nosso problema é então descrito por:

$$
\label{notacao}
(v, u_t(t)) + \kappa(v, u(t)) + (v, g(u(t))) = (v, f(t)), \quad \forall v \in V. \tag{3}
$$

Note que no problema definido em (\ref{notacao}), as funções transientes estão fixadas em um determinado tempo $t$. O próximo passo é aproximar o espaço $V$ por um espaço $V_m \triangleq \{\sum^m_{j=1}C_j\varphi_j(x): C_j \in \mathbb{R}, j = 1,\ 2,\ \ldots,\ m\}$, em que as funções $\varphi_i$, $\forall i \in \{1,\ 2,\ \ldots,\ m\}$, formam uma base para um subespaço finito das funções teste $V$. Vamos considerar a seguinte função construída como combinação linear das funções da base do subespaço finito:

$$
\tilde{u}(x) = \sum^m_{j=1}C_j\varphi_j(x).
$$

Essa será a função que utilizaremos para resolver numericamente o problema no espaço. Agora podemos passar para um problema aproximado que será definido logo abaixo.

### 3. Problema Semi-discreto - Via Método de Galerkin na Variável Espacial:

Dados $\alpha > 0$, $\beta, \gamma \geq 0$ constantes reais e funções funções $f: [0, 1] \times [0, T] \rightarrow \mathbb{R}$ e $\tilde{u}_0 \in V_m$ uma aproximação para $u_0$, determine $\tilde{u}(t) \in V_m, t \in [0,\ T]$ tal que:

$$
\label{galerkin}
\begin{cases}
(\tilde{v}, \tilde{u}_t(t)) + \kappa(\tilde{v}, \tilde{u}(t)) + (\tilde{v}, g(\tilde{u}(t))) = (\tilde{v}, f(t)), \quad \forall \tilde{v} \in V_m \\
\tilde{u}(0) = \tilde{u}_0
\end{cases}. \tag{4}
$$

Queremos determinar $\tilde{u}(t)$, isto é, queremos determinar os coeficientes $C_i$, $\forall i \in \{1,\ 2,\ \ldots,\ m\}$ tal que $\tilde{u}(t)$ satisfaça (\ref{galerkin}). Porém, a nossa definição de $\tilde{u}(x)$ não depende do tempo. A estratégia é definir um conjunto de funções $\tilde{u}^{(n)}(x)$, em que teremos que determinar as componentes do vetor $C^{(n)}$ para cada tempo $t_n$. Primeiro, vamos fazer uma aproximação de ordem quadrática para a derivada $\tilde{u}_t(t_{n-\frac{1}{2}}) \approx \frac{\tilde{u}(t_{n}) - \tilde{u}(t_{n-1})}{\tau}$ (Método Crank-Nicolson), em que $t_{n-\frac{1}{2}} = \frac{t_n + t_{n-1}}{2}$, e discretizar o domínio do tempo uniformemente a partir de um passo $\tau \in \mathbb{R^+}$ de modo que $0 = t_0 < t_1 < \ldots < t_N \leq T, t_n = n\tau\ \forall n \in \{0,\ \ldots,\ N\}$. Nosso problema totalmente discreto é definido a seguir.

### 4. Problema Totalmente Discreto - Galerkin no Espaço e Crank-Nicolson no Tempo:

Dados $\alpha > 0$, $\beta, \gamma \geq 0$ constantes reais e funções funções $f: [0, 1] \times [0, T] \rightarrow \mathbb{R}$ e $\tilde{u}_0 \in V_m$ uma aproximação para $u_0$, determine $U^{(n)} \in V_m$ tal que:

$$
\begin{cases}\label{crank_nic}
\left(\tilde{v}, \frac{U^{(n)} - U^{(n-1)}}{\tau}\right) + \kappa\left(\tilde{v}, \frac{U^{(n)} + U^{(n-1)}}{2}\right) + (\tilde{v}, g(U^{(n-\frac{1}{2})})) = \left(\tilde{v}, f(t_{n - \frac{1}{2}})\right), \quad \forall \tilde{v} \in V_m \\
U^{(0)} = \tilde{u}_0
\end{cases}. \tag{5}
$$

Para $U^{(n)} := \tilde{u}(t_n)$, $\tilde{u}(t_{n - \frac{1}{2}}) \approx \frac{U^{(n)} + U^{(n-1)}}{2}$, $0 = t_0 < t_1 < \ldots < t_N \leq T, t_n = n\tau\ \forall n \in \{0,\ \ldots,\ N\}$ e $t_{n-\frac{1}{2}} = \frac{t_n + t_{n-1}}{2}$. Tendo em vista que $\tilde{v} \in V_m$, em particular, $\tilde{v}$ pode ser uma das funções da base. Essa escolha irá ser de grande ajuda para determinar os coeficientes de $U^{(n)}(x) = \sum^m_{j=1}C^{(n)}_j\varphi_j(x)$, como veremos a seguir:

$$
\left(\tilde{v}, \frac{U^{(n)} - U^{(n-1)}}{\tau}\right) + \kappa\left(\tilde{v}, \frac{U^{(n)} + U^{(n-1)}}{2}\right) + \left(\tilde{v}, g(U^{(n-\frac{1}{2})})\right) = \left(\tilde{v}, f(t_{n - \frac{1}{2}})\right) \Rightarrow (\text{Definição de } U) \\
\left(\tilde{v}, \frac{\sum^m_{j=1}C^{(n)}_j\varphi_j - \sum^m_{j=1}C^{(n-1)}_j\varphi_j}{\tau}\right) + \kappa\left(\tilde{v}, \frac{\sum^m_{j=1}C^{(n)}_j\varphi_j + \sum^m_{j=1}C^{(n-1)}_j\varphi_j}{2}\right) + \left(\tilde{v}, g(U^{(n-\frac{1}{2})})\right) = \left(\tilde{v}, f(t_{n - \frac{1}{2}})\right) \Rightarrow (\text{Escolha de função da base}) \\
\left(\varphi_i, \frac{\sum^m_{j=1}C^{(n)}_j\varphi_j - \sum^m_{j=1}C^{(n-1)}_j\varphi_j}{\tau}\right) + \kappa\left(\varphi_i, \frac{\sum^m_{j=1}C^{(n)}_j\varphi_j + \sum^m_{j=1}C^{(n-1)}_j\varphi_j}{2}\right) + \left(\varphi_i, g(U^{(n-\frac{1}{2})})\right) = \left(\varphi_i, f(t_{n - \frac{1}{2}})\right) \Rightarrow (\text{Rearranjando o somatório}) \\
\left(\varphi_i, \frac{\sum^m_{j=1}(C^{(n)}_j - C^{(n-1)}_j)\varphi_j}{\tau}\right) + \kappa\left(\varphi_i, \frac{\sum^m_{j=1}(C^{(n)}_j + C^{(n-1)}_j)\varphi_j}{2}\right) + \left(\varphi_i, g(U^{(n-\frac{1}{2})})\right) = \left(\varphi_i, f(t_{n - \frac{1}{2}})\right) \Rightarrow (\text{Operador } \kappa \text{ bilinear}) \\
\sum^m_{j=1}\left(\varphi_i, \varphi_j\right)\frac{(C^{(n)}_j - C^{(n-1)}_j)}{\tau} + \sum^m_{j=1}\kappa\left(\varphi_i, \varphi_j\right)\frac{(C^{(n)}_j + C^{(n-1)}_j)}{2} + \left(\varphi_i, g(U^{(n-\frac{1}{2})})\right) = \left(\varphi_i, f(t_{n - \frac{1}{2}})\right), \quad \forall i \in \{1,\ 2,\ \ldots,\ m\}.
$$

### 5. Lidando com a Função Não-linear g(s) (Linearizando com Ordem Quadrática):

Vamos utilizar a aproximação $U^{(n-\frac{1}{2})} \approx \frac{3U^{(n-1)} - U^{(n-2)}}{2}$ para obter a seguinte equação:

$$
\sum^m_{j=1}\left(\varphi_i, \varphi_j\right)\frac{(C^{(n)}_j - C^{(n-1)}_j)}{\tau} + \sum^m_{j=1}\kappa\left(\varphi_i, \varphi_j\right)\frac{(C^{(n)}_j + C^{(n-1)}_j)}{2} + \left(\varphi_i, g(\frac{3U^{(n-1)} - U^{(n-2)}}{2})\right) = \left(\varphi_i, f(t_{n - \frac{1}{2}})\right) \Rightarrow \\
\sum^m_{j=1}\left(\varphi_i, \varphi_j\right)\frac{(C^{(n)}_j - C^{(n-1)}_j)}{\tau} + \sum^m_{j=1}\kappa\left(\varphi_i, \varphi_j\right)\frac{(C^{(n)}_j + C^{(n-1)}_j)}{2} = \left(\varphi_i, f(t_{n - \frac{1}{2}})\right) - \left(\varphi_i, g(\frac{3U^{(n-1)} - U^{(n-2)}}{2})\right), \quad \forall i \in \{1,\ 2,\ \ldots,\ m\}.
$$

Podemos montar um sistema linear $m\times m$ para um vetor $C^{(n)}_{1\times m}$ de incógnitas descrito a seguir.

### 6. Formulação Matricial:

Vamos considerar a última equação acima para $i = 1,\ 2,\ \ldots,\ m$:

$$
\begin{cases}
\sum^m_{j=1}\left(\varphi_1, \varphi_j\right)\frac{(C^{(n)}_j - C^{(n-1)}_j)}{\tau} + \sum^m_{j=1}\kappa\left(\varphi_1, \varphi_j\right)\frac{(C^{(n)}_j + C^{(n-1)}_j)}{2} = \left(\varphi_1, f(t_{n - \frac{1}{2}})\right) - \left(\varphi_1, g(\frac{3U^{(n-1)} - U^{(n-2)}}{2})\right) \\
\sum^m_{j=1}\left(\varphi_2, \varphi_j\right)\frac{(C^{(n)}_j - C^{(n-1)}_j)}{\tau} + \sum^m_{j=1}\kappa\left(\varphi_2, \varphi_j\right)\frac{(C^{(n)}_j + C^{(n-1)}_j)}{2} = \left(\varphi_2, f(t_{n - \frac{1}{2}})\right) - \left(\varphi_2, g(\frac{3U^{(n-1)} - U^{(n-2)}}{2})\right) \\
\vdots \\
\sum^m_{j=1}\left(\varphi_m, \varphi_j\right)\frac{(C^{(n)}_j - C^{(n-1)}_j)}{\tau} + \sum^m_{j=1}\kappa\left(\varphi_m, \varphi_j\right)\frac{(C^{(n)}_j + C^{(n-1)}_j)}{2} = \left(\varphi_m, f(t_{n - \frac{1}{2}})\right) - \left(\varphi_m, g(\frac{3U^{(n-1)} - U^{(n-2)}}{2})\right)
\end{cases}.
$$

Sendo assim, matricialmente, temos que:

$$
\frac{1}{\tau}
\begin{bmatrix}
(\varphi_1, \varphi_1) & (\varphi_1, \varphi_2) & \cdots & (\varphi_1, \varphi_m) \\
(\varphi_2, \varphi_1) & (\varphi_2, \varphi_2) & \cdots & (\varphi_2, \varphi_m) \\
\vdots & \cdots & \ddots & \vdots \\
(\varphi_m, \varphi_1) & (\varphi_m, \varphi_2) & \cdots & (\varphi_m, \varphi_m)
\end{bmatrix}
\begin{bmatrix}
C^{(n)}_1 - C^{(n-1)}_1 \\
C^{(n)}_2 - C^{(n-1)}_2 \\
\vdots \\
C^{(n)}_m - C^{(n-1)}_m
\end{bmatrix}
+
\frac{1}{2}
\begin{bmatrix}
\kappa(\varphi_1, \varphi_1) & \kappa(\varphi_1, \varphi_2) & \cdots & \kappa(\varphi_1, \varphi_m) \\
\kappa(\varphi_2, \varphi_1) & \kappa(\varphi_2, \varphi_2) & \cdots & \kappa(\varphi_2, \varphi_m) \\
\vdots & \cdots & \ddots & \vdots \\
\kappa(\varphi_m, \varphi_1) & \kappa(\varphi_m, \varphi_2) & \cdots & \kappa(\varphi_m, \varphi_m)
\end{bmatrix}
\begin{bmatrix}
C^{(n)}_1 + C^{(n-1)}_1 \\
C^{(n)}_2 + C^{(n-1)}_2 \\
\vdots \\
C^{(n)}_m + C^{(n-1)}_m
\end{bmatrix}
= \\
\begin{bmatrix}
(\varphi_1, f(t_{n - \frac{1}{2}})) \\
(\varphi_2, f(t_{n - \frac{1}{2}})) \\
\vdots \\
(\varphi_m, f(t_{n - \frac{1}{2}}))
\end{bmatrix}
-
\begin{bmatrix}
\left(\varphi_1, g(\frac{3U^{(n-1)} - U^{(n-2)}}{2})\right) \\
\left(\varphi_2, g(\frac{3U^{(n-1)} - U^{(n-2)}}{2})\right) \\
\vdots \\
\left(\varphi_m, g(\frac{3U^{(n-1)} - U^{(n-2)}}{2})\right)
\end{bmatrix}
$$

Reescrevendo a equação acima, temos o seguinte sistema linear:

$$
\label{matriz_vetor}
\frac{1}{\tau}M(C^{(n)} - C^{(n-1)}) + \frac{1}{2}K(C^{(n)} + C^{(n-1)}) = F(n - \frac{1}{2}) - G(\frac{3U^{(n-1)} - U^{(n-1)}}{2}). \tag{6}
$$

Para $M_{ij} = (\varphi_i, \varphi_j)$, $K_{ij} = \kappa(\varphi_i, \varphi_j)$, $F_i(n - \frac{1}{2}) = (\varphi_i, f(t_{n - \frac{1}{2}}))$ e $G_i(\frac{3U^{(n-1)} - U^{(n-1)}}{2}) = \left(\varphi_i, g(\frac{3U^{(n-1)} - U^{(n-2)}}{2})\right)$.

### 7. Implementação:

Estamos utilizando uma função linear para $\varphi_i(x)$ definida por:

$$
\varphi_i(x) =
\begin{cases}
\frac{x - x_{i-1}}{x_{i} - x_{i-1}}, x_{i-1} \leq x \leq x_{i}  \\
\frac{x_{i+1} - x}{x_{i+1} - x_{i}}, x_{i} \leq x \leq x_{i+1}  \\
0, \text{ caso contrário}
\end{cases} \quad \text{e} \quad
\frac{d\varphi_i(x)}{dx} =
\begin{cases}
\frac{1}{x_{i} - x_{i-1}}, x_{i-1} \leq x \leq x_{i} \\
\frac{-1}{x_{i+1} - x_{i}}, x_{i} \leq x \leq x_{i+1} \\
0, \text{caso contrário}
\end{cases}, \text{ para } i = 1,\ 2,\ \ldots,\ m.
$$

Em que $x_i = a + i(x_{i} - x_{i-1}) = a + ih_i$ e $h_i = x_{i} - x_{i-1}$, para $i = 1,\ 2,\ \ldots,\ m$. Como as partições são uniformes, $h_i = h$ é constante. **Uma consideração importante é que estaremos utilizando uma extensão dos vetores e das matrizes para $m+1$ ao invés de $m$ para manter a escrita das equações mais limpa e facilitar o uso das estruturas locais e globais**. As matrizes $M$ e $K$, e os vetores de termos independentes $F(n-\frac{1}{2})$ e $G(U)$ são calculados a partir da aproximação das respectivas integrais por quadratura gaussiana (notação local):

$$
\bullet K^e_{ab} = \alpha\int^{x^e_2}_{x^e_1}\frac{d\varphi^e_a(x)}{\,dx}\frac{d\varphi^e_b(x)}{dx}\,dx + \gamma\int^{x^e_2}_{x^e_1}\varphi^e_a(x)\frac{d\varphi^e_b(x)}{dx}\,dx + \beta\int^{x^e_2}_{x^e_1}\varphi^e_a(x)\varphi^e_b(x)\,dx = \\
\frac{2\alpha}{h}\int^{1}_{-1}\frac{d\phi_a(\xi)}{d\xi}\frac{d\phi_b(\xi)}{d\xi}\,d\xi + \gamma\int^{1}_{-1}\phi_a(\xi)\frac{d\phi_b(\xi)}{d\xi}\,d\xi + \frac{\beta h}{2}\int^{1}_{-1}\phi_a(\xi)\varphi_b(\xi)\,d\xi \approx \\
\frac{2\alpha}{h}\sum^{N_{PG}}_{j=1}W_j\frac{d\phi_a(X_j)}{d\xi}\frac{d\phi_b(X_j)}{d\xi} + \gamma\sum^{N_{PG}}_{j=1}W_j\phi_a(X_j)\frac{d\phi_b(X_j)}{d\xi} + \frac{\beta h}{2}\sum^{N_{PG}}_{j=1}W_j\phi_a(X_j)\varphi_b(X_j), \quad \text{para } e \in \{1,\ 2,\ \ldots,\ m+1\}, a, b \in \{1, 2\}.
$$

$$
\bullet M^e_{ab} = \int^{x^e_2}_{x^e_1}\varphi^e_a(x)\varphi^e_b(x)\,dx = \frac{h}{2}\int^{1}_{-1}\phi_a(\xi)\varphi_b(\xi)\,d\xi \approx \frac{h}{2}\sum^{N_{PG}}_{j=1}W_j\phi_a(X_j)\varphi_b(X_j), \quad \text{para } e \in \{1,\ 2,\ \ldots,\ m+1\}, a, b \in \{1, 2\}.
$$

Em que $\phi_1(\xi) = \frac{1 - \xi}{2}$ e $\phi_2(\xi) = \frac{1 + \xi}{2}$, $X_j$ e $W_j$, para $j = 1,\ 2,\ \ldots,\ N_{PG}$ são os pontos e pesos de Gauss, respectivamente, $N_{Pg}$ é o número de pontos de Gauss, e $x^e_{1} = x_{e-1}$ e $x^e_{2} = x_{e}$ (pontos começando de $0$ até $m+1$, com $m$ pontos internos).

$$
\bullet F^e_{a}(n - \frac{1}{2}) = \int^{x^e_2}_{x^e_1}\varphi^e_a(x)f(x, t_{n - \frac{1}{2}})\,dx = \frac{h}{2}\int^1_{-1}\phi_a(\xi)f(\frac{h}{2}(1 + \xi) + x_{e-1}, t_{n - \frac{1}{2}})d\xi \approx \\
\frac{h}{2}\sum^{N_{PG}}_{j=1}W_j\phi_a(X_j)f(\frac{h}{2}(X_j+1) + x_{e-1}, t_{n - \frac{1}{2}}), \quad \text{para } e \in \{1,\ 2,\ \ldots,\ m+1\}, a \in \{1, 2\}.
$$

$$
\bullet G^e_a(U) = \int^{x^e_2}_{x^e_1}\varphi_a(x)g(\sum^m_{j=1}C_j\varphi_j(x))\,dx = \int^{x^e_2}_{x^e_1}\varphi_a(x)g(C_{e-1}\varphi_{e-1}(x) + C_{e}\varphi_e(x))\,dx = \int^{1}_{-1}\phi_a(\xi)g(C_{EQ[LG[1, e]]}\phi_1(\xi) + C_{EQ[LG[2, e]]}\phi_2(\xi))\frac{h}{2}\,d\xi \approx \\
\frac{h}{2}\sum^{N_{PG}}_{j=1}W_j\phi_a(X_j)g(C_{EQ[LG[1, e]]}\phi_1(X_j) + C_{EQ[LG[2, e]]}\phi_2(X_j)), \quad \text{para } e \in \{1,\ 2,\ \ldots,\ m+1\}, a \in \{1, 2\}.
$$

Vamos considerar as seguintes estruturas para mapear os elementos locais para os elementos globais:

$$
LG =
\begin{bmatrix}
1 & 2 & \cdots & m+1 \\
2 & 3 & \cdots & m + 2
\end{bmatrix}_{2\times(m+1)} \quad \text{ e } \quad
EQ =
\begin{bmatrix}
m+1 & 1 & 2 & \cdots & m & m+1
\end{bmatrix}_{1\times(m+2)}.
$$

Portanto:
\begin{align}
M*{EQ[LG[a, e]], EQ[LG[b, e]]} &\gets M^e*{ab} \\
K*{EQ[LG[a, e]], EQ[LG[b, e]]} &\gets K^e*{ab} \\
F*{EQ[LG[a, e]]} &\gets F^e_a \\
G*{EQ[LG[a, e]]} &\gets G^e_a
\end{align}

O erro entre a solução aproximada e a solução exata é dado por:

$$
||u(x) - \tilde{u}(x)||^2_{L^2} = \int^1_0|u(x) - \tilde{u}(x)|^2dx = \sum^{m+1}_{i=1}\int^{x_i}_{x_{i-1}}\left(u(x) - \sum^m_{j=1}C_j\varphi_j(x)\right)^2dx = \\
\int^{x_1}_{x_0}(u(x) - C_1\varphi_1(x))^2dx + \sum^{m}_{i=2}\int^{x_i}_{x_{i-1}}\left(u(x) - C_{i-1}\varphi_{i-1}(x) - C_{i}\varphi_{i}(x)\right)^2dx + \int^{x_{m+1}}_{x_m}(u(x) - C_m\varphi_1(m))^2dx = \\
\frac{h}{2}\left[\int^{1}_{-1}(u(\frac{h}{2}(\xi + 1) + x_0) - C_1\frac{1 + \xi}{2})^2d\xi + \sum^{m}_{i=2}\int^{1}_{-1}\left(u(\frac{h}{2}(\xi + 1) + x_{i-1}) - C_{i-1}\frac{1 - \xi}{2} - C_{i}\frac{1 + \xi}{2}\right)^2d\xi + \int^{1}_{-1}(u(\frac{h}{2}(\xi + 1) + x_m) - C_m\frac{1 - \xi}{2})^2d\xi\right] \approx \\
\frac{h}{2}\left[\sum^{N_{PG}}_{j=1}W_j\left(u(\frac{h}{2}(X_j + 1) + x_{0}) - C_1\frac{1 + X_j}{2}\right)^2 + \sum^m_{i=2}\sum^{N_{PG}}_{j=1}W_j\left(u(\frac{h}{2}(X_j + 1) + x_{i-1}) - C_{i-1}\frac{1 - X_j}{2} - C_{i}\frac{1 + X_j}{2}\right)^2 + \sum^{N_{PG}}_{j=1}W_j\left(u(\frac{h}{2}(X_j + 1) + x_{m}) - C_m\frac{1 - X_j}{2}\right)^2\right].
$$

### 8. Aproximação de $u_0$ - Escolhas de $U^{(0)} \in V_m$:

#### 8.1 $U^0$ como interpolante de $u_0$:

A primeira escolha é dada pela interpolante de $u_0$. Temos que:

$$
U^0(x_i) = \sum^m_{j=1}C^0_j\varphi_j(x_i) = C^0_i.
$$

Sendo assim:

$$
C^0 =
\begin{bmatrix}
u_0(x_1) \\
u_0(x_2) \\
\vdots \\
u_0(x_m)
\end{bmatrix}.
$$

#### 8.2 $U^0$ como projeção $L^2$ de $u_0$:

Seja $U^0 \in V_m$ tal que:

$$
(U^0 - u_0, v_h) = 0 \Rightarrow (U^0, v_h) - (u_0, v_h) = 0, \quad \forall v_h \in V_m.
$$

Tomando $U^0(x) = \sum^m_{j=1}C^0_j\varphi_j(x)$ e $v_h = \varphi_i$, para $i = 1,\ 2,\ \ldots,\ m$, na equação acima, temos que:

$$
MC^0 =
\begin{bmatrix}
(u_0, \varphi_1) \\
(u_0, \varphi_2) \\
\vdots \\
(u_0, \varphi_m)
\end{bmatrix}.
$$

#### 8.3 $U^0$ como projeção $H^1_0$ de $u_0$:

Seja $U^0 \in V_m$ tal que:

$$
\left(\frac{d}{dx}\left(U^0 - u_0\right), \frac{dv_h}{dx}\right) = 0 \Rightarrow \left(\frac{dU^0}{dx}, \frac{dv_h}{dx}\right) - \left(\frac{du_0}{dx}, \frac{dv_h}{dx}\right) = 0, \quad \forall v_h \in V_m.
$$

Tomando $U^0(x) = \sum^m_{j=1}C^0_j\varphi_j(x)$ e $v_h = \varphi_i$, para $i = 1,\ 2,\ \ldots,\ m$, na equação acima, temos que:

$$
\begin{bmatrix}
(\frac{d\varphi_1}{dx}, \frac{d\varphi_1}{dx}) & (\frac{d\varphi_1}{dx}, \frac{d\varphi_2}{dx}) & \cdots & (\frac{d\varphi_1}{dx}, \frac{d\varphi_m}{dx}) \\
(\frac{d\varphi_2}{dx}, \frac{d\varphi_1}{dx}) & (\frac{d\varphi_2}{dx}, \frac{d\varphi_2}{dx}) & \cdots & (\frac{d\varphi_2}{dx}, \frac{d\varphi_m}{dx}) \\
\vdots & \vdots & \ddots & \vdots \\
(\frac{d\varphi_m}{dx}, \frac{d\varphi_1}{dx}) & (\frac{d\varphi_m}{dx}, \frac{d\varphi_2}{dx}) & \cdots & (\frac{d\varphi_m}{dx}, \frac{d\varphi_m}{dx})
\end{bmatrix}C^0
=
\begin{bmatrix}
(\frac{du_0}{dx}, \frac{d\varphi_1}{dx}) \\
(\frac{du_0}{dx}, \frac{d\varphi_2}{dx}) \\
\vdots \\
(\frac{du_0}{dx}, \frac{d\varphi_m}{dx})
\end{bmatrix}.
$$

O lado direito pode ser obtido a partir das seguintes equações:

$$
\begin{cases}
(\frac{du_0}{dx}, \frac{d\varphi^{(e)}_1}{dx}) = \int^1_0\frac{du_0}{dx}\frac{d\varphi_1}{dx}\,dx = -\frac{1}{h}\int^1_0\frac{du_0}{dx}\,dx = (u_0(a + (e-1)h) - u_0(a + eh))/h \\
(\frac{du_0}{dx}, \frac{d\varphi^{(e)}_1}{dx}) = \int^1_0\frac{du_0}{dx}\frac{d\varphi_1}{dx}\,dx = \frac{1}{h}\int^1_0\frac{du_0}{dx}\,dx = (u_0(a + (e-1)h) - u_0(a + eh))/h
\end{cases}, \quad \text{para } e \in \{1,\ 2,\ \ldots,\ m+1\}, a \in \{1, 2\}. \\
V_{EQ[LG[a, e]]} \gets V^{(e)}_a,\quad \text{em que } V \text{ representa o vetor do lado direito}.
$$

#### 8.4 Utilizando o operador $\kappa(\cdot, \cdot)$ para realizar a projeção de $u_0$:

Seja $U^0 \in V_m$ tal que:

$$
\kappa\left(\left(U^0 - u_0\right), v_h\right) = 0 \Rightarrow \kappa\left(U^0, v_h\right) - \kappa\left(u_0, v_h\right) = 0, \quad \forall v_h \in V_m.
$$

Tomando $U^0(x) = \sum^m_{j=1}C^0_j\varphi_j(x)$ e $v_h = \varphi_i$, para $i = 1,\ 2,\ \ldots,\ m$, na equação acima, temos que:

$$
KC^0
=
\begin{bmatrix}
\kappa(u_0, \varphi_1) \\
\kappa(u_0, \varphi_2) \\
\vdots \\
\kappa(u_0, \varphi_m)
\end{bmatrix}.
$$

O lado direito pode ser obtido a partir das seguintes equações:

$$
\begin{cases}
\kappa(u_0, \varphi^{(e)}_1) = \alpha\int^1_0\frac{du_0(x)}{\,dx}\frac{d\varphi^{(e)}_1(x)}{dx}\,dx + \gamma\int^1_0u_0(x)\frac{d\varphi^{(e)}_1(x)}{dx}\,dx + \beta\int^1_0u_0(x)\varphi^{(e)}_1(x)\,dx = \frac{\alpha}{h}(u_0(a+(e-1)h)-u_0(a+eh)) + \frac{\gamma h}{2}\sum^{N_{pg}}_{j=1}W_j\left(\frac{du_0}{\,dx}\left(\frac{h}{2}(X_j + 1) + (e-1)h\right)\phi_1(X_j)\right) + \frac{\beta h}{2}\sum^{N_{pg}}_{j=1}W_j\left(u_0\left(\frac{h}{2}(X_j + 1) + (e-1)h\right)\frac{1 - X_j}{2}\right) \\
\kappa(u_0, \varphi^{(e)}_2) = \alpha\int^1_0\frac{du_0(x)}{\,dx}\frac{d\varphi^{(e)}_2(x)}{dx}\,dx + \gamma\int^1_0u_0(x)\frac{d\varphi^{(e)}_2(x)}{dx}\,dx + \beta\int^1_0u_0(x)\varphi^{(e)}_2(x)\,dx = \frac{\alpha}{h}(u_0(a+eh) - u_0(a+(e-1)h)) + \frac{\gamma h}{2}\sum^{N_{pg}}_{j=1}W_j\left(\frac{du_0}{\,dx}\left(\frac{h}{2}(X_j + 1) + (e-1)h\right)\phi_2(X_j)\right) + \frac{\beta h}{2}\sum^{N_{pg}}_{j=1}W_j\left(u_0\left(\frac{h}{2}(X_j + 1) + (e-1)h\right)\frac{1 + X_j}{2}\right)
\end{cases}, \quad \text{para } e \in \{1,\ 2,\ \ldots,\ m+1\}, a \in \{1, 2\}. \\
V_{EQ[LG[a, e]]} \gets V^{(e)}_a,\quad \text{em que } V \text{ representa o vetor do lado direito}.
$$

### 9. Determinando $U^{(1)}$ - Método Preditor-Corretor:

Precisamos determinar $U^{(1)}$ para poder utilizar a relação de recorrência da forma matricial do problema. Para isso, vamos utilizar um método preditor-corretor que em um passo nos dá uma aproximação no tempo $t_1$ a partir das duas seguintes etapas:

#### (i) Etapa 1:

Definimos $\tilde{U}^{(1)} \in V_m$ como solução do seguinte problema:

$$
\left(\tilde{v}, \frac{\tilde{U}^{(1)} - U^{(0)}}{\tau}\right) + \kappa\left(\tilde{v}, \frac{\tilde{U}^{(1)} + U^{(0)}}{2}\right) = \left(\tilde{v}, f(t_{n - \frac{1}{2}})\right) - (\tilde{v}, g(U^{(0)})), \quad \forall \tilde{v} \in V_m.
$$

#### (ii) Etapa 2:

Determinamos uma aproximação para $U^{(1)} \in V_m$ que é a solução do seguinte problema:

$$
\left(\tilde{v}, \frac{U^{(1)} - U^{(0)}}{\tau}\right) + \kappa\left(\tilde{v}, \frac{U^{(1)} + U^{(0)}}{2}\right) = \left(\tilde{v}, f(t_{n - \frac{1}{2}})\right) - (\tilde{v}, g(\frac{\tilde{U}^{(1)} + U^{(0)}}{2})), \quad \forall \tilde{v} \in V_m.
$$

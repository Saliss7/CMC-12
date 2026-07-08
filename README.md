# CMC-12 — Projeto de Exame

**Estudo Comparativo de Estimadores de Estado num Pêndulo Invertido sob Malha de Controle**

Instituto Tecnológico de Aeronáutica (ITA) — Sistemas de Controle Contínuos e Discretos (CMC-12)
Repositório: <https://github.com/Saliss7/CMC-12>

---

## 1. Objetivos do Projeto

Implementar e comparar, dentro de uma malha de controle, diferentes técnicas de
**estimação de estados** aplicadas a um **pêndulo invertido sobre carro (cart-pole)**
— um sistema dinâmico não-linear.

O curso de CMC-12 abordou apenas **filtragem passa-baixas clássica** e nenhuma
estimação de estados por espaço de estados. Este projeto avança sobre isso ao
implementar e confrontar quatro níveis de estimadores, do visto em aula ao estado
da arte:

```
Diferenciação + passa-baixas  →  KF linear  →  EKF  →  UKF
   (baseline do curso)          (linearizado)  (Jacobiano online)  (sigma points)
```

### Por que o pêndulo invertido
O controle opera em **duas fases** que exercitam os filtros em regimes distintos:

- **Swing-up** (energy-based): leva o pêndulo por **ângulos grandes** → regime
  fortemente não-linear, onde o KF linear degrada e EKF/UKF se destacam.
- **Estabilização** (LQR na vertical): regime quase-linear, onde todos os
  filtros convergem de forma semelhante.

Assim, um único experimento evidencia **onde cada filtro compensa**.

### Escopo da comparação
- **Estados:** posição do carro `x`, velocidade `ẋ`, ângulo `θ`, velocidade angular `θ̇`.
- **Medições:** apenas as posições `x` e `θ` com ruído (encoders). As velocidades
  `ẋ` e `θ̇` são estimadas.
- **Métricas:** RMSE por estado, tempo de convergência, robustez a ruído alto,
  a baixa taxa de amostragem / perda de medida, a erro de sintonia de `Q`/`R`,
  sensibilidade à inicialização e custo computacional.
- **Análise avançada:** consistência estatística (NEES/NIS — teste qui-quadrado).
- **Malha fechada:** desempenho do controle com estado estimado por cada filtro
  vs estado ideal (princípio da separação na prática).

### Cenários de teste
| ID | Cenário |
|----|---------|
| S1 | Estabilização na vertical (quase-linear) |
| S2 | Swing-up (não-linearidade forte) |
| S3 | Ruído de medição alto |
| S4 | Baixa taxa de amostragem / perda de medidas |
| S5 | Sensibilidade a erro de sintonia de `Q`/`R` |
| S6 | Inicialização ruim (covariância inicial grande) |

Cada cenário é avaliado em Monte Carlo (múltiplas sementes) para RMSE estatístico.

---

## 2. Divisão de Tarefas

> Preencher os nomes/handles do GitHub de cada integrante.

| | **Pessoa A — Planta, Controle & Harness** | **Pessoa B — Filtros clássicos + KF linear** | **Pessoa C — Filtros não-lineares (EKF/UKF)** |
|---|---|---|---|
| **Integrante** | _(a definir)_ | _(a definir)_ | _(a definir)_ |
| **Dono de** | Dinâmica não-linear do cart-pole; simulador da "verdade"; modelo de sensores (ruído, taxa, dropout, outliers); swing-up + LQR; harness de métricas/gráficos | Baseline (diferenciação + passa-baixas); filtro complementar; KF linear; linearização (Jacobiano na vertical) e discretização | EKF (Jacobianos online); UKF (sigma points / unscented transform); sintonia de `Q`/`R` dos não-lineares |
| **Entregáveis** | `plant_cartpole.m`, `sensor_model.m`, `controller_lqr.m`, `swingup.m`, `run_scenario.m`, `compute_metrics.m` | `est_lowpass.m`, `est_complementary.m`, `est_kf.m` + doc da linearização | `est_ekf.m`, `est_ukf.m` |
| **Seção do relatório** | Modelagem, controle e metodologia de simulação | KF linear, linearização e filtragem clássica | EKF/UKF e análise de consistência |

**Compartilhado (todos):** resultados comparativos, discussão, conclusão e manual do usuário.

### Contrato de interfaces (fixado no Dia 1)
Para permitir trabalho paralelo, todo o grupo fixa três interfaces antes de codar:

1. **Modelo de dados do experimento** — log padrão `[t, x_true(4), u, z(2)]`.
2. **Assinatura comum de estimador** (drop-in para qualquer filtro):
   ```matlab
   [xhat, P] = estimator_update(xhat_prev, P_prev, z_k, u_k, params)
   ```
3. **Harness de avaliação** — recebe qualquer estimador + um cenário e devolve
   métricas e gráficos.

Com o contrato fixo, cada pessoa testa seu módulo com dados sintéticos antes de a
planta final estar pronta.

### Fluxo em fases (paralelismo)
```
Dia 1  ── TODOS: fixar contrato + parâmetros físicos do pêndulo
Fase 1 ── A: planta+sensores+controle | B: KF linear+baseline | C: EKF+UKF     (PARALELO)
Marco 1 ─ A entrega datasets reais + harness → B e C plugam seus filtros
Fase 2 ── Integração: rodar os 4 filtros nos cenários S1..S6                    (semi-paralelo)
Fase 3 ── TODOS: análise + relatório + manual do usuário
```

---

## 3. Convenção de Commits e Branches

### Branches
- `main` — sempre estável; só recebe merge via Pull Request revisado.
- Uma branch por pessoa/módulo, prefixada pela área:
  - `plant/*` — Pessoa A (planta, controle, harness)
  - `filter-linear/*` — Pessoa B (KF linear, clássicos)
  - `filter-nonlinear/*` — Pessoa C (EKF, UKF)
  - `report/*` — texto do relatório
- Exemplo: `filter-nonlinear/ekf-jacobian`, `plant/swingup`.

### Mensagens de commit (Conventional Commits)
Formato: `<tipo>(<escopo>): <descrição no imperativo>`

- **Tipos:** `feat`, `fix`, `refactor`, `docs`, `test`, `chore`.
- **Escopos sugeridos:** `plant`, `sensor`, `control`, `kf`, `ekf`, `ukf`,
  `lowpass`, `harness`, `metrics`, `report`.

Exemplos:
```
feat(plant): adiciona dinâmica não-linear do cart-pole
feat(ekf): implementa predição com Jacobiano analítico
fix(kf): corrige sinal na matriz de discretização
docs(report): escreve fundamentação teórica do UKF
test(harness): adiciona cenário S4 (perda de medidas)
```

### Regras de convivência (para paralelismo sem conflito)
- **Não editar o contrato de interfaces sozinho.** Mudança nas assinaturas =
  aviso ao grupo + PR revisado por todos.
- Commits pequenos e frequentes; cada commit compila/roda.
- Um arquivo `est_*.m` tem um único dono → evita conflito de merge.
- `main` só avança por PR; ninguém dá `push` direto em `main`.
- Sincronizar (`git pull --rebase`) antes de abrir PR.

---

## 4. Estrutura de Pastas (proposta)
```
CMC-12/
├── README.md
├── src/
│   ├── plant_cartpole.m
│   ├── sensor_model.m
│   ├── controller_lqr.m
│   ├── swingup.m
│   ├── est_lowpass.m
│   ├── est_complementary.m
│   ├── est_kf.m
│   ├── est_ekf.m
│   ├── est_ukf.m
│   ├── run_scenario.m
│   └── compute_metrics.m
├── scenarios/        # definições de S1..S6
├── results/          # gráficos e tabelas gerados
└── report/           # relatório (LaTeX/IEEE) e manual do usuário
```

---

## 5. Como Executar
_(a preencher — manual do usuário: versão do MATLAB/Simulink, script de entrada
e como reproduzir cada cenário.)_

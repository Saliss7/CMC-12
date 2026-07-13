# CMC-12 — Projeto de Exame

**Estudo Comparativo de Estimadores de Estado num Pêndulo Invertido sob Malha de Controle**

Instituto Tecnológico de Aeronáutica (ITA) — Sistemas de Controle Contínuos e Discretos (CMC-12)
## Integrantes

- Pietro Maragno Trindade Marcucci
- Bernardo Affonso Cheron Vendramini
- Matheus Felipe Ramos Borges

## O que o projeto faz

Simula um pêndulo invertido sobre carro (cart-pole) em malha fechada (swing-up + LQR) e
compara, sob o mesmo controle, cinco estimadores de estado — passa-baixas, filtro
complementar, KF linear, EKF e UKF — em seis cenários de ruído/perda de medida/sintonia,
avaliando RMSE, consistência estatística (NEES/NIS) e custo computacional.

---

## Estrutura de pastas

```
CMC-12/
├── README.md
├── CHANGELOG.md      # registro da integração da Fase 2 e bugs corrigidos
├── src/                # implementação (planta, controle, estimadores)
│   ├── plant_cartpole.m     # dinâmica não-linear do cart-pole
│   ├── sensor_model.m       # modelo de sensores (ruído, taxa, dropout)
│   ├── controller_lqr.m     # controlador LQR (estabilização na vertical)
│   ├── swingup.m             # controlador de swing-up (energy shaping)
│   ├── linearise_upright.m  # linearização + discretização em theta=0
│   ├── jacobian_f.m          # Jacobiano analítico da dinâmica (usado pelo EKF)
│   ├── est_lowpass.m         # baseline: diferenciação + passa-baixas
│   ├── est_complementary.m  # filtro complementar
│   ├── est_kf.m               # Kalman Filter linear
│   ├── est_ekf.m              # Extended Kalman Filter
│   ├── est_ukf.m              # Unscented Kalman Filter
│   ├── params.m               # parâmetros físicos e de sintonia
│   ├── run_scenario.m         # roda um estimador num cenário, gera log
│   ├── compute_metrics.m     # RMSE, NEES/NIS, custo computacional
│   └── main.m                  # ponto de entrada: roda tudo (5 estimadores x 7 cenários)
├── scenarios/          # definição dos cenários de teste S0..S6
│   ├── S0_ideal.m              # referência ideal: sem ruído de medição
│   ├── S1_stabilisation.m    # estabilização na vertical (quase-linear)
│   ├── S2_swingup.m           # swing-up (não-linearidade forte)
│   ├── S3_high_noise.m        # ruído de medição alto
│   ├── S4_low_rate.m          # baixa taxa de amostragem / perda de medidas
│   ├── S5_mistuned_QR.m       # sensibilidade a erro de sintonia de Q/R
│   └── S6_bad_init.m          # inicialização ruim (covariância inicial grande)
├── tests/               # testes unitários e geração das tabelas/figuras de resultados
│   ├── test_classical_filters.m    # baseline, complementar, KF linear
│   ├── test_nonlinear_filters.m    # EKF, UKF
│   ├── build_report_results.m      # gera tabelas/figuras a partir dos .mat em results/
│   └── animate_cartpole.m           # animação MATLAB (verdade vs. estimativa)
└── results/             # gráficos, tabelas e animações gerados (saída, não editar à mão)
```

---

## Como executar

### Requisitos
- MATLAB (sem toolboxes especiais — apenas funções base; `chi2inv` requer a
  **Statistics and Machine Learning Toolbox**, usada em `build_report_results.m`).

### Rodar tudo (simulação + métricas)
Da raiz do repositório:
```matlab
run('src/main.m')
```
Isso executa os 5 estimadores (`lowpass`, `complementary`, `kf`, `ekf`, `ukf`) em cada um
dos 7 cenários (`S0`..`S6`, sendo `S0` a referência ideal sem ruído), calcula as métricas
(`compute_metrics.m`) e salva um `.mat` por combinação cenário/estimador em
`results/<cenario>_<estimador>.mat`.

### Gerar tabelas e figuras de resultados
Após rodar `main.m` (os `.mat` em `results/` precisam existir):
```matlab
addpath('src', 'tests');
build_report_results
```
Gera em `results/` a tabela-resumo (RMSE, NEES/NIS, tempo de convergência, custo
computacional por par cenário/estimador) e as figuras de resultados.

### Rodar os testes unitários
Da raiz do repositório:
```matlab
addpath('src', 'tests');
test_classical_filters   % baseline, complementar, KF linear
test_nonlinear_filters   % EKF, UKF
```
Ou via linha de comando (headless):
```bash
matlab -batch "addpath('src','tests'); test_classical_filters"
matlab -batch "addpath('src','tests'); test_nonlinear_filters"
```

### Reproduzir um cenário isolado
Cada cenário em `scenarios/` retorna uma struct consumida por `run_scenario.m`; para
rodar só um filtro num só cenário:
```matlab
addpath('src', 'scenarios');
p  = params();
sc = S2_swingup();
log     = run_scenario(@est_ukf, sc, p);
metrics = compute_metrics(log, p);
```

### Visualizar a animação do cart-pole
```matlab
addpath('src', 'tests');
animate_cartpole
```
Mostra a haste do pêndulo (verdade vs. estimativa) em movimento para inspecionar
qualitativamente o desempenho de um estimador.

### Visualizar os resultados (painel HTML)
Após rodar `main.m` e `build_report_results`, abra `results/dashboard.html` diretamente
no navegador (não requer servidor) para explorar interativamente as trajetórias, o
RMSE, o NEES/NIS e o custo computacional dos cinco estimadores em cada cenário S1–S6.

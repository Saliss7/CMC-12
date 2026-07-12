# Changelog — Integração da Fase 2 e Seção C no relatório

Registro da sessão em que o harness completo (Pessoa A) e os filtros
não-lineares (Pessoa C) foram rodados de ponta a ponta pela primeira vez, dois
bugs de integração encontrados nesse processo foram corrigidos, e `main.tex`
foi atualizado para refletir os resultados reais.

## Bugs de integração corrigidos

1. **`src/run_scenario.m` — sentinela de *dropout* inconsistente.**
   `sensor_model.m` sinaliza medição perdida com `[]` (vazio), mas todos os
   `est_*.m` (lowpass, complementar, KF, EKF, UKF) verificam
   `any(~isfinite(z_k))`, esperando o sentinela `[NaN; NaN]`. Como
   `any(~isfinite([]))` é `false`, o `[]` passava direto: os filtros baseados
   em covariância silenciosamente pulavam o tratamento de *dropout*, e
   `est_lowpass.m` quebrava com erro de indexação (`z_k(1)` sobre vetor
   vazio) assim que o cenário S4 (baixa taxa/perda de medidas) era executado.
   Corrigido convertendo `[]` para `nan(2,1)` em `run_scenario.m` antes de
   chamar o estimador.

2. **`scenarios/S3_high_noise.m` — ruído do sensor não refletido no `R` do
   estimador.** O cenário multiplicava `sigma_x`/`sigma_theta` por 10× (ruído
   real injetado pelo sensor), mas nunca definia `scenario.R` — e
   `run_scenario.m` só lê `scenario.R`, não os campos `sigma_*`. Resultado: os
   estimadores continuavam confiando na covariância de medição nominal
   (pequena) enquanto o ruído real era 10× maior, produzindo uma divergência
   espúria por *R* mal informado (NEES do KF ~87000) que se confundia com o
   teste de ruído alto propriamente dito — esse é o papel do cenário S5
   (`mistuned_QR`), não do S3. Corrigido definindo
   `s.R = diag([s.sigma_x^2, s.sigma_theta^2])` em `S3_high_noise.m`, para que
   o filtro seja "justamente" informado do ruído real.

## O que foi rodado

- `tests/test_classical_filters.m` e `tests/test_nonlinear_filters.m` — todos
  os testes passam (12/12).
- `src/main.m` — os 5 estimadores (lowpass, complementar, KF, EKF, UKF) nos 6
  cenários (S1–S6), malha fechada real (`swingup.m` → `controller_lqr.m`)
  alimentada pela estimativa, não pelo estado verdadeiro. Saída em
  `results/*.mat` (log + métricas por cenário/estimador).
- `tests/build_report_results.m` (novo script) — agrega os 30 `.mat` em
  `results/summary_table.csv` e gera 4 figuras novas (tema claro, prontas
  para impressão): trajetória de $\theta$ no *swing-up* com os 5 estimadores,
  NEES do *swing-up* (KF vs. EKF vs. UKF), RMSE de $\theta$ S1 vs. S3 (5
  estimadores) e custo computacional médio por estimador.

## O que mudou em `main.tex`

- **Nova Seção 8 — "Estimadores de Estado Não-Lineares: EKF e UKF"**
  (`sec:ekf-ukf`, responsabilidade da Pessoa C), antes inexistente. Cobre o
  Jacobiano analítico (`jacobian_f.m`), a predição do EKF (RK4 + linearização
  local via `expm`), a transformada *unscented* escalonada do UKF (*sigma
  points*, pesos de Van der Merwe), verificação numérica
  (Jacobiano vs. diferença finita, erro $2{,}51\times10^{-9}$), os resultados
  do *swing-up* (S2) mostrando a divergência otimista do KF linear
  (NEES $\sim\!2{,}5\times10^5$ contra o limite $\chi^2_4$ de $11{,}14$) e a
  robustez de EKF/UKF, e a síntese comparativa final com os 5 estimadores nos
  6 cenários.
- **Seção "Estimadores de Estado Clássicos" (Seção B):** corrigidos o
  *footnote* que dizia que a Seção C "ainda não existe", a nota de que o NIS
  estava "pendente de implementação" (já implementado), e a descrição de
  `linearise_upright.m` como tendo *placeholders* (já preenchida). Adicionado
  parágrafo "Sensibilidade a $Q/R$ mal ajustados (S5)" com números reais, e
  parágrafo "Da demonstração à Fase 2" encaminhando para os resultados reais
  na Seção 8.
- **Seção "Arquitetura da malha e discussão":** pequenas correções e nova
  **Conclusão** sintetizando os achados de todo o projeto e listando trabalhos
  futuros (recalibração de $Q$ por Van Loan, forma de Joseph uniforme no UKF,
  anomalia de $\dot\theta$ do UKF em S6).
- Bibliografia: adicionadas as referências de Julier & Uhlmann (1997) e
  Wan & van der Merwe (2000), base teórica do UKF.
- `main.pdf` recompilado (29 páginas), sem referências ou citações quebradas.

## Arquivos novos/gerados

- `tests/build_report_results.m` — script de agregação de resultados.
- `results/*.mat` (30 arquivos) — log e métricas de cada par
  (cenário, estimador).
- `results/summary_table.csv`, `results/S2_swingup_theta.png`,
  `results/S2_nees.png`, `results/S1_S3_rmse_theta.png`,
  `results/cpu_time.png`.

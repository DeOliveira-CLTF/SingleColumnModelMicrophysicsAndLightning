# mp_offline — driver offline (SCM) para microfísica MONAN/MPAS

Física extraída **sem alteração** (MPAS-Model, NCAR/MMM-physics): Thompson, WSM6, MYNN.
Drivers, pontes de dados e saídas: código nosso. Ver o relatório INPE (RPQ) para detalhes.

## Drivers
- `mp_offline_driver.F90` — caso GOAMAZON (DEPHY-SCM), Thompson ou WSM6, MYNN opcional.
- `ensemble_driver.F90` — ensemble ERA5 x perturbação de saturação (precisa dos dados ERA5;
  NÃO exercitado desde as correções abaixo, por falta dos dados).
- `src/unified_edmf/` — protótipo exploratório de EDMF unificado (Tan/Cohen/Lopez-Gomez);
  NÃO é parte do driver validado e não foi validado contra observações/LES.

## Chaves do namelist (mp_offline_driver)
`use_mynn`; `couple_mynn_condensate` (aplica rqcblten/rqiblten antes da microfísica);
`mynn_condensate_local_only` (aplica só a parte local, calculada por uma chamada-sombra do
bl_mynn_run com EDMF off a partir do mesmo estado); `mynn_edmf` (0 desliga o EDMF por
completo — CUIDADO: remove também a termodinâmica das plumas, confunde a comparação).
Controle (GOAMAZON, 35 h, ruído ~1%): acoplamento total -46% (Thompson) / -15% (WSM6) na chuva;
mistura local reduz (-1.20 / -0.90 mm), plumas EDMF devolvem (+0.41 / +0.43 mm), mesmo sinal
nos dois esquemas. Autoconversão medida: mistura local a reduz 91%/67% (espalha qc por mais células).
Ver relatório, Secs. 'Experimento de controle' e 'Teste da explicação'.

## Compilar

```
make              # compila os dois drivers (build/mp_offline_driver, build/ensemble_driver)
make -j4          # em paralelo
make DEBUG=1      # -O0 -g -fbacktrace em vez de -O2 (mais rapido pra recompilar, mais lento pra rodar)
make clean        # remove build/
```

Requer `gfortran` e uma instalacao de NetCDF-Fortran com `nf-config` no PATH.
Nao baixa nada: toda a fisica extraida (Thompson, WSM6, MYNN...) ja esta em `src/`.

Depois de compilar:
```
cd build
# copie ou edite um namelist.mp_offline (ver exemplos nas secoes acima e em ../*.nml)
./mp_offline_driver
```

ATENCAO (achado desta sessao): o Makefile foi adicionado tarde neste projeto; antes disso as
compilacoes eram manuais, arquivo a arquivo, e um `build/` reaproveitado por dezenas de patches
seguidos ficou com objetos desatualizados entre si (um .F90 mudava mas nem todo `.o` que dependia
dele era recompilado). Um build limpo com este Makefile é bit a bit igual ao build manual mais
recente nos casos conferidos, mas se algo parecer nao bater com o relatorio, o primeiro passo é
`make clean && make` antes de desconfiar da fisica.

## Duas opções de pluma dentro do MYNN (`mynn_edmf_type`)
- `0` (padrão): DMP_mf original do MPAS. Com `0` o resultado é IDÊNTICO bit a bit ao MYNN original.
- `1`: núcleo novo de plumas dentro de `src/mynn/bl_mynn_subroutines.F90` (sub-rotina interna
  `uedmf_core` do DMP_mf): fechamento de entranhamento/detranhamento de Cohen et al. (2020) com a
  TKE do MYNN, sem teto de 3 m/s, dM/dz=M(eps-delta). Condensação = a do MYNN; NÃO chama Thompson/WSM6
  (o condensado chega à microfísica da grade pelo caminho usual s_awqc/rqcblten).
  Diff vs. upstream: +569 linhas, 1 linha original alterada. bl_mynn.F90 e mynn_shared.F90 intactos.
  Tunáveis (NÃO calibrados): `mynn_uedmf_exc_mult`, `mynn_uedmf_kfac_max`; `mynn_uedmf_orig_gate=1`
  exige também os critérios de ativação do MPAS (controle).
  Termodinâmica/sumidouro (padrões = comportamento anterior, bit a bit):
  `mynn_uedmf_thermo`: 0 = condensation_edmf original; 1 = corrige o Exner ausente no THV latente;
  2 = 1 + calor latente de FUSAO (fracao de gelo como o xl_blend do MYNN). `mynn_uedmf_c0_auto` (1/m) +
  `mynn_uedmf_qc_crit`: sumidouro de condensado da pluma (a agua volta a GRADE no mesmo nivel; a pluma NAO
  produz chuva de superficie; MYNN nao tem canal de chuva).
- Profundidade medida pelo FLUXO DE MASSA (z onde M cai a 10%/1% do valor a 500 m), nao pelo topo da pluma mais alta.
  ATENCAO: o arquivo de forcante tem theta constante (367.7 K) de 16.5 a 20 km (camada neutra artificial): plumas
  com energia cinetica sobem ate o topo do dominio. GOAMAZON 35 h, Thompson: z10 mediano 3.2 km (original) ->
  4.4-8.5 km (nova, thermo=2, c0 0..1e-2); Exner e o fator dominante, fusao muda pouco, sumidouro aprofunda.
  Sem validacao observacional; c0 e o lancamento NAO calibrados. Hipotese TESTADA E REFUTADA: TKE~0 na alta troposfera
  NAO causa plumas profundas demais. Chaves: `mynn_uedmf_e0_floor` (piso uniforme de e0), `mynn_uedmf_econv=1`
  (TKE convectiva da fonte por detranhamento, Lopez-Gomez et al. 2020 Eq. 8, c_d=0.22) + `mynn_uedmf_econv_l`.
  Piso de e0 de 1e-4 a 3 m2/s2: z10 mediano 5.3->4.7-5.6 km, M>8km 8.4->9.8% (nao cai). Acima de 5 km:
  eps_dyn~1-4e-5, eps_turb~2-6e-6, detranhamento 5e-5..1e-4 (> entranhamento): pluma quase nao diluida.
  Termo dinamico (secao 'O termo dinamico'): eps_dyn ~ c_eps*b/w^2 (escala de Gregory), 5 ordens de grandeza;
  pluma nublada em ambiente seco: delta/eps ~ 3-4 (beta controla). Sensibilidade (um a um, Thompson): dominam
  c_eps (resposta ingreme: x2.5 -> x3 muda o regime de fluxo de massa alto para <=8 km), c_delta e beta; mu0, chi,
  c_lambda ~ irrelevantes; c_gamma e H_i=topo da pluma (definicao do artigo, `mynn_uedmf_hmode=1`) efeito pequeno.
  Parametros da Tabela 2 foram ajustados para OUTRO esquema (area/TKE prognosticas): nao ha razao para esperar
  que se transfiram sem recalibracao.
  Chaves: `mynn_uedmf_c_eps|c_delta|c_lambda|beta|mu0|chi|c_gamma`, `mynn_uedmf_hmode`.
  Testes: `src/unified_edmf/test_uedmf_thermo.F90`, `test_uedmf_entdet_map.F90`, `test_uedmf_entdet_moist.F90`.

## Microfisica explicita do Thompson DENTRO da pluma (teste de viabilidade, `mynn_uedmf_mp=1`)
- A pluma carrega theta,qv,qc,qr,qi,qs,qg e chama uma microfisica de UMA camada por `uedmf_mp_step`, um PONTEIRO DE
  PROCEDIMENTO que o driver liga a `step_plume_thompson` (bl_mynn_subroutines.F90 nao tem `use` de nenhum modulo de
  microfisica). Ao MYNN volta o mesmo de antes: qt = qv + todo o condensado; TODO o condensado pelo canal qc.
  Ambiente: so qc (DMP_mf nao recebe qi/qr/qs/qg do ambiente). ni,nr zerados a cada chamada; dt=dz/w limitado a 300 s.
- mp=0 e bit a bit igual a antes. Diff vs upstream MYNN: +466 linhas, 1 original alterada.
- Resultado (GOAMAZON 35 h): profundidade do fluxo de massa ~igual a saturacao do MYNN (z01 mediano 14824 m nos dois),
  chuva +15% (1.16->1.33 mm), qg=0 e raios=0 na grade (particao nao chega a grade: so canal qc).
- Particao do condensado da pluma (por fluxo de massa): <8 km 93-100% qc liquido; 8-12 km 68% qc,14% qi... >12 km 93% neve.
  => canais novos (qr,qs,qg) so mudariam a parte da pluma acima de 8 km (~8% do fluxo de massa).
- Custo (-O2, 1 coluna, 3 repeticoes): MYNN original 0.97 s; pluma nova+saturacao MYNN 2.07 s; pluma nova+Thompson 3.23 s
  (533140 chamadas: 79% ativas ~4.7 us, 21% ociosas ~2.8 us). O Thompson SUBSTITUI a iteracao de saturacao (~1.1 s).
  Passo com plumas ativas: ~1.3 ms bruto/coluna (~6x o custo medio de um passo do SCM original). NAO medido em 3D.
- Nao testado: dupla contagem com a microfisica da grade; efeitos colaterais do estado de modulo compartilhado do
  Thompson sobre a grade com a opcao ligada; energia do calor de fusao ao entregar gelo como qc.

### Canais de especies COM BALANCO (mp=3)  [substitui o mp=2]
- `mynn_uedmf_mp=3`: hand-off ao solver = qt com TODA a agua, thl = theta_il de TODO o condensado, qc = so agua liquida de nuvem.
  A pluma exporta, por camada, a fonte da microfisica (variacao por especie dentro da camada, ponderada pelo fluxo de massa medio) e o
  aquecimento latente; a queda da precipitacao = -(soma das variacoes) (fechamento exato de agua da chamada), repartida em
  chuva/neve+gelo/graupel pelos diagnosticos do proprio Thompson (RAINNCV=total, SNOWNCV=neve+gelo, GRAUPELNCV=graupel) e depositada
  na mesma camada. Tudo escalado pelo fator do limitador de fluxo do MYNN.
- Driver (`apply_plume_balance`): especies externas += -dF_env/dz + producao; vapor += fonte - B_ext (devolve o que o solver pos no
  canal de vapor por carregar as especies em qt); qc += fonte; theta += aquecimento latente + calor que o theta_il escondeu nas especies.
- Verificacao (emulando o solver, forma de fluxo explicita): agua fecha na coluna a 7e-17; theta_il a 0.55%; calor latente do
  Thompson vs constantes do MYNN difere 2.4% (limite de exatidao da energia); Thompson conserva agua na camada a 5e-16; massa criada
  por corte em zero 0.062 kg/m2 (0.5% dos termos; 3.5% da chuva) contra 0.49 kg/m2 (14%) no mp=2.
- Resultado (GOAMAZON 35 h, Thompson na grade): qg(-40..0C) 2.6e-5 (WSM6 na grade: 2.4e-5); raios 32.0 em 2903 instantes seguindo o ciclo
  diurno (rajadas 11-17h dia 1, 08-20h dia 2; nenhuma a noite; 52% com pluma ativa, resto = graupel residual); chuva 1.77 mm.
  MAGNITUDE DOS RAIOS SEM SIGNIFICADO: coeficientes do ramo EP do modulo de raios sao PLACEHOLDERS ('calibrar depois').
- Nao feito: solver implicito real (a verificacao emula), calibracao dos raios, dupla contagem com a microfisica da grade,
  reconciliar constantes de calor latente, queda pelas camadas abaixo (hoje: mesma camada).
- mp=2 (canais sem fonte, cria massa) fica no codigo SO para reproduzir o relatorio.

## Saídas (NetCDF)
Perfil: theta, qv..qg, pressure, tke, mixing_length, exch_h/m, cldfra_bl, qc_bl, qi_bl,
mf_a, mf_w. Escalares: rainnc(v), snownc(v) [equivalente líquido, NÃO profundidade],
flash_rate(_cg/_ic), pblh. Resumo do ensemble: NetCDF com dimensão `member`.

## Erros de integração já corrigidos (não reintroduzir)
1. **Thompson nunca executava**: `ide=jde=1` => `i_end=j_end=0` no laço de colunas do
   `mp_gt_driver`. Usar `ide=jde=2` (laço vai até `ide-1`).
2. `mp_gt_driver` em modo MPAS exige `rainprod`/`evapprod`.
3. **WSM6 com calor latente errado**: passar `xlv0 = xlv` (2.5e6), NÃO `xlv0` (3.15e6).
4. MYNN: chamar `bl_mynn_init` (constantes compartilhadas, ex. karman) e usar
   `bl_mynn_output=.true.` para que os diagnósticos do EDMF sejam copiados.
Conclusão retratada: o atalho `no_micro` do Thompson NÃO bloqueia parcelas quentes
supersaturadas (acima de 0 C, qvsi=qvs).

## Checklist ao integrar/alterar chamadas de física
- Auditar TODOS os argumentos contra a chamada de referência do MPAS
  (`src/mpas_atmphys_driver_microphysics_REFERENCE.F`).
- Um resultado nulo de um esquema = suspeita de erro de integração até provar que a
  rotina executou (contador/print dentro do laço).
- Consistência de saturação: RH sobre água em células quentes deve ficar ~<= 101% após o passo.

## Pendências
- Instrumentar os outros caminhos de chuva (coleta, derretimento); a autoconversão explica a direção, não a magnitude.
- Rodar ensemble/caso andino (dados ERA5, chave CDS do usuário).
- Validar o protótipo de EDMF contra LES (seco/raso/profundo).

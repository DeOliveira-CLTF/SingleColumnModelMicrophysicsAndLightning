# Esquemas e adaptações

Referências e requisitos originais: `configs/lightning/schemes.yaml`.
As adaptações têm sufixos LS/EDMF e `1D_adaptado`, sem substituir entradas CRM.

Exportam-se PR92-wmax, PW10, IFluxT, McCaul F1/F2/F3, LPI 1D e CAPE×P.
No EDMF, wmax é o máximo vertical da média de plumas, não o máximo entre plumas.
Limiares estritos: PW10 `w>10 m/s`, IFluxT `w>1 m/s` a 260 K. McCaul F1
interpola velocidade e graupel separadamente a 258,15 K, com graupel em g/kg.
F2 integra gelo total com densidade seca, sem fator 1000.

LPI omite a condição horizontal de vizinhança; usa média em profundidade e
máscara térmica no centro das camadas. Razões 0/0 recebem zero; ausência da
camada térmica recebe NaN. Hidrometeoros da coluna são apenas uma hipótese
de representatividade das plumas nas versões exploratórias EDMF.

CAPE usa [MetPy surface_based_cape_cin](https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.surface_based_cape_cin.html),
parcela pseudoadiabática líquida na célula simulada inferior. CAPE no início
da janela multiplica a precipitação total dos 15 min seguintes. Hipótese e
janela diferem de Romps: identificador `CAPExP_SBCAPE15_1D_adaptado`.

Números PR92/McCaul são escores antes da calibração SCM; não se atribui taxa
absoluta pelo coeficiente publicado. Calibração exploratória de controles fica
separada da calibração principal baseada na média dos 40 membros.

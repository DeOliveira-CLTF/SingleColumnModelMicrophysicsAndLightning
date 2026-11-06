# Interfaces diagnósticas

A saída `.vertical.nc` só lê argumentos, após MYNN e antes da aplicação das
tendências. O instante registrado é o início do passo físico.

```
rhoz(k) = [rho(k)*dz(k+1)+rho(k+1)*dz(k)]/[dz(k+1)+dz(k)]
w_EDMF(k) = mf_w(k)/rhoz(k)
```

`rho_MYNN` é exatamente a densidade passada pelo driver, com sua aproximação
de temperatura virtual. Não é a densidade seca usada para concentrações de
hidrometeoros. `mf_w` já foi dividido por `sum(UPA)` no MYNN: não dividir
novamente por `mf_a`, que contém `Psig_w`. Não se inventam velocidades de pluma
na superfície ou no topo do modelo. O MYNN recebe `w2=0`; `w_LS` é mantido
separado, sem soma de velocidades.

Isotermas: menor cruzamento descendente da temperatura da coluna. O alcance
exige interfaces ativas adjacentes; não afirma que temperatura de pluma e coluna
sejam iguais.

`uedmf_sp_flux`: canais 1=chuva, 2=gelo de nuvem, 3=neve, 4=graupel,
em `kg m-2 s-1`; índice k+1 corresponde à interface física k. O preenchimento
exige `edmf_type=1`, callback associado e `uedmf_mp_mode>=2`. No experimento
principal type=0/mode=0, zeros de inicialização não são fluxos disponíveis.
O modo 3 não foi ativado.

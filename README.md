# IPS: N corpos + FFT (versão enxuta)

Versão em Fortran 90 dos programas `ips_thalassa.f`, `ips_n_naiad.f`,
`ips_j_am.f` e `ssatfftm-2SEMFFT-2.f`. A física é a mesma: coordenadas
heliocêntricas canônicas, J2/J4, integrador RA15 e FFT das séries temporais.
Tudo o que antes mudava dentro do código fica agora em dois arquivos de
entrada:

| arquivo       | o que tem                                                        |
|---------------|------------------------------------------------------------------|
| `ips.in`      | planeta central, passo, nº de pontos, **variação de a e e**, ângulos ressonantes, limiares de saída |
| `planetas.in` | um satélite por linha: `nome massa(kg) a(km) e I w Om M` (ângulos em graus) |

## Compilar e rodar

```bash
gfortran -O2 -o ips ips.f90
./ips                 # lê ips.in do diretório atual
./ips outro.in        # ou outro arquivo de configuração
```

Versão paralela (opcional): cada ponto da grade vai para um núcleo.

```bash
gfortran -O2 -fopenmp -o ips ips.f90
OMP_NUM_THREADS=4 ./ips
```

## Variação (aa e ee)

No `ips.in`:

```
1                                 ! corpo variado (ordem no planetas.in)
48260.0d0   48300.0d0   300       ! a_ini  a_fim (km)   n. de intervalos
1.40d-3     1.40d-3     0         ! e_ini  e_fim        n. de intervalos
```

* Com `0` intervalos, só o valor inicial é usado.
* Se os dois tiverem intervalos, a grade é 2D: para cada `e`, varre todos os `a`.
* Os valores de `a` e `e` do corpo variado no `planetas.in` são ignorados.

## Ângulos ressonantes

```
1  2                     ! corpo A  corpo B
-69  73   0  0   -4   0  ! res1
-69  73   0  0    0  -4  ! res2
```

`phi = k1*lamA + k2*lamB + k3*varpiA + k4*varpiB + k5*OmA + k6*OmB`,
com `lam = M + w + Om` e `varpi = w + Om`.

## Saída

São 5 séries: `res1`, `res2`, `inc`, `a` e `e`. As três últimas são do corpo
indicado em "SAIDA". Para cada série e cada limiar há um arquivo
`<serie>_<limiar>.dat`, por exemplo `res1_8.dat` (equivale ao antigo
`edh8RES1.dat`) e `e_0.5.dat` (equivale ao `edh05e2.dat`). As colunas são:

```
a(km)   e   periodo(dias)
```

Cada linha é um pico do espectro com amplitude maior que `limiar`% da
amplitude máxima. A máxima é procurada entre os picos com período maior que
o "período mínimo". Os arquivos são gravados na ordem da grade à medida que
cada ponto termina, então dá para plotar com a rodada ainda em andamento.
Exemplo no gnuplot:

```
plot 'res1_8.dat' u 1:3 w d
```

## Gráficos

```bash
pip install numpy matplotlib
python3 graf_ips.py                  # todos os .dat da pasta atual
python3 graf_ips.py pasta            # .dat de outra pasta
python3 graf_ips.py . res1_8 e_0.5   # só alguns arquivos
```

Os PNGs vão para `graficos_ips/`. O eixo x é escolhido sozinho: `a` se só
`a` variou, `e` se só `e` variou, e na grade 2D sai um gráfico por valor
de `e` (`res1_8_e0.002.png`, ...). O nome do corpo e o número de clones do
título vêm do `ips.in`. Também lê os `edh*.dat` antigos (2 colunas).

## Exemplos

`exemplos/` tem os três sistemas dos códigos antigos. Para rodar um deles:

```bash
cd exemplos/saturno_anthe && ../../ips
```

* `netuno_naiad`: Naiad–Thalassa 73:69 (igual ao `ips.in` da raiz)
* `jupiter_amalthea`: Amalthea (a ressonância ainda precisa ser ajustada, veja o comentário no arquivo)
* `saturno_anthe`: Anthe–Mimas 11:10

## Diferenças em relação aos códigos antigos

* Compila com `gfortran` sem `-std=legacy`. O original com `-O2` travava por
  causa de variáveis não inicializadas.
* A resolução para os momentos (antes `gaussj` com `big`, `dum` e `pivinv` em
  precisão simples) virou fórmula fechada em dupla precisão. Com o `gaussj`
  corrigido para dupla, o original e este código dão as mesmas séries
  (diferença ~1e-14 nos elementos) e os mesmos picos da FFT.
* A FFT usa os fatores `wr` e `wi` em dupla (o original fazia `sngl`).
* A `FORCE` não chama mais a `ENTRE` a cada passo. Ficou ~20x mais rápido.
* `N` e o tamanho da série vêm dos arquivos de entrada, sem dimensões fixas
  (o `ips_j_am.f` tinha `N=7` com vetores de tamanho 5).

"""
Plota os IPS gerados pelo ips.f90.

Le todos os arquivos <serie>_<limiar>.dat (colunas: a(km)  e  periodo(dias))
e salva um PNG para cada um em graficos_ips/.

  - so a variou          -> eixo x = a (km)
  - so e variou          -> eixo x = e
  - a e e variaram (2D)  -> um grafico por valor de e, com eixo x = a

Uso:
  python3 graf_ips.py                 # arquivos .dat do diretorio atual
  python3 graf_ips.py pasta_da_rodada # ou de outra pasta
  python3 graf_ips.py . res1_8 e_0.5  # so alguns arquivos
"""
import os
import re
import sys
from glob import glob

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import StrMethodFormatter

NOMES = {
    'res1': 'res1',
    'res2': 'res2',
    'inc': 'inclinação',
    'a': 'semieixo maior',
    'e': 'excentricidade',
}


def ler_linhas_uteis(arquivo):
    """Linhas de dados de um .in (sem comentarios e linhas vazias)."""
    with open(arquivo) as f:
        for linha in f:
            linha = linha.strip()
            if linha and linha[0] not in '!#':
                yield linha


def info_rodada(pasta):
    """Nome do corpo variado e tamanho da grade (n_a, n_e), lidos do
    ips.in + planetas.in. Devolve (None, None, None) se nao achar."""
    try:
        linhas = list(ler_linhas_uteis(os.path.join(pasta, 'ips.in')))
        arq_pla = linhas[0].split()[0]
        ivar = int(linhas[8].split()[0])
        n_a = int(linhas[9].split()[2]) + 1
        n_e = int(linhas[10].split()[2]) + 1
        corpos = list(ler_linhas_uteis(os.path.join(pasta, arq_pla)))
        return corpos[ivar - 1].split()[0], n_a, n_e
    except (OSError, IndexError, ValueError):
        return None, None, None


def ler_dat(arquivo):
    """Devolve a, e, periodo. Aceita tambem o formato antigo (a, periodo)."""
    try:
        d = np.loadtxt(arquivo, ndmin=2)
    except ValueError:
        return None
    if d.size == 0:
        return None
    if d.shape[1] == 2:                       # edh*.dat antigo
        return d[:, 0], np.zeros(len(d)), d[:, 1]
    return d[:, 0], d[:, 1], d[:, 2]


def plotar(x, per, xlabel, titulo, saida, xfmt):
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(x, per, 'ko', markersize=2, alpha=0.6)
    ax.set_yscale('log')
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel('Período (dias)', fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    ax.set_title(titulo, fontsize=18)
    ax.grid(True, linestyle='-', alpha=0.5)
    if x.min() < x.max():
        ax.set_xlim(x.min(), x.max())
        ax.set_xticks(np.linspace(x.min(), x.max(), 6))
    ax.set_ylim(per.min(), per.max())
    ax.xaxis.set_major_formatter(StrMethodFormatter(xfmt))
    plt.tight_layout()
    plt.savefig(saida, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('  ->', saida)


def main():
    pasta = sys.argv[1] if len(sys.argv) > 1 else '.'
    escolhidos = sys.argv[2:]
    saida_dir = os.path.join(pasta, 'graficos_ips')
    os.makedirs(saida_dir, exist_ok=True)

    corpo, grade_a, grade_e = info_rodada(pasta)
    de_quem = f' de {corpo}' if corpo else ''

    arquivos = sorted(glob(os.path.join(pasta, '*.dat')))
    if escolhidos:
        arquivos = [f for f in arquivos
                    if os.path.splitext(os.path.basename(f))[0] in escolhidos]
    if not arquivos:
        print('Nenhum arquivo .dat encontrado em', pasta)
        return

    for arquivo in arquivos:
        nome = os.path.splitext(os.path.basename(arquivo))[0]
        dados = ler_dat(arquivo)
        if dados is None:
            print(nome, ': vazio, pulando')
            continue
        a, e, per = dados

        m = re.match(r'(.+)_([0-9.]+)$', nome)
        serie, lim = (m.group(1), m.group(2)) if m else (nome, None)
        rotulo = NOMES.get(serie, serie)
        if lim:
            rotulo += f'  (picos > {lim}% da amp. máx.)'

        na, ne = len(np.unique(a)), len(np.unique(e))
        print(nome, f': {len(a)} picos, {na} valores de a, {ne} de e')

        if ne == 1 or (na > 1 and ne > 1):
            # eixo x = a  (um grafico por e se a grade for 2D)
            for ev in np.unique(e):
                sel = e == ev
                titulo = f'IPS de {grade_a or na} clones{de_quem} – {rotulo}'
                arq = nome
                if ne > 1:
                    titulo += f'\ne = {ev:.6g}'
                    arq += f'_e{ev:.6g}'
                plotar(a[sel], per[sel], 'Semieixo maior inicial (km)', titulo,
                       os.path.join(saida_dir, arq + '.png'), '{x:.0f}')
        else:
            # so e variou: eixo x = e
            titulo = f'IPS de {grade_e or ne} clones{de_quem} – {rotulo}\na = {a[0]:.6g} km'
            plotar(e, per, 'Excentricidade inicial', titulo,
                   os.path.join(saida_dir, nome + '.png'), '{x:.4g}')


if __name__ == '__main__':
    main()

import numpy as np
import matplotlib.pyplot as plt

from genbank import get_orfs

def plot_orfs(var_start: int = 143509, var_end: int = 143509, size=15):
    if var_start-var_end == 0:
        bt_lim = var_start-size*10
        tp_lim = var_start+size*10
    else:
        bt_lim = var_start-(var_end-var_start)*0.05
        tp_lim = var_end+(var_end-var_start)*0.05

    orfs = get_orfs()
    fig, ax = plt.subplots(figsize=(15, 2))
    max_end = max([end for start, end in orfs.values()])
    zeros = np.zeros(max_end)
    global_matrix = [zeros.copy()]
    for orf, (start, end) in orfs.items():
        i = 0
        while global_matrix[i][start:end].any():
            i += 1

            if i == len(global_matrix)-1: break

        if i == len(global_matrix)-1:
            global_matrix.append(zeros.copy())

        global_matrix[i][start:end] = 1

        if end>bt_lim and start<tp_lim:
            ax.plot([start, end], [i, i], label=orf)

    ax.set_xlabel('Position')
    ax.set_yticks([])
    ax.set_title(f'{["Déletion", "Insertion"][var_start == var_end]} de {[var_end - var_start, size][var_start == var_end]}pb à la position {var_start}')

    if var_start != var_end:
        ax.vlines(var_start, 0, len(global_matrix), color='red', linestyle='--')
        ax.vlines(var_end, 0, len(global_matrix), color='red', linestyle='--')
    else:
        ax.vlines(var_start, 0, len(global_matrix), color='green', linestyle='--')

    plt.xlim(bt_lim, tp_lim)
    plt.legend( loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    plt.tight_layout()
    plt.show(dpi=300)

if __name__ == '__main__':
    plot_orfs()
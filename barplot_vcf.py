import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps
import pandas as pd

from vcf import is_above
import scipy.stats as stats

def print_sv():
    # Chargement des données
    data = pd.read_csv(f'data/vcf/bilanfiltered+notpassing+notprecise.csv')

    hot = []
    cold = []
    steps = []

    for experiment in data['experiment'].unique():
        hot.append([])
        cold.append([])
        steps.append(experiment)

        for replicate in data['replicate'].unique():
            hot[-1].append(0)
            cold[-1].append(0)
            choc = ["Choc froid", "Choc chaud"][replicate>5]
            subset = data[(data['experiment'] == experiment) & (data['replicate'] == replicate)]
            for indel in subset['indel'].unique():
                for i in range(len(subset)):
                    if choc == 'Choc froid':
                        cold[-1][-1] += 1
                    else:
                        hot[-1][-1] += 1

    hot = [x[5:] for x in hot]
    cold = [x[:5] for x in cold]

    fig, ax = plt.subplots()
    # boxplots with cold and hot close together
    ax.boxplot(cold[0]+hot[0], positions=[0], widths=0.6, patch_artist=True, boxprops=dict(facecolor='green'))
    ax.boxplot([0]+cold[1:], positions=np.array(range(len(cold)))*2.0-0.3, widths=0.6, patch_artist=True, boxprops=dict(facecolor='blue'))
    ax.boxplot([0]+hot[1:], positions=np.array(range(len(hot)))*2.0+0.3, widths=0.6, patch_artist=True, boxprops=dict(facecolor='red'))
    ax.set_xticks(range(len(steps)*2))
    ax.set_xticklabels([item for pair in zip(steps, [''] * len(steps)) for item in pair])
    ax.set_xlim(-2, len(steps)*2)
    ax.set_ylabel('Nombre de variants')
    ax.set_title(f'Répartition des variants chauds et froids')

    plt.ylim(0)
    plt.xlim(-1, len(steps)*2-1)
    # remove half the x gradations
    ax.set_xticks(np.arange(0, len(steps)*2, 2))

    plt.tight_layout()
    plt.savefig(f'image/distribution_variants-sv.png', dpi=300)

    for h, c, step in zip(hot[1:], cold[1:], list(data.keys())[1:]):
        if step == "P30":
            print(step, stats.mannwhitneyu(h, c)[1])
        else:
            print(step, stats.ttest_ind(h, c)[1])

        fig, axes = plt.subplots(2, 1, sharex=True)
        for plot, ax in enumerate(axes):
            x = [i for i in range(5)]
            points = []
            for n in range(5):
                if plot == 1:
                    points.append(hot[n])
                else:
                    points.append(cold[n])

            cmap = colormaps['Dark2']
            for n, replicata in enumerate(points):
                ax.plot(x, replicata, label=f"Réplicat {5*plot+n+1}", c=cmap(n/5))

            ax.set_xticks(x)
            ax.set_xticklabels(range(1, 6))
            ax.set_ylabel('Nombre de variants')
            ax.set_title(f'Répartition des variants{" chauds" if plot else " froids"}')
            ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.xlim(0, len(steps)-1)
        plt.ylim(0, 50)
        plt.tight_layout()
        plt.savefig('image/cinetique_variants.png', dpi=300)

############################
# DELETION FUZE DEUX ORFS? #
############################

def print_snp():
    data = pd.read_csv('data/vcf/bilanfiltered_snp.csv')

    hot = []
    cold = []
    steps = []

    for experiment in data['experiment'].unique():
        hot.append([])
        cold.append([])
        steps.append(experiment)

        for replicate in data['replicate'].unique():
            hot[-1].append(0)
            cold[-1].append(0)
            choc = ["Choc froid", "Choc chaud"][replicate > 5]
            subset = data[(data['experiment'] == experiment) & (data['replicate'] == replicate)]

            for i in range(len(subset)):
                if choc == 'Choc froid':
                    cold[-1][-1] += 1
                else:
                    hot[-1][-1] += 1

    hot = [x[5:] for x in hot]
    cold = [x[:5] for x in cold]

    fig, ax = plt.subplots()
    # boxplots with cold and hot close together
    ax.boxplot(cold[0] + hot[0], positions=[0], widths=0.6, patch_artist=True, boxprops=dict(facecolor='green'))
    ax.boxplot([0] + cold[1:], positions=np.array(range(len(cold))) * 2.0 - 0.3, widths=0.6, patch_artist=True,
               boxprops=dict(facecolor='blue'))
    ax.boxplot([0] + hot[1:], positions=np.array(range(len(hot))) * 2.0 + 0.3, widths=0.6, patch_artist=True,
               boxprops=dict(facecolor='red'))
    ax.set_xticks(range(len(steps) * 2))
    ax.set_xticklabels([item for pair in zip(steps, [''] * len(steps)) for item in pair])
    ax.set_xlim(-2, len(steps) * 2)
    ax.set_ylabel('Nombre de variants')
    ax.set_title(f'Répartition des variants chauds et froids')

    plt.ylim(0)
    plt.xlim(-1, len(steps) * 2 - 1)
    # remove half the x gradations
    ax.set_xticks(np.arange(0, len(steps) * 2, 2))

    plt.tight_layout()
    plt.savefig(f'image/distribution_variants-snp.png', dpi=300)

    for h, c, step in zip(hot[1:], cold[1:], list(data["experiment"].unique())[1:]):
        if step in ["P50", "P30"]:
            print(step, stats.mannwhitneyu(h, c)[1])
        else:
            print(step, stats.ttest_ind(h, c)[1])

if __name__ == '__main__':
    print_sv()
    print_snp()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from vcf import neg_exp_func, is_above
from scipy.optimize import curve_fit


def make_flat_plots(dataset="bilanfiltered+notpassing+notprecise+notSAF"):
    data = pd.read_csv(f'data/vcf/{dataset}.csv')

    # Forme des points
    forms = ['o', 's', 'D', 'x', '+']

    for experiment in data['experiment'].unique():
        step_stats = {
            'Choc froid': {
                'INS': {
                    'above': [],
                    'below': []
                },
                'DEL': {
                    'above': [],
                    'below': []
                }
            },
            'Choc chaud': {
                'INS': {
                    'above': [],
                    'below': []
                },
                'DEL': {
                    'above': [],
                    'below': []
                }
            }
        }
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))
        for replicate in data['replicate'].unique():
            choc = ["Choc froid", "Choc chaud"][replicate > 5]
            choice_index = ['Choc froid', 'Choc chaud'].index(choc)

            subset = data[(data['experiment'] == experiment) & (data['replicate'] == replicate)]

            for indel_ax in range(2):
                indel_subset = subset[data['indel'] == ['INS', "DEL"][indel_ax]]
                x = indel_subset['dv']
                y = indel_subset['af']
                if x.empty or y.empty: continue
                if sum(x) == 0 or sum(y) == 0: continue
                axes[choice_index][indel_ax].scatter(x, y,
                                                     marker=forms[replicate%5],
                                                     label=f"Réplicat {replicate}",
                                                     c='b', alpha=0.5)
                step_stats[choc][['INS', 'DEL'][indel_ax]]['above'].append(sum([is_above(s, a) for s, a in zip(x, y)]))
                step_stats[choc][['INS', 'DEL'][indel_ax]]['below'].append(len(y) - step_stats[choc][['INS', 'DEL'][indel_ax]]['above'][-1])

        x = np.linspace(1, 10000, 10000)
        x_data = [30, 100, 1000]
        y_data = [1.0, 0.1, 0.05]
        popt = curve_fit(neg_exp_func, x_data, y_data, p0=(1, 0.01, 0))[0]

        for n, ax in enumerate(axes.flatten()):
            if n>1: ax.set_xlabel('Support')
            if n<2: ax.set_title(["Insertion", 'Délétions'][n], fontsize=20)
            if n%2 == 0: ax.set_ylabel('Fréquence allélique')
            if n%2 == 1: ax.set_ylabel(["CHOC\nFROID", "CHOC\nCHAUD"][n//2],
                                       labelpad=30, fontsize=20, rotation=0,
                                       color=['blue', 'red'][n//2])


            ax.plot(x, neg_exp_func(x, *popt), 'r--')

            ax.set_xscale('log')
            ax.set_yscale('log')

            ax.set_xlim(10, 10000)
            ax.set_ylim(1e-3, 1)
            ax.grid(True)

            ax.legend(loc='lower right')

        plt.tight_layout()
        plt.savefig(f"image/{experiment}_sv.png", dpi=300)

        print(f"Nombre de variants pour {experiment}:")
        for choc, indel_dict in step_stats.items():
            print(f"  {choc}:")
            for indel, above_below_dict in indel_dict.items():
                print(f"    {indel}:")
                for above_below, values in above_below_dict.items():
                    print(f"      {above_below}: {sum(values)}")
                # ratio
                print(f"      Ratio: {above_below_dict['above'][0]/(above_below_dict['above'][0]+above_below_dict['below'][0])}")


if __name__ == '__main__':
    make_flat_plots()
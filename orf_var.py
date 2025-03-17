import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colormaps
from genbank import get_targeted_orfs


def hist_of_var_on_orfs(path=''):
    data = pd.read_csv('data/vcf/bilanfiltered.csv')
    colors = colormaps['tab10'].colors  # 10 colors

    for experiment in data['experiment'].unique():  # For each experiment
        rows = 1 if experiment == "P15" else 2  # P15 does not differentiate between cold and hot shocks
        fig, axes = plt.subplots(rows, 1, figsize=(20, 5))

        if isinstance(axes, np.ndarray): axes = [ax for ax in axes]  # Convert to list
        if type(axes) != list: axes = [axes]  # Convert to list

        ax_to_chose = [int(n>4) for n in range(10)] if experiment != "P15" else [0]*10  # 1-5 -> 0, 6-10 -> 1 except for P15 where all are 0
        bar_heights = [[0]*156, [0]*156] if experiment != "P15" else [[0]*156]  # Makes 2 lists of 156 zeros if experiment is not P15, else 1 list of 156 zeros

        for replicate in data['replicate'].unique():  # For each replicate
            orfs = [0]*156  # 156 ORFs
            subset = data[(data['experiment'] == experiment) & (data['replicate'] == replicate)]  # Get the subset of the data
            indexation = ax_to_chose[replicate-1]  # Get the index of the axes to use

            for i in range(len(subset)):
                orf = get_targeted_orfs(subset.iloc[i]['start'], subset.iloc[i]['end'])  # Get the ORFs targeted by the variant
                for name, start, end in orf:
                    name = name.replace('ORF', '').split('_')[0]
                    orfs[int(name)-1] += 1  # Increment the number of variants for the ORF

            axes[indexation].bar(range(1, 157),
                                 orfs,
                                 color=colors[replicate-1],
                                 label=f"Replicate {replicate}",
                                 bottom=bar_heights[indexation])

            bar_heights[indexation] = [bar_heights[indexation][i] + orfs[i] for i in range(156)]  # Update the bar heights
            axes[indexation].set_ylabel('Number of variants')
            axes[indexation].legend(loc='center left', bbox_to_anchor=(1, 0.5))

        axes[0].set_title(experiment)
        axes[-1].set_xlabel('ORF')
        max_bar_heights = 0
        for bar_height in bar_heights: max_bar_heights = max(max_bar_heights, max(bar_height))
        for ax in axes: ax.set_ylim(0, max_bar_heights)

        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        if path == '': plt.show(dpi=300)
        else: plt.savefig(f"{path}/{experiment}-barplot.png", dpi=300)



def repr_orfs():
    data = pd.read_csv('data/vcf/bilanfiltered.csv')

    sns.catplot(x='replicate', hue='indel', data=data, kind='count', col="experiment", col_wrap=3)

    plt.show()


if __name__ == '__main__':
    hist_of_var_on_orfs("./orf_analyses")

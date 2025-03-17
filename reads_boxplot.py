import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

def get_flagstats(path: str) -> dict:
    with open(path, 'r') as f:
        lines = f.readlines()
    return {
        'n_reads': int(lines[0].split(' ')[0]),  # Ex: 103741 + 0 in total (QC-passed reads + QC-failed reads)
        'n_mapped': int(lines[6].split(' ')[0]), # Ex: 88566 + 0 mapped (85.37% : N/A)
    }
experiments = ["P90", "P65", "P50", "P30", "P15"]
n_reads = []
n_mapped = []

for experiment in experiments:
    n_reads.append([])
    n_mapped.append([])
    for replicate in range(1, 11):
        flagstat = get_flagstats(f'data/flagstats/{experiment}/{experiment}.{replicate}.trimed1000.flagstat')
        n_reads[-1].append(flagstat['n_reads'])
        n_mapped[-1].append(flagstat['n_mapped'])

fig, ax = plt.subplots()
ax.boxplot(n_mapped, positions=np.array(range(len(experiments))))

ax.set_xticks(range(len(experiments)))
ax.set_xticklabels([item for item in experiments])
ax.set_ylabel('Nombre de reads mappés')
plt.ylim(0)
plt.tight_layout()
plt.savefig('image/reads_boxplot.png', dpi=300)

# Différence entre P90 et P65, P50, P30, P15
print(mannwhitneyu(n_mapped[0], [item for sublist in n_mapped[1:] for item in sublist])[1])
# >>> 3.545955911799725e-05
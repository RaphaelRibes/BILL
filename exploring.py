import pandas as pd
from genbank import get_targeted_orfs


def export_mut(orf, experiments=['P90', 'P65', 'P50', 'P30', 'P15'], replicates=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]):
    dt = pd.read_csv("data/vcf/bilanfiltered+notpassing+notprecise.csv")
    final = ""
    for experiment in experiments:
        subset = dt[(dt['experiment'] == experiment)]
        for _, line in subset.iterrows():
            if f"ORF{orf}" in [o for o, _, _ in get_targeted_orfs(int(line['start']), line['end'])] and line['replicate'] in replicates:
                header = f">ORF{orf}_{experiment}"
                header += f"({line['replicate']})"
                header += ['+', '-'][line['indel'] == 'DEL']
                header += str(line['start'])
                header += f"<>{line['end']}" if line['end'] != line['start'] else ''

                final += (header + '\n' + line['vseq'] + '\n')

    with open(f'{'_'.join(experiments)}_ORF{orf}.fasta', 'w') as f:
        f.write(final)

if __name__ == "__main__":
    export_mut(136)  # Exporte les mutations de l'ORF136 des expériences P15 à P90 pour tous les réplicats
    export_mut(136, ['P90'])  # Exporte les mutations de l'ORF136 de l'expérience P90 pour tous les réplicats
    export_mut(136, ['P90'], [1, 2, 3])  # Exporte les mutations de l'ORF136 de l'expérience P90 pour les réplicats 1, 2 et 3
    export_mut(136, ['P90', 'P65'], [1, 2, 3])  # Exporte les mutations de l'ORF136 des expériences P90 et P65 pour les réplicats 1, 2 et 3
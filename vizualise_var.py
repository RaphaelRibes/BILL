from genbank import get_orfs
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import SeqIO, Seq
from tqdm.auto import tqdm
import os

def check_lists(x, y):
    # Check if x == y in any order
    if sorted(x) == sorted(y):
        return True

    # Check if x contains y
    if all(item in y for item in x):
        return True

    return False

def reverse_complement(seq):
    return seq.translate(str.maketrans('ATCG', 'TAGC'))[::-1]

def make_mut(header, orf_seq, line, snp, real_start, real_end, l_indel, l_len):
    insertion = ""
    if snp:
        header += f'~{real_start}'
        assert orf_seq[real_start] == line['ref'], f"Reference sequence does not match the reference sequence in the VCF file: {orf_seq[real_start]} != {line['ref']}"

        mutated_seq = orf_seq[:real_start - 1] + line['alt'] + orf_seq[real_start:]
    elif l_indel == 'INS':
        header += f'+{l_len}>{real_start}'
        insertion = line['vseq']
        if insertion == '<INS>':
            insertion = 'N'*l_len

        if real_start <= 0:
            mutated_seq = insertion + orf_seq
        elif real_end >= len(orf_seq):
            mutated_seq = orf_seq + insertion
        else:
            mutated_seq = orf_seq[:real_start] + insertion + orf_seq[real_end:]

    else:
        header += f'-{real_start}<>{real_end}'
        if real_start <= 0:
            mutated_seq = orf_seq[real_end:]
        elif real_end >= len(orf_seq):
            mutated_seq = orf_seq[:real_start]
        else:
            mutated_seq = orf_seq[:real_start] + orf_seq[real_end:]

    return header, mutated_seq, insertion

def get_subsequence(refseq, start: int, end: int):
    if end < start:
        start, end = end, start
    return str(refseq[start-1:end])

class AutoVar:
    def __init__(self, orf: str, snp = False, filters = None):
        self.snp = snp

        self.dataset = pd.read_csv('data/vcf/bilanfiltered_snp.csv') if self.snp else pd.read_csv('data/vcf/bilanfiltered+notpassing+notprecise.csv')
        self.reference_seq = SeqIO.read('data/KHV-U_trunc.fasta', 'fasta')


        orfs = get_orfs()
        orf = [orf] + list(orfs[orf])

        self.start = int(orf[1])
        self.end = int(orf[2])
        self.orf = orf[0]
        if orf[3]:
            self.orf_seq = reverse_complement(get_subsequence(self.reference_seq.seq, self.start, self.end))
        else:
            self.orf_seq = get_subsequence(self.reference_seq.seq, self.start, self.end)

        self.mutated = self._get_mutated_orfs()
        self.export_path = f"./mutated_orfs/{"snp" if self.snp else "sv"}/{self.orf}"

        if filters is not None:
            if filters['presence'] is not None:
                if not check_lists(filters['presence'], self.mutated['experiment'].unique()):
                    raise AssertionError(f'Variation not found in all experiments ({', '.join(filters["presence"])})')

        # create export path if it doesn't exist
        if not os.path.exists(self.export_path):
            os.makedirs(self.export_path)

    def _get_mutated_orfs(self):
        results = {
            'experiment': [],
            'replicate': [],
            'indel': [],
            'start': [],
            'end': [],
            'header': [],
            'mut seq': [],
            'vseq': [],
            'length': []
        }

        for experiment in ['P90', 'P65', 'P50', 'P30', 'P15']:
            for replicate in range(1, 11):
                subset = self.dataset[(self.dataset['experiment'] == experiment) & (self.dataset['replicate'] == replicate)]

                for _, line in subset.iterrows():
                    if self.snp:
                        mut_start = mut_end = line['pos']
                        l_len = 1
                        l_indel = 'SNP'
                    else:
                        mut_start = line['start']
                        mut_end = line['end']
                        l_len = line["length"]
                        l_indel = line['indel']

                    if mut_end < self.start or mut_start > self.end:
                        continue

                    header = f">{self.orf}_{experiment}({replicate})"
                    real_start = mut_start - self.start
                    real_end = mut_end - self.start

                    assert real_start >=0 or real_end >= 0, f"Real start or end is negative: {real_start}, {real_end}"

                    try:
                        header, mutated_seq, insertion = make_mut(header, self.orf_seq, line, self.snp, real_start, real_end, l_indel, l_len)
                    except:
                        continue

                    results['experiment'].append(experiment)
                    results['replicate'].append(replicate)
                    results['indel'].append(l_indel)
                    results['start'].append(mut_start)
                    results['end'].append(mut_end)
                    results['header'].append(header)
                    results['mut seq'].append(''.join([c for c in mutated_seq if c.isalpha()]))
                    results['vseq'].append(insertion)
                    results['length'].append(l_len)

        return pd.DataFrame(results)

    def plot_mutations(self, save=False):
        for indel in self.mutated['indel'].unique():
            subset = self.mutated[(self.mutated['indel'] == indel)]
            fig, ax = plt.subplots(1, 1, figsize=(12, 12))

            if subset['length'].unique().size > 1 or subset['start'].unique().size > 1:
                sns.catplot(subset, x='experiment', hue='replicate', col='start', kind='count', row='length')
            else:
                sns.countplot(subset, x='experiment', hue='replicate', ax=ax)
                ax.set_title(subset['indel'].unique()[0])

            if save:
                plt.savefig(f'{self.export_path}/{self.orf}_{indel}.png')
            else:
                plt.show()
            plt.close(fig)

    def to_fasta_DNA(self):
        fastafile = f">{self.orf}\n{self.orf_seq}\n"
        for _, line in self.mutated.iterrows():
            fastafile += f"{line['header']}\n{line['mut seq']}\n"

        return fastafile

    def to_fasta_AA(self):
        fastafile = f">{self.orf}\n{Seq.Seq(self.orf_seq).translate()}\n"
        for _, line in self.mutated.iterrows():
            # remove every character that is not a nucleotide
            try:
                mut_seq = Seq.Seq(line["mut seq"]).translate()
            except:
                mut_seq = "COULD NOT TRANSLATE"
            fastafile += f"{line['header']}\n{mut_seq}\n"

        return fastafile

    def export(self):
        if self.mutated.empty: return
        self.mutated.to_csv(f'{self.export_path}/mutations.csv')
        self.plot_mutations(save=True)

        with open(f'{self.export_path}/mutated_orfs_seq.fasta', 'w') as outfile:
            outfile.write(self.to_fasta_DNA())

        with open(f'{self.export_path}/mutated_orfs_seq_proteins.fasta', 'w') as outfile:
            outfile.write(self.to_fasta_AA())

        with open(f'{self.export_path}/summary.txt', 'w') as outfile:
            txt = f"ORF {self.orf}\nStart: {self.start}\nEnd: {self.end}\n\nn mutations: {len(self.mutated)}\n"
            if not self.snp:
                txt += f"n DEL: {len(self.mutated[self.mutated['indel'] != 'DEL'])}\n"
                txt += f"n DEL sizes: {len(self.mutated[self.mutated['indel'] != 'DEL']['length'].unique())}\n"
                txt += f"n INS: {len(self.mutated[self.mutated['indel'] != 'INS'])}\n"
                txt += f"n INS sizes: {len(self.mutated[self.mutated['indel'] != 'INS']['length'].unique())}\n"
            outfile.write(txt)

def make_real_seqs(orf, start, end, reverse):
    try:
        sv = pd.read_csv(f'mutated_orfs/sv/{orf}/mutations.csv')
    except:
        sv = None
    snp = pd.read_csv('data/vcf/bilanfiltered_snp.csv')

    virus_refseq = SeqIO.read('data/KHV-U_trunc.fasta', 'fasta').seq
    if reverse:
        refseq = reverse_complement(''.join([c for c in get_subsequence(virus_refseq, start, end) if c.isalpha()]))
    else:
        refseq = ''.join([c for c in get_subsequence(virus_refseq, start, end) if c.isalpha()])

    dna = f">ORF{orf}\n{refseq}\n"
    aa = f">ORF{orf}\n{Seq.Seq(refseq).translate()}\n"

    for experiment in ['P90', 'P65', 'P50', 'P30', 'P15']:
        for replicate in range(1, 11):
            mut_seq = refseq

            for mut in snp[(snp['experiment'] == experiment) & (snp['replicate'] == replicate) &
                           (snp['pos'] <= end) & (start <= snp['pos'])]['pos']:
                snp_val = snp[(snp['experiment'] == experiment) &
                                (snp['replicate'] == replicate) &
                                (snp['pos'] == mut)]['alt'].values[0].replace(',', '')
                mut_seq = mut_seq[:mut-start] + snp_val + mut_seq[mut-start+1:]

            if sv is not None:
                stvar = sv[
                    (sv['experiment'] == experiment) &
                    (sv['replicate'] == replicate) &
                    (sv['end'] <= end) &
                    (sv['start'] >= start)
                    ]
                if not stvar.empty:
                    for _, line in stvar.iterrows():
                        _, mut_seq, _ = make_mut("", mut_seq, line, False, line['start']-start, line['end']-start, line['indel'], line['length'])

            header = f">ORF{orf}_{experiment}({replicate})"
            dna += f"{header}\n{mut_seq}\n"
            aa += f"{header}\n{Seq.Seq(mut_seq).translate()}\n"

    with open(f'mutated_orfs/combined/{orf}_dna.fasta', 'w') as f:
        f.write(dna)

    with open(f'mutated_orfs/combined/{orf}_aa.fasta', 'w') as f:
        f.write(aa)

def generate_true_sequences():
    for orf, (start, end, reverse) in get_orfs().items():
        make_real_seqs(orf, start, end, reverse)

def generate_var():
    for snp in [False, True]:
        filters = {
            'presence': ['P15', 'P30', 'P50', 'P65', 'P90'],
        }

        orfs = [o for o in get_orfs()]
        for orf in tqdm(range(1, 157)):
            if orf < 9:
                for subgene in ["_1", "_2"]:
                    av = AutoVar(f"ORF{orf}{subgene}", snp)
                    av.export()
            else:
                if f"ORF{orf}" not in orfs: continue
                av = AutoVar(f"ORF{orf}", snp)
                av.export()
        sv_snp = 'snp' if snp else 'sv'

        for folder in os.listdir(f'mutated_orfs/{sv_snp}'):
            if len(os.listdir(f'mutated_orfs/{sv_snp}/{folder}/')) == 0: os.rmdir(f'mutated_orfs/{sv_snp}/{folder}/')

if __name__ == '__main__':
    generate_var()
    generate_true_sequences()
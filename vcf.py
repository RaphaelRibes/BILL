import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from tqdm.auto import tqdm

def vcf_to_df(path: str, save: bool = False):
    if not path.endswith('.vcf'):
        raise ValueError('The file must be a VCF file')

    with open(path, 'r') as f:
        for l in f:
            if l.startswith('##'):  # Commentaires
                continue
            elif l.startswith('#'):  # Headers
                header = l.strip().split('\t')
            else:
                break

    vcf = pd.read_csv(path, sep='\t', comment='#', header=None)
    vcf.columns = header

    path = path.replace('.vcf', '.csv')
    if save: vcf.to_csv(path, index=False)

    return vcf

def get_svinfo(df: pd.DataFrame, pos: int):
    dt = df.iloc[pos]['INFO'].split(';')
    final = {'precise': dt[0] == 'PRECISE',
             'PASS': df.iloc[pos]['FILTER'] == 'PASS',
             'POS': int(df.iloc[pos]['POS']),
             'AF': 0,
             'QUAL': int(df.iloc[pos]['QUAL'])}

    for d in dt[1:]:
        k, v = d.split('=')
        if k == "SUPPORT": continue  # On utilise DV qui est la meme chose
        try:
            if '.' in v:  # Si c'est un float
                final[k] = float(v)
            else:
                final[k] = int(v)
        except ValueError:
            final[k] = v

    summary = df.iloc[pos][df.columns[-1]].split(':')

    # Genotype
    if summary[0] == '0/0':
        final['genotype'] = 'Homozygous reference'
    elif summary[0] == '1/1':
        final['genotype'] = 'Homozygous alternative'
    elif summary[0] == '0/1' or summary[0] == '1/0':
        final['genotype'] = 'Heterozygous'
    elif summary[0] == './.':
        final['genotype'] = 'No call'
    else:
        print(summary[0])
        raise ValueError('Unknown genotype')

    final['GQ'] = int(summary[1])  # Genotype Quality
    final['DR'] = int(summary[2])  # Depth of reference
    final['DV'] = int(summary[3])  # Depth of variant

    final['SEQ'] = df.iloc[pos]['ALT']

    return final

def make_bilan_sv(name, SAF=False, precise=False, passing=False, duplicates=False, save=True):
    experiments_names = ["P15", "P30", "P50", "P65", "P90"]
    replicates_names = [n for n in range(1, 11)]

    df = pd.DataFrame(columns=["experiment", "replicate", "indel", "start", "end", "length", "vseq"])

    for experiment in experiments_names:
        for replicate in replicates_names:
            sv = vcf_to_df(f"data/vcf/{experiment}/{experiment}-{replicate}.trimed1000.sv_sniffles.vcf")

            l = 0
            while l < sv.shape[0]:
                svinfo = get_svinfo(sv, l)
                if not passing_filters(svinfo['DV'], svinfo['AF'], svinfo['precise'], svinfo["PASS"], svinfo["QUAL"],
                                       filter_SAF=SAF, filter_precise=precise, filter_passing=passing):
                    l+=1
                    continue
                buffer: list[(np.int64, dict)] = [(sv.iloc[l]['POS'], svinfo)]
                if l != sv.shape[0] - 1:
                    while buffer[-1][0] == sv.iloc[l + 1]['POS']:
                        l += 1
                        buffer.append((sv.iloc[l]['POS'], get_svinfo(sv, l)))

                        if l == sv.shape[0] - 1:
                            buffer.append((sv.iloc[l]['POS'], get_svinfo(sv, l)))
                            break

                if len(buffer) > 2:
                    if duplicates:
                        candidate = buffer[0][1]
                        for b in buffer[1:]:
                            if b[1]['DV'] > candidate['DV']:
                                candidate = b[1]
                    else:
                        for b in buffer:
                            df = df._append({
                                "experiment": experiment,
                                "replicate": replicate,
                                "indel": "DEL" if svinfo["SVTYPE"] == "DEL" else "INS",  # Allow to consider DUP as INS
                                "start": b[1]["POS"],
                                "end": b[1]["END"],
                                "length": b[1]["SVLEN"],
                                "af": b[1]["AF"],
                                "dv": b[1]["DV"],
                                "vseq": b[1]["SEQ"]
                            }, ignore_index=True)
                else:
                    candidate = buffer[0][1]
                    svinfo = candidate
                    df = df._append({
                        "experiment": experiment,
                        "replicate": replicate,
                        "indel": "DEL" if svinfo["SVTYPE"] == "DEL" else "INS",
                        "start": svinfo["POS"],
                        "end": svinfo["END"],
                        "length": svinfo["SVLEN"],
                        "af": svinfo["AF"],
                        "dv": svinfo["DV"],
                        "vseq": svinfo["SEQ"]
                    }, ignore_index=True)
                l += 1


    if save:
        df.to_csv(f"data/vcf/{name}.csv", index=False)
    return df

def make_bilan_snp(name, min_qual=30, save=True):
    experiments_names = ["P15", "P30", "P50", "P65", "P90"]
    replicates_names = [n for n in range(1, 11)]

    df = pd.DataFrame(columns=["experiment", "replicate", "pos", "ref", "alt"])

    for experiment in experiments_names:
        for replicate in replicates_names:
            sv = vcf_to_df(f"data/vcf/{experiment}/{experiment}-{replicate}.trimed1000.snp.vcf")
            l = 0
            while l < sv.shape[0]:
                cur_line = sv.iloc[l]
                if passing_filters(0, 0., True, True, cur_line["QUAL"],
                                   filter_SAF=False, filter_precise=False, filter_passing=False,
                                   qual_threshold=min_qual):
                    df = df._append({
                        "experiment": experiment,
                        "replicate": replicate,
                        "pos": cur_line["POS"],
                        "ref": cur_line["REF"],
                        "alt": cur_line["ALT"]
                    }, ignore_index=True)
                l += 1

    if save:
        df.to_csv(f"data/vcf/{name}.csv", index=False)
    return df


# Définir la fonction exponentielle négative
def neg_exp_func(x, a, b, c):
    return a * np.exp(-b * x) + c

def is_above(support, af):
    x_data = [30, 100, 1000]   # Support
    y_data = [1.0, 0.1, 0.05]  # Fréquence Allélique
    popt = curve_fit(neg_exp_func, x_data, y_data, p0=(1, 0.01, 0))[0]
    return af > neg_exp_func(support, *popt)

def passing_filters(support: int, af: float, precise: bool, passing: bool, qual: int,
                    filter_SAF=True, filter_precise=True, filter_passing=True, qual_threshold=30) -> bool:
    if qual < qual_threshold:                        return False  # Quality
    if filter_precise and not precise:               return False  # Precise
    if filter_passing and not passing:               return False  # Passing pipeline filters
    if filter_SAF     and not is_above(support, af): return False  # Support on Allelic Frequency above the threshold

    return True


if __name__ == '__main__':
    # Support on Allelic Frequency = SAF
    make_bilan_sv(SAF=False, precise=False, passing=False, duplicates=True, name='bilanfiltered+notpassing+notprecise+notSAF')
    make_bilan_sv(SAF=True, precise=False, passing=False, duplicates=True, name='bilanfiltered+notpassing+notprecise')
    make_bilan_snp(f"bilanfiltered_snp", min_qual=33)

    """quals = range(0, 60, 5)
    n_var = []
    for qual in tqdm(quals):
        df = make_bilan_snp(f"bilanfiltered_snp_{qual}", min_qual=qual, save=False)
        n_var.append(len(df))

    plt.plot(quals, n_var)
    plt.hlines(y=sum(n_var)/len(n_var), xmin=0, xmax=60)

    between = ()
    l = 0
    while len(between) < 2 and l < len(n_var):
        if n_var[l] > sum(n_var)/len(n_var):
            between += (l+9,l+10)
        l += 1

    plt.xlabel('Qualité minimale')
    plt.ylabel('Nombre de SNPs')

    plt.xticks(range(0, 60, 5))

    plt.xlim(0, 60)
    plt.ylim(min(n_var), max(n_var))

    plt.tight_layout()
    plt.savefig("image/qualite_minimale_snp.png", dpi=300)"""
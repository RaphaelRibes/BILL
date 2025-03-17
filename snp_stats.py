import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

dt = pd.read_csv('data/vcf/bilanfiltered_snp.csv')

sns.catplot(x="replicate", kind="count", data=dt, col="experiment")
plt.show()
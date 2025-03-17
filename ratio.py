import pandas
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize

nano_drop = pd.read_csv('data/BILL - NanoDrop.csv')
qubit = pd.read_csv('data/BILL - Qubit.csv')

ratio = pd.DataFrame()
ratio['Groupe'] = nano_drop['Groupe']

ratio_calc = []

for n, q in zip(nano_drop['ng/µl'], qubit['ng/µl']):
    if q == 0 or n == 0:
        ratio_calc.append(-1000)
    else:
        ratio_calc.append(1-(q/n))

ratio['1-Ratio'] = ratio_calc

ratio_260_280 = []
for r in nano_drop['260/280']:
    if r < 1.8:
        ratio_260_280.append(-1)
    elif r > 2.0:
        ratio_260_280.append(1)
    else:
        ratio_260_280.append(0)

ratio['260/280'] = ratio_260_280

ratio_260_230 = []
for r in nano_drop['260/230']:
    if r < 2.0:
        ratio_260_230.append(-1)
    elif r > 2.2:
        ratio_260_230.append(1)
    else:
        ratio_260_230.append(0)

ratio['260/230'] = ratio_260_230

def heatmap(df: pandas.DataFrame):
    # for 'NaN' values, we will use black color
    # Define custom colormap
    colors = [(0, "red"), (0.5, "green"), (1, "red")]
    custom_cmap = LinearSegmentedColormap.from_list("RedGreenMap", colors)

    # Normalization to fix the range from -1 to 1
    norm = Normalize(vmin=-1, vmax=1)

    # Plot the matrix
    fig, ax = plt.subplots(figsize=(3, len(df["Groupe"])*0.25))
    cax = ax.matshow(df.iloc[:, 1:], cmap=custom_cmap, norm=norm)

    # Set the labels
    ax.set_xticks(range(len(df.columns) - 1))
    ax.set_xticklabels([r"$1-\frac{[\text{Qubit}]}{[\text{NanoDrop}]}$", r"$\frac{A_{260}}{A_{280}}$", r"$\frac{A_{260}}{A_{230}}$"])
    ax.set_yticks(range(len(df["Groupe"])))
    ax.set_yticklabels(df["Groupe"])

    # Rotate the tick labels and set their alignment
    plt.setp(ax.get_xticklabels(), rotation=90, rotation_mode="anchor", ha="left", va="center")

    # Add text annotations for the values
    for i in range(len(df["Groupe"])):
        for j in range(len(df.columns) - 1):
            value = df.iloc[i, j + 1]
            if value == -1000:
                ax.text(j, i, "NaN", ha="center", va="center", color="black", fontsize=8)
            else:
                if j == 0:
                    ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="black", fontsize=6)
                else:
                    if value == 1:
                        ax.text(j, i, "H", ha="center", va="center", color="black", fontsize=8)
                    elif value == -1:
                        ax.text(j, i, "L", ha="center", va="center", color="black", fontsize=8)
                    else:
                        ax.text(j, i, "G", ha="center", va="center", color="black", fontsize=8)

    # Show that L is low, H is high and G is good
    ax.text(5, len(df["Groupe"])//2-1, "H : High", ha="center", va="center", color="red", fontsize=8)
    ax.text(5, len(df["Groupe"])//2+1, "L : Low", ha="center", va="center", color="red", fontsize=8)
    ax.text(5, len(df["Groupe"])//2, "G : Good", ha="center", va="center", color="green", fontsize=8)
    fig.tight_layout()
    plt.show(dpi=300)
    plt.close(fig)

heatmap(ratio)

ratio['1-Ratio'] = ratio['1-Ratio'].replace(-1000, 0)
# barplot
fig, ax = plt.subplots(figsize=(5, 5))
ax.barh(ratio['Groupe'], ratio['1-Ratio'], color='blue')
ax.set_xlabel('Ratio')
ax.set_ylabel('Groupe')
plt.yticks(fontsize=8)
ax.set_title(r'$1-\frac{[\text{Qubit}]}{[\text{NanoDrop}]}$')
# set subtitle "Better is close to 0"
ax.text(0.5, 1.1, "(Better is close to 0)", ha='center', va='center', transform=ax.transAxes)
ax.axvline(0, color='black', lw=0.5)
plt.xlim(-1, 1)
plt.ylim(-0.5, len(ratio['Groupe'])-0.5)
fig.tight_layout()
for n, r in enumerate(ratio['1-Ratio']):
    if r == 0:
        ax.text(0.05, n, 'NaN', ha='center', va='center', color='black', fontsize=5)
plt.show()

ratio['1-Ratio'] = ratio['1-Ratio'].replace(0, -1000)
valid_groups = ["1.1-2", "1.3-4", "2.1-2", "2.3-4", "3.7-8", "4.7-8", "5.1-2",
                "5.3-4", "6.3-4", "6.7-8", "7.3-4", "8.7-8", "9.1-2", "10.5-6"]
valid_groups = ["P90." + group for group in valid_groups]
ratio = ratio[ratio['Groupe'].isin(valid_groups)]
heatmap(ratio)

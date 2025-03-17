import pandas as pd
import matplotlib.pyplot as plt


def boxplot(qubit = pd.read_csv('data/BILL - Qubit.csv'),
            nano_drop=pd.read_csv('data/BILL - NanoDrop.csv')):
    fig, ax = plt.subplots()
    ax.boxplot([qubit['ng/µl'], nano_drop['ng/µl']])
    plt.xticks([1, 2], ['Qubit', 'NanoDrop'])
    plt.ylabel(r'$ng \cdot \mu l^{-1}$')
    plt.tight_layout()
    plt.savefig('image/boxplot.png', dpi=300)

boxplot()

def concentration_on_260280(nano_drop_path='data/BILL - NanoDrop.csv',
                            sequenced_path='data/sequenced.csv'):
    # Load the data
    nano_drop = pd.read_csv(nano_drop_path)
    sequenced = pd.read_csv(sequenced_path)

    # Fix potential mismatches in the 'Groupe' column
    sequenced['Groupe'] = sequenced['Groupe'].str.replace('_', '-')  # Replace underscores with hyphens

    # Filter nano_drop to only include groups present in sequenced
    filtered = []
    for group in sequenced['Groupe']:
        filtered.append(nano_drop[nano_drop['Groupe'] == group]['260/280'].values.tolist()[0])
    plt.scatter(filtered, sequenced['ng/µl'])  # Create a scatter plot
    plt.xlabel(r'$A_{260}/A_{280}$')
    plt.ylabel(r'$ng \cdot \mu l^{-1}$')


    plt.axhline(y=20, color='blue', linestyle='--')  # Add a horizontal line
    plt.axvline(x=1.8, color='red', linestyle='--')  # Add a vertical line

    # plt.text(1.47, 10, 'Low concentration\nproteins contamination', color='red', fontsize=12, ha='center', va='center')
    # plt.text(1.47, 39, 'Good concentration\nproteins contamination', color='orangered', fontsize=12, ha='center', va='center')
    # plt.text(1.85, 39, 'Good concentration\nno proteins contamination', color='green', fontsize=12, ha='center', va='center', rotation=270)

    for i, txt in enumerate(sequenced['Groupe']):
        plt.annotate(txt, (filtered[i], sequenced['ng/µl'][i]), textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8)
    plt.ylim(0, 60)
    plt.xlim(1.1,1.9)
    plt.xticks([1.8], ['1.8'])
    plt.tight_layout()
    plt.savefig('image/concentration_on_260280.png', dpi=300)

concentration_on_260280()
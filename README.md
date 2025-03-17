# Bioinformatic Leanring Lab (BILL) - 2025
## Introduction
This is a display repository for the Bioinformatic Leanring Lab (BILL) in 2025.

## Members
- Raphaël Ribes - Data analyses, code, writing
- Farah Moukachar - Writing, code
- Lucas Boulet - Writing, bibliography
- Aicha El Jai - Writing

You can check the commit history in `log.csv`.

## Requirements

- Python 3.12 or higher
    ```bash
    sudo apt-get install python3.12-full
    ```
  
## Installation
0. If you don't have git installed, you can install it using:
    ```bash
    sudo apt-get install git
    ```
   
1. **Clone the Repository**:

   ```bash
   git clone https://github.com/RaphaelRibes/samReader.git
   ```

2. **Navigate to the Directory**:

   ```bash
   cd samReader
   ```

3. **Install the Virtual Environment (venv)**:

   ```bash
   python3.12 -m venv .venv
   ```

4. **Activate the venv**:

   ```bash
   source .venv/bin/activate
   ```
   Or for Windows:
   ```bash
    .venv\Scripts\activate
    ```
   
5. **Install Requirements**:

   ```bash
    pip install -r requirements.txt
    ```
   
## What is what

- Figure 1: [concentrations.py](concentrations.py) with the boxplot() function
- Figure 2: [concentrations.py](concentrations.py) with the concentration_on_260280() function
- Figure 4: [reads_boxplot.py](reads_boxplot.py)
- Figure 5: commented part at the end of [vcf.py](vcf.py)
- Figure 6: [sv_visualization.py](sv_visualization.py) with the make_flat_plots() function
- Figure 7: [barplot_vcf.py](barplot_vcf.py) with the print_sv() function
- Figure 8 and 9: [PCA.py](PCA.py)
- Figure 10 to 13: [vizualise_var.py](vizualise_var.py) with the plot_mutations() function from the AutoVar class

If you need the gitlab working gitlab repo, ask me by mail: raphael.ribes@etu.umontpellier.fr
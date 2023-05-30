import os
import gepia
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn3
from IPython.display import IFrame
from scipy import stats
import sqlite3
import numpy as np

def get_file_path():
    return os.getcwd()

def plot_venn_diagram(data_sets, labels):
    venn3(data_sets, labels)
    plt.show()

def find_overlapping_genes(data_sets):
    overlapping_genes = []
    for i in range(len(data_sets) - 1):
        for j in range(i + 1, len(data_sets)):
            genes = list(set(data_sets[i]) & set(data_sets[j]))
            overlapping_genes.append(genes)
    return overlapping_genes

def boxplot_analysis(gepia_obj, params, signature):
    gepia_obj.setParams(params)
    gepia_obj.setParam('signature', signature)
    gepia_obj.query()


def survival_analysis(gepia_obj, params, signature):
    gepia_obj.setParams(params)
    gepia_obj.setParam('signature', signature)
    gepia_obj.query()


def plot_histogram(data, title):
    t = data['Tumor']
    n = data['Normal']
    plt.hist([t, n], density=True, bins=7)  # Blue Tumor, Orange Normal
    plt.ylabel('Frequency')
    plt.xlabel('Log2(TPM)')
    plt.title(title)
    xt = plt.xticks()[0]
    xmin, xmax = min(xt), max(xt)
    lnspc = np.linspace(xmin, xmax, len(xt))
    m, s = stats.norm.fit(xt)  # get mean and standard deviation
    pdf_g = stats.norm.pdf(lnspc, m, s)  # theoretical values
    plt.plot(lnspc, pdf_g)

def main():
    # Set the working directory
    mydir = get_file_path()

    # Define the signatures
    sigs_df = pd.read_csv(os.path.join(mydir, 'signatures.csv'))
    sigs = [sigs_df[column].dropna().tolist() for column in sigs_df]
    # Convert the signatures into sets
    signatures = [set(genes) for genes in sigs]

    # Plot Venn diagram
    labels = ('MammaPrint', 'OncotypeDx', 'GGI')
    plot_venn_diagram(signatures, labels)

    # Find overlapping genes
    overlapping_genes = find_overlapping_genes(sigs)
    for i, genes in enumerate(overlapping_genes):
        if i+1 < len(labels):
            print(f"The genes that intersect at {labels[i]} and {labels[i+1]} are: {genes}")

    # Initialize GEPIA objects
    bp = gepia.boxplot()
    sur = gepia.survival()

    # Read expression data from CSV files
    OncoExp = pd.read_csv(os.path.join(mydir, 'OncoExp.csv'))
    MammaExp = pd.read_csv(os.path.join(mydir, 'MammaExp.csv'))
    GGIExp = pd.read_csv(os.path.join(mydir, 'GGIExp.csv'))

    # Pulling up and down regulated genes
    OncoUp = OncoExp[(OncoExp['Diff']>0)]
    OncoDown = OncoExp[(OncoExp['Diff']<0)]

    MammaUp = MammaExp[(MammaExp['Diff']>0)]
    MammaDown = MammaExp[(MammaExp['Diff']<0)]

    GGIUp = GGIExp[(GGIExp['Diff']>0)]
    GGIDown = GGIExp[(GGIExp['Diff']<0)]

    # Convert the signatures into sets
    up_regulated = [set(df['gene']) for df in [OncoUp, MammaUp, GGIUp]]
    down_regulated = [set(df['gene']) for df in [OncoDown, MammaDown, GGIDown]]

    # Plot Venn diagrams
    labels = ('MammaPrint', 'OncotypeDx', 'GGI')
    plot_venn_diagram(up_regulated, labels)
    plot_venn_diagram(down_regulated, labels)

    # Perform boxplot and survival analysis for up-regulated genes
    for up_genes in up_regulated:
        boxplot_analysis(bp, bp.params, list(up_genes))
        survival_analysis(sur, sur_params, list(up_genes))

    # Perform boxplot and survival analysis for down-regulated genes
    for down_genes in down_regulated:
        boxplot_analysis(bp, bp.params, list(down_genes))
        survival_analysis(sur, sur_params, list(down_genes))
        
        # Perform boxplot and survival analysis for up-regulated genes for each signature
    for i, up_genes in enumerate(up_regulated):
        signature_genes = signatures[i] & up_genes
        print(f"Boxplot and survival analysis for up-regulated genes in {labels[i]} signature:")
        boxplot_analysis(bp, bp.params, list(signature_genes))
        survival_analysis(sur, sur.params, list(signature_genes))

    # Perform boxplot and survival analysis for down-regulated genes for each signature
    for i, down_genes in enumerate(down_regulated):
        signature_genes = signatures[i] & down_genes
        print(f"Boxplot and survival analysis for down-regulated genes in {labels[i]} signature:")
        boxplot_analysis(bp, bp.params, list(signature_genes))
        survival_analysis(sur, sur.params, list(signature_genes))


if __name__ == '__main__':
    main()
    

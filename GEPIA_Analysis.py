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

mydir= '/Users/NM/Documents/University/GMU/Spring2021/Baranova/'
sur= gepia.survival()

bp = gepia.boxplot()
bp_params = bp.params
bp.setOutDir(mydir)
sur_params = sur.params
sim = gepia.similar()

cor = gepia.correlation()
sur = gepia.survival()
sur_map = gepia.survival_map()

sur.setParam('dataset', ['BRCA'])
sur_params = sur.params
sur.setOutDir(mydir)

sur_map.setParam('dataset', ['BRCA'])
sur_mapparams = sur_map.params
sur_map.setOutDir(mydir)

outdir = mydir


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

    # Perform boxplot analysis for up-regulated genes
    for signature in sigs:
        boxplot_analysis(bp, bp_params, signature)

    # Perform survival analysis for up-regulated genes
    for signature in sigs:
        survival_analysis(sur, sur_params, signature)

    # Read data from GEPIA2 for log2TPM values
    gene_lists_file = os.path.join(get_file_path(), '.gene_lists.csv')

    # Load data into a pandas DataFrame
    data = pd.read_csv(gene_lists_file)

    # Plot histogram for log2(TPM) values
    plot_histogram(data, 'Log2(TPM) Distribution')

    # Save the histogram as a PDF file
    plt.savefig(os.path.join(get_file_path(), 'log2tpm_histogram.pdf'))

    # Display the saved PDF file
    IFrame(os.path.join(get_file_path(), 'log2tpm_histogram.pdf'), width=400, height=300)


if __name__ == '__main__':
    main()

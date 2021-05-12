import gzip
import pandas as pd


class Gene:
    def __init__(self, name, label):
        self.name = name
        self.label = label

    def markGene(self):
        self.label = 1

    def unmarkGene(self):
        self.label = 0


# Create dictionary to map all gene alias to official gene name
def createHumanGeneAliasDict():
    humanREF = 'D:/Storage/Lab/sleep/sleep_gene/workDIR/feature_selection_batch/REFERENCE/Homo_sapiens.gene_info.gz'
    geneDict = {}
    with gzip.open(humanREF, "rb") as f:
        next(f)
        for line in f:
            line = str(line, 'utf-8')
            info = line.strip().split("\t")
            geneSymbol = info[2]
            geneDict[geneSymbol] = geneSymbol
            geneAlias = info[4].split("|")
            ## if an alias mapped to multiple official gene symbol, the one with smallest entrez id is chosen
            if geneAlias != ['-']:
                for gene in geneAlias:
                    try:
                        geneDict[gene]  ## test if gene already in dictionary
                    except KeyError:
                        geneDict[gene] = geneSymbol
    return geneDict


def readSleepGeneList():
    sleep_gene_file = 'D:/Storage/Lab/sleep/sleep_gene/data/sleep_regulating_gene.xlsx'
    sleepTable = pd.read_excel(sleep_gene_file,
                               sheet_name='Full_SRGdata',
                               usecols=['Gene', 'Tag', 'Tier', 'HGNC', 'mouse', 'Phenotype_curated', 'Model_for_genename'])
    return sleepTable

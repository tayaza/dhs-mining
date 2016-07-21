#usr/bin/python
import sys, sqlite3, csv




def prep_gene_list(gene_file):
    gene_list = []
    with open(gene_file, 'rb') as infile:
        reader = csv.reader(infile, delimiter = " ")
        for row in reader:
            gene = row[0].split()[0]
            if gene not in gene_list:
                gene_list.append(gene)
    return gene_list
    




gene_file = 'MWC/obesity/obesity_genes.txt'
genes = prep_gene_list(gene_file)
print genes

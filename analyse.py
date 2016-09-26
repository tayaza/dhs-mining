#! /usr/bin/env python

import sys
import csv
import os
import argparse


def set_filepath(input_fp):
    output_fp = input_fp
    if ('/' in output_fp):
        while not (output_fp.endswith('/')):
            output_fp = output_fp[:len(output_fp)-1]
        output_fp = output_fp[:len(output_fp)-1]
    else:
        output_fp = ''
    if not os.path.isdir(output_fp):
        os.mkdir(output_fp)
    return output_fp
                        


    """
    with open(gwas_dir+'/snps.txt', 'wb') as all_snps:
        all_writer = csv.writer(all_snps)
        for trait in snp_list:
            for snps in snp_list[trait]:
                with open(output_fp + '/' + trait + '_snps_nt.txt', 'a') as trait_snps:
                    trait_writer = csv.writer(trait_snps)
                    all_writer.writerow([snps])
                    trait_writer.writerow([snps])
    print output_fp + '/' + trait + 'snps_.txt'
    """

def find_genes(gene_file, gwas_file):
    """ See if CODE3D genes are mapped in GWAS records. """
    gene_list = []
    output = []
    if os.path.isfile(gene_file):
        with open(gene_file, 'r') as infile:
            for line in infile:
                if line.startswith('rs'):
                    line = line.strip().split('\t')
                    e_snp = line[0]
                    e_gene = line[4]
                    gene_list.append([e_snp, e_gene])
        infile.close()
    #test_genes = ['NOTCH2', 'PPARG', 'MADD', 'FTO', 'ACP2']
    if os.path.isfile(gwas_file):
        gene_counter = 0
        mapped_counter = 0
        for pair in gene_list:
            snp = pair[0]
            gene = pair[1]
            with open(gwas_file, 'rb') as gwas:
                mapped = False
                for row in gwas:
                    row = row.strip().split('\t')
                    gwas_snp = row[21]
                    reported_genes = []
                    mapped_genes = []
                    for r_gene in row[13].split():
                        r_gene = r_gene.replace(',', '')
                        reported_genes.append(r_gene)
                    if '-' in reported_genes:
                        reported_genes.remove('-')
                    for m_gene in row[14].split():
                        m_gene = m_gene.replace(',', '')
                        mapped_genes.append(m_gene)
                    if '-' in mapped_genes:
                        mapped_genes.remove('-')
                    author = row[2]
                    pubmed_id = row[1]
                    if gene in reported_genes:
                        result = snp, gene, gwas_snp, reported_genes, \
                                   mapped_genes, author, pubmed_id
                        output.append(result)
                        if mapped == False:
                            mapped_counter += 1
                            mapped = True
            gwas.close()
            gene_counter += 1
        with open(filepath + '/mapped.txt', 'wb') as mapped:
            writer = csv.writer(mapped, delimiter = '\t')
            writer.writerow(['SNP', 'GENE', 'GWAS_SNP', 'REPORTED_GENES', \
                                'MAPPED_GENES', 'AUTHOR', 'PUBMED_ID'])
            writer.writerows(output)
        print mapped_counter,  'out of',  gene_counter,  'genes have been previously', \
                   'mapped in GWAS'
        print 'Find results at '+filepath
def eqtls_by_tissues(eqtls_file):
    tissue_pool = {}
    with open(eqtls_file, 'rb') as infile:
        eqtls = csv.reader(infile, delimiter = '\t')
        next(eqtls, None)
        for line in eqtls:   #TODO remove header
            k = []
            epair = tuple(line)
            k.append(epair)
            tissue =  line[7]
            if tissue not in tissue_pool.keys():
                tissue_pool[tissue] = k
            else:
                rec = tissue_pool[tissue]
                rec.append(epair)
                tissue_pool[tissue] = rec
       
    with open(filepath + '/eqtls_by_tissues.txt', 'wb') as mapped:
        writer = csv.writer(mapped, delimiter = '\t')
        writer.writerow(['TISSUE', 'SNP-GENES#'])
        for tissue in sorted(tissue_pool):
            output = [tissue, len(tissue_pool[tissue])]
            writer.writerow(output)
        #TODO Finer details of SNPs and genes

def genes_by_snps(matched_file):
    gene_list = {}
    gene_details = {}
    with open(matched_file, 'rb') as mfile:
        reader = csv.reader(mfile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            snps = []
            details = []
            snp = line[0]
            gene = line[4]
            snps.append(snp)
            details.append(line)
            if gene not in gene_list:
                gene_list[gene] = snps
                gene_details[gene] = details
            else:
                k = gene_list[gene]
                k.append(snp)
                gene_list[gene] = k
                d = gene_details[gene]
                d.append(details)
    with open(filepath + '/genes_by_snps_details.txt', 'wb') as outfile:
        writer = csv.writer(outfile, delimiter = '\t')
        writer.writerow(['GENE', 'SNP', 'INTERACTION', \
                                     'INTERACTION_TYPE', 'HiC_CELLTYPES', \
                                     'GTEx_TISSUE_REGULATED'])
        for gene in gene_details:
            if len(gene_details[gene]) > 1:
                trans = []
                for row in gene_details[gene]:
                    line = []
                    cis = 'Intra'
                    if len(row) == 15:
                        line = row
                    elif len(row) == 1:
                        line = row[0]
                    if line[1] != line[5]:
                        cis = 'Inter'
                    if cis == 'Inter' or line[6] == 'Trans':
                        trans.append(gene)
                    to_file = line[4], line[0], line[6], cis, line[7], line[13]
                    writer.writerow(to_file)

    with open(filepath + '/genes_by_snps.txt', 'wb') as outfile:
        writer = csv.writer(outfile, delimiter = '\t')
        writer.writerow(['GENE', 'SNPS#', 'SNPs'])

        for gene in gene_list:
            snp_number = len(gene_list[gene])
            if snp_number > 1:
                output = [gene, snp_number, gene_list[gene]]
                writer.writerow(output)

def snps_by_genes(matched_file):
    snp_list = {}
    snp_details = {}
    with open(matched_file, 'rb') as mfile:
        reader = csv.reader(mfile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            genes = []
            snp = line[0]
            gene = line[4]
            genes.append(gene)
            if snp not in snp_list:
                snp_list[snp] = genes
            else:
                k = snp_list[snp]
                k.append(gene)
                snp_list[snp] = k
                #print snp_list[snp]
    with open(filepath + '/snps_by_genes.txt', 'wb') as outfile:
        writer = csv.writer(outfile, delimiter = '\t')
        writer.writerow(['SNP', 'GENES#', 'GENES'])
        for snp in snp_list:
            gene_number = len(snp_list[snp])
            if gene_number > 1:
                output = [snp, gene_number, snp_list[snp]]
                writer.writerow(output)

def sort_trans(matched_file):
    trans = {}
    trans_list = []
    exist = []
    multiple = []
    with open(matched_file, 'rb') as tfile:
        reader = csv.reader(tfile, delimiter ='\t')
        next(reader, None)
        for line in reader:
            gene = line[4]
            interaction = line[6]
            if interaction == 'Trans':
                to_trans = []
                to_trans.append((line[0], line[1], line[2], gene, line[5], line[6], \
                                line[13]))
                if gene not in trans.keys():
                    trans[gene] = to_trans
                else:
                    update_trans = [] 
                    for row in trans[gene]:
                        update_trans.append(row)
                    update_trans.append(to_trans[0])
                    trans[gene] = update_trans
    with open(filepath + '/trans_genes.txt', 'wb') as outfile:
        writer = csv.writer(outfile, delimiter = '\t')
        writer.writerow(['SNP', 'SNP_CHR', 'SNP_DHS_ID', 'GENES', 'GENE_CHR', \
                             'INTERACTION', 'eQTL_TISSUE'])
        
        for gene in trans.keys():
            record = trans[gene]
            for row in record:
                writer.writerow(row)


if __name__ == "__main__":
    #parser = argparse.ArgumentParser(description = "")
    #parser.add_argument("-i", "--input", required = True, help = "Filepath of folder containing downloaded GWAS association from GWAS Catalogue.")
    #args = parser.parse_args()
    #gwas_dir = args.input

    CODE3D_file = '/mnt/3dgenome/projects/tfad334/hvn/obesity/codes3d_output/analysis/match.txt'
    gwas_file = '/mnt/3dgenome/projects/tfad334/hvn/obesity/data/gwas-association-downloaded_2016-07-13-obesity.tsv'
    sig_eqtls = '/mnt/3dgenome/projects/tfad334/hvn/obesity/codes3d_output/analysis/sig_snp-gene_eqtls.txt'

    matched = '/mnt/3dgenome/projects/tfad334/hvn/obesity/codes3d_output/analysis/match.txt'

    
    filepath = set_filepath(sig_eqtls)    
    find_genes(CODE3D_file, gwas_file)    
    eqtls_by_tissues(sig_eqtls)
    genes_by_snps(matched)
    snps_by_genes(matched)
    sort_trans(matched)

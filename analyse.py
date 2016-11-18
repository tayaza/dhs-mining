#! /usr/bin/env python

import sys
import csv
import os
import argparse


def set_filepath(input_fp):
    output_fp = input_fp
    if output_fp.endswith('/'):
        output_fp = output_fp[:len(output_fp)-1]
    return output_fp

def find_genes(gene_file, gwas_file):
    """ See if CoDeS3D genes are mapped in GWAS records. """
    gene_list = {}
    gwas_associations = []
    output = []
    all_mapped_genes = []
    all_reported_genes = []
    gwas_mapped = []
    gwas_unmapped = []
    if os.path.isfile(gene_file):
        with open(gene_file, 'r') as infile:
            for line in infile:
                if line.startswith('rs'):
                    line = line.strip().split('\t')
                    snp = line[0]
                    gene = line[4]
                    snp_list = []
                    snp_list.append(snp)
                    if gene not in gene_list.keys():
                        gene_list[gene] = snp_list
                    else:
                        for s in gene_list[gene]:
                            snp_list.append(s)
                        gene_list[gene] = snp_list
        infile.close()
    if os.path.isfile(gwas_file):
        gene_counter = 0
        mapped_counter = 0
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
                for g in reported_genes:
                    if g not in all_reported_genes:
                        all_reported_genes.append(g)
                for m in mapped_genes:
                    if m not in all_mapped_genes:
                        all_mapped_genes.append(m)
                row[13] = reported_genes
                row[14] = mapped_genes
                gwas_associations.append(row)
        gwas.close()
    in_line = []
    for gene in gene_list.keys():
        if gene in all_mapped_genes or gene in all_reported_genes:
            output.append(gene)
            author_snp = []
            mapped_counter += 1
            for line in gwas_associations:
                if gene in line[13] or gene in line[14]:
                    same_snp = []
                    reported_genes = line[13]
                    mapped_genes = line[14]
                    snp = line[21]
                    author = line[2] + ' ' + line[3][:4]
                    pubmed_id = line[1]
                    tester = author+snp
                    if tester not in author_snp:
                        author_snp.append(tester)
                        for eqtl_snp in gene_list[gene]:
                            if eqtl_snp in snp:
                                same_snp.append('True')
                                gline = gene, eqtl_snp
                                if gline not in in_line:
                                    in_line.append(gline)
                            else:
                                same_snp.append('False')
                        to_gwas_mapped = gene, gene_list[gene], reported_genes, \
                            mapped_genes, snp, author, pubmed_id, same_snp
                        gwas_mapped.append(to_gwas_mapped)
        else:
            to_unmapped = gene, gene_list[gene]
            gwas_unmapped.append(to_unmapped)

    with open(filepath + '/gwas_mapped_genes.txt', 'wb') as mapped:
        writer = csv.writer(mapped, delimiter = '\t')
        writer.writerow(['GENE', 'eQTL_SNPS', 'REPORTED_GENES', 'MAPPED_GENES', \
                             'GWAS_SNP', 'AUTHOR', 'PUBMED_ID', 'SAME_SNP?'])
        writer.writerows(gwas_mapped)
    mapped.close()

    with open(filepath + '/gwas_unmapped_genes.txt', 'wb') as unmapped:
        writer = csv.writer(unmapped, delimiter = '\t')
        writer.writerow(['GENE', 'eQTL_SNPS'])
        writer.writerows(gwas_unmapped)
    unmapped.close()


    with open(filepath + '/gwas_in_line.txt', 'wb') as inline:
        writer = csv.writer(inline, delimiter = '\t')
        writer.writerow(['GENE', 'eQTL_SNPS'])
        writer.writerows(in_line)
    inline.close()

    with open(filepath + '/mapped_gwas_genes.txt', 'wb') as mapped:
        for gene in all_mapped_genes:
            mapped.write(gene + '\n')
    mapped.close()

    with open(filepath + '/reported_gwas_genes.txt', 'wb') as mapped:
        for gene in all_reported_genes:
            mapped.write(gene + '\n')
    mapped.close()

   # print '\t',  mapped_counter,  'out of',  len(gene_list.keys()), \
   #     'genes have been previously mapped or reported in GWAS'

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


def spatial_pairs(summary):
    pair_list = []
    check_pairs = []
    with open(summary, 'rb') as sfile:
        reader = csv.reader(sfile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            snp = line[0]
            gene = line[3]
            test = snp+gene
            if test not in check_pairs:
                check_pairs.append(test)
                pair_list.append((snp, gene))
                
    print 'Total spatial SNP-gene pairs : \t', len(pair_list)
    return pair_list


def gwas_significant(gwas_bed, sig_eqtls, match, spatials):
    gwas_snps_08 = {}
    sig_eqtls_08 = []
    match_08 = []
    genes_08 = []
    spatial_08 = []
    snps_08 = []
    with open(gwas_bed, 'rb') as bed_file:
        reader = csv.reader(bed_file, delimiter = '\t')
        next(reader, None)
        for line in reader:
            snp = line[3]
            pval = float(line[4])
            if pval < 5.0e-08:
                if snp not in gwas_snps_08.keys():
                    gwas_snps_08[snp] = pval
                else:
                    if pval <= gwas_snps_08[snp]:  # use SNP's lowest pval
                        gwas_snps_08[snp] = pval
    bed_file.close()


    for snp in gwas_snps_08.keys():
        with open(sig_eqtls, 'rb') as eqtl_file:
            reader = csv.reader(eqtl_file, delimiter = '\t')
            next(reader, None)
            for line in reader:
                eqtl_snp = line[0]
                if eqtl_snp == snp:
                    sig_eqtls_08.append(line)

        with open(match, 'rb') as match_file:
            reader = csv.reader(match_file, delimiter = '\t')
            next(reader, None)
            for line in reader:
                eqtl_snp = line[0]
                gene = line[4]
                if eqtl_snp == snp:
                    match_08.append(line)
                    if gene not in genes_08:
                        genes_08.append(gene)
                    if eqtl_snp not in snps_08:
                        snps_08.append(eqtl_snp)

        for pair in spatials:
            eqtl_snp = pair[0]
            if eqtl_snp == snp:
                spatial_08.append(pair)



    print 'GWAS significant spatial SNP-gene pairs : \t', len(spatial_08)
    print 'Number of GWAS significant SNPs : \t', len(gwas_snps_08.keys())
    print 'Number of GWAS significant SNP eQTL : \t', len(snps_08)
    print 'Number of SNP-gene eQTLs : \t', len(sig_eqtls_08) 
    print 'Number of unique significant SNP-gene eQTLs : \t', len(match_08)
    print 'Number of significant eGenes : \t', len(genes_08)
    dir_08 = filepath + '/gwas_08'
    if not os.path.isdir(dir_08):
        os.mkdir(dir_08)
    with open(dir_08 + '/sig_eqtls_08.txt', 'wb') as sig_file:
        writer = csv.writer(sig_file, delimiter = '\t')
        file_header = ('SNP', 'SNP_Chromosome', 'SNP-Locus', 'Gene_Name', \
                           'Gene_Chromosome', 'Gene_Start', 'Gene_End', \
                           'Tissue', 'p-value', 'q-value', 'Cell_Lines', \
                           'GTEx_cis_p_Threshold', 'cis_SNP_Gene_Interaction' \
                           'SNP-Gene_Distance', 'Expression_Level_in_eQTL_Tissue'\
                           'Max_Expressed_Tissue', 'Maximum_Expression_Level', \
                           'Min_Expressed Tissue', 'Min_Expression_Level')
        writer.writerow(file_header)
        writer.writerows(sig_eqtls_08)
    sig_file.close()
    with open(dir_08 + '/snp-gene_pairs_08.txt', 'wb') as unique_file:
        writer = csv.writer(unique_file, delimiter = '\t')
        file_header = ('SNP', 'SNP_CHR', 'SNP_DHS_ID', 'SNP_OPEN_CELLTYPES', \
                           'GENE', 'GENE_CHR', 'CIS', 'HIC_CELLTYPE' 'HIC_DHS%', \
                           'HIC_TISSUE', 'MAX_DHS_CELL', 'MAX_DHS_CELL%', \
                           'MAX_DHS_TISSUE', 'eQTL_TISSUE', 'GTEx_MAX_TISSUE')

        writer.writerow(file_header)
        writer.writerows(match_08)
    unique_file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-d", "--dir", required = True, \
                            help = "Directory with dhsquery results.")
    parser.add_argument("-g", "--gwas", required = False, \
                            help = "Downloaded GWAS association from GWAS \
                            Catalogue.")
    args = parser.parse_args()
    gwas_file = args.gwas
    filepath = set_filepath(args.dir)    
    sig_eqtls = filepath + '/sig_SNP-gene_eqtls.txt'
    match = filepath + '/match.txt'
    gwas_bed = filepath + '/gwas_snps.bed'
    summary = filepath + '/summary.txt'
    if gwas_file:
        find_genes(match, gwas_file)
    else:
        print 'Warning: \t \
            No GWAS Catalog associations file for the trait is specified or \
            found!'
    eqtls_by_tissues(sig_eqtls)
    genes_by_snps(match)
    snps_by_genes(match)
    sort_trans(match)
    spatials = spatial_pairs(summary)
    gwas_significant(gwas_bed, sig_eqtls, match, spatials) 

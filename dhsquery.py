#!usr/bin/env python
import sys
import sqlite3
import csv
import os
import argparse

    

def sig_SNP_genes(CODE3D_summary):
    """
    Extract from CoDeS3D summary only SNP-gene eQTLs with FDR less than or equal 
    to 0.05. 
    """
    snp_eqtls = []
    if os.path.isfile(CODE3D_summary):
        print ('\t Finding significant SNP-gene eQTLs...')
        with open(CODE3D_summary, 'rb') as infile:
            for line in infile:
                line = line.strip().split()
                if len(line) == 19:
                    qvalue = line[9]
                    if line[0].startswith('rs') and float(qvalue) <= 0.05:
                        snp_eqtls.append(line)
                        #print '\t\t', line[0], line[3], qvalue, line[7]
        infile.close()
    with open(filepath + '/sig_SNP-gene_eqtls.txt', 'wb') as sig:
        writer = csv.writer(sig, delimiter = '\t')
        writer.writerow(['SNP', 'SNP_Chromosome', 'SNP_Locus', 'Gene_Name', 
                         'Gene_Chromosome', 'Gene_Start', 'Gene_End', 'Tissue',
                         'p-value', 'q-value', 'Cell_Lines', 'GTEx_cis_p_Threshold',
                         'cis_SNP_Gene_Interaction', 'SNP-gene_Distance', 
                         'Expression_Level_in_eQTL_Tissue', 'Max_Expressed_Tissue',
                         'Maximum_Expression_Level', 'Min_Expressed_Tissue', 
                         'Min_Expression_Level'])
        writer.writerows(snp_eqtls)
    sig.close()
    return snp_eqtls



def set_filepath(input_fp):
    """ Set output directory filepath."""
    filepath = input_fp
    if ('/' in filepath):
        while not (filepath.endswith('/')):
            filepath = filepath[:len(filepath)-1]
        filepath = filepath[:len(filepath)-1]
    else:
        filepath = ''
    output_fp = filepath + '/dhs_results'
    if not os.path.isdir(output_fp):
        os.mkdir(output_fp)
    return output_fp



def get_openCellTypes(cluster_id, dhsDB):
    conn = sqlite3.connect(dhsDB)
    conn.text_factory = str
    cur = conn.cursor()
    cellTypes = []
    
    output = [cluster_id, 'NA', 'NA', 'NA']
    if isinstance(cluster_id, int):
        cur.execute("SELECT * FROM openCellTypes WHERE cluster_id = ?;", \
                        (cluster_id,))
        data = cur.fetchone()
        if isinstance(data, tuple):
            cellTypes = data
        else:
            cellTypes = output
    else:
        cellTypes = output
 
    return cellTypes
    cur.close()
    conn.close()
 

def get_openSamples(cluster_id, dhsDB):
    conn = sqlite3.connect(dhsDB)
    conn.text_factory =str
    cur = conn.cursor()
    openSamples = []

    output = [cluster_id, 'NA', 'NA', 'NA','NA', 'NA', 'NA', 'NA','NA', \
                  'NA', 'NA', 'NA']
    if isinstance(cluster_id, int):
        cur.execute("SELECT * FROM openSamples WHERE cluster_id = ?;", \
                        (cluster_id,))
        data = cur.fetchone()
        if isinstance(data, tuple):
            openSamples = data
        else:
            openSamples = output
    else:
        openSamples = output
 
    return openSamples
    cur.close()
    conn.close()
 
def get_overlaps(cluster_id, dhsDB):
    conn = sqlite3.connect(dhsDB)
    conn.text_factory =str
    cur = conn.cursor()
    overlaps = []

    output = [cluster_id, 'NA', 'NA', 'NA', 'NA']
    if (cluster_id != 'NA'):
        cur.execute("SELECT * FROM overlap WHERE cluster_id = ?;", (cluster_id,))
        data = cur.fetchone()
        if isinstance(data, tuple):
            overlaps=data
        else:
            overlaps = output
    else:
        overlaps = output

    return overlaps

    cur.close()
    conn.close()


def get_motifs(cluster_id, dhsDB):
    conn = sqlite3.connect(dhsDB)
    conn.text_factory = str
    cur = conn.cursor()
    motifs = []

    # TODO: Modify to retrieve all significant (e < 10-6) TFs?
    if isinstance(cluster_id, int):
        cur.execute("SELECT * FROM motifJaspar WHERE cluster_id = ?;", \
                        (cluster_id,))
        data = cur.fetchall()
        tf_list = []
        if len(data) is 0:
            tf_list = 'NA'
        else:
            for row in data:
                tf_list.append(row[4] + ':' + str(row[5]))
        motifs = [cluster_id,tf_list]
    else:
        motifs = [cluster_id, 'NA']
    return motifs

    cur.close()
    conn.close()

def get_sampleDHS_signal(dhs_id, dhsDB):
    conn = sqlite3.connect(dhsDB)
    conn.text_factory = str
    cur = conn.cursor()
    sig_Samples = []
    head = []

    if not isinstance(dhs_id, int):
        sig_Samples = [dhs_id, ['NA']]
    else:
        cur.execute("SELECT * FROM dhs112 WHERE rowid =?;", (dhs_id+1,))
        dhs_signal = cur.fetchone()
        head = [description[0] for description in cur.description]
        sig_dhs = []
        for i in xrange(3, len(dhs_signal)):
            if (dhs_signal[i] >= 0.1):
                sig_dhs.append(head[i] + ':' + str(dhs_signal[i]))
            i = i + 1
        if len(sig_dhs) is 0:
            sig_Samples = [dhs_id, ['NA']]
        else:
            sig_Samples = [dhs_id, sig_dhs]
    return sig_Samples

def get_snpDHS(snp, snp_chr, snp_pos, dhsDB):
    conn = sqlite3.connect(dhsDB)
    conn.text_factory =str
    cur = conn.cursor()
    dhs_data = []
    snp_chr = "chr" + snp_chr
    cur.execute("SELECT rowid, chr, start, end FROM dhs112 WHERE chr = ? \
                 AND start <=? AND end >=?;", (snp_chr, snp_pos, snp_pos))
    data = cur.fetchone()
    if (data == None):
        dhs_data = [snp, snp_chr, snp_pos, 'NA', 'NA', 'NA']
    else:
        dhs_data = [snp, snp_chr, snp_pos, data[0]-1, data[2], data[3]]
    dhs_id = dhs_data[3]
    if isinstance(dhs_id, int):
        cur.execute("SELECT chr, start, end, refined_cluster FROM dhsCluster LIMIT 1 \
                    OFFSET "+ str(dhs_id - 1) + ";")
        data = cur.fetchone()
        dhs_data.append(data[3])
    else:
        dhs_data.append('NA')
    cluster_id = dhs_data[6]
    open_Samples = get_openSamples(cluster_id, dhsDB)
    overlaps = get_overlaps(cluster_id, dhsDB)
    motifs = get_motifs(cluster_id, dhsDB)
    sig_dhsSamples = get_sampleDHS_signal(dhs_id, dhsDB)
    
    snp_dhs = (snp, snp_chr, snp_pos, dhs_id, dhs_data[4], dhs_data[5], sig_dhsSamples[1], 
               cluster_id, overlaps[1], overlaps[2], overlaps[3], motifs[1], 
               open_Samples[1], open_Samples[6], open_Samples[7], open_Samples[8], 
               open_Samples[9], open_Samples[10], open_Samples[11]
               )

    return snp_dhs
    
    cur.close()
    conn.close()

def get_geneDHS(gene, gene_start, gene_end, dhsDB):
    conn = sqlite3.connect(dhsDB)
    conn.text_factory = str
    cur = conn.cursor()
    gene_dhs = []
    dhs_data = []
    open_Samples = []
    overlaps = []
    motifs = []
    sig_dhsSamples = []
    cur.execute("SELECT gene_name, dhs_id, dhs_start, dhs_end, pval from \
                    gene_correlations_p05 WHERE pval >= 0.95 AND gene_name = ?;",\
                    (gene,))
    cor_query = cur.fetchall()
    if not cor_query:
        cur.execute("SELECT rowid, chr, start, end, refined_cluster FROM \
                        dhsCluster WHERE start >= ? AND end <= ?;", \
                        (int(gene_start), int(gene_end)))
        query = cur.fetchall()
        for row in query:
            dhs_data.append((gene, int(row[0])-1, row[2], row[3], 'NA', row[4]))
        if not query:
            dhs_data.append((gene, 'NA', 'NA', 'NA', 'NA', 'NA'))
    else:
        for item in cor_query:
            dhs_id = item[1]
            dhs_start = int(item[2])
            dhs_end = int(item[3])
            pval = 1 - item[4]
            cluster_id = 0
            cur.execute("SELECT chr, start, end, refined_cluster FROM dhsCluster LIMIT 1 \
                        OFFSET "+ str(dhs_id-1)+";")
            cluster_query = cur.fetchone()
            if cluster_query is None:
                dhs_data.append((gene, dhs_id, dhs_start, dhs_end, pval,'NA'))
            else:
                dhs_data.append((gene, dhs_id, dhs_start, dhs_end, pval, cluster_query[3]))
            
    for dhs in dhs_data:
        cluster_id = dhs[5]
        dhsID = dhs[1]
        open_Samples = get_openSamples(cluster_id, dhsDB)
        overlaps = get_overlaps(cluster_id, dhsDB)
        motifs = get_motifs(cluster_id, dhsDB)
        sig_dhsSamples = get_sampleDHS_signal(dhsID, dhsDB)
        output = (dhs[0], dhs[1], dhs[2], dhs[3], dhs[4], sig_dhsSamples[1], \
                      dhs[5], overlaps[1], overlaps[2], overlaps[3], motifs[1],\
                      open_Samples[1], open_Samples[6], open_Samples[7], \
                      open_Samples[8], open_Samples[9], open_Samples[10], \
                      open_Samples[11])
         gene_dhs.append(output)
    return gene_dhs


def get_SNPGeneDHS(snp_eqtls, dhsDB):
    """ Find DHS for SNP-gene pairs. """
    conn = sqlite3.connect(dhsDB)
    conn.text_factory = str
    cur = conn.cursor()
     
    gene_cor = []
    open_samples = []
    dhs_result = []
    snps_done = []
    snps_dhs_data = []
    genes_done = []
    summary_data = []
    detail_Header = False
    placeholder = []
    print '\t Finding DHS info for...'
    for row in xrange(0, len(snp_eqtls)):
        snp = snp_eqtls[row][0]
        snp_chr = snp_eqtls[row][1]
        snp_pos = snp_eqtls[row][2]
        gene = snp_eqtls[row][3]
        gene_chr = snp_eqtls[row][4]
        gene_start = snp_eqtls[row][5]
        gene_end = snp_eqtls[row][6]
        pval = snp_eqtls[row][8]
        qval = snp_eqtls[row][9]
        cis = snp_eqtls[row][12]
        snp_gene_distance = snp_eqtls[row][13]
        genes_data = []
        tester = snp+gene
        if tester not in placeholder:
            print ('\t\t' + snp + ' : ' + gene)
            placeholder.append(tester)
            if snp not in snps_done:
                snps_done.append(snp)
                snp_dhs = get_snpDHS(snp, snp_chr, snp_pos, dhsDB)
                snps_dhs_data.append(snp_dhs)
            if gene not in genes_done:
                gene_dhs = get_geneDHS(gene, gene_start, gene_end, dhsDB)
                genes_data.append(gene_dhs)
        for gene_rows in genes_data:
            this_snp = []
            for done_snp in snps_dhs_data:
                if snp == done_snp[0]:
                    this_snp = done_snp
            for  dhs in gene_rows:
                to_summary = [snp, snp_chr, snp_pos, this_snp[3], this_snp[6], 
                              this_snp[7], this_snp[11], this_snp[18], dhs[0], 
                              gene_chr, gene_start, gene_end, dhs[1], dhs[2], 
                              dhs[3], dhs[4], dhs[5], dhs[6], dhs[10], dhs[16], 
                              dhs[17]]
                summary_data.append(to_summary)     
                with open(filepath + "/dhs_details.txt", 'a') as detail_file:
                    detail = csv.writer(detail_file, delimiter = '\t')
                    if detail_Header == False:
                        detail_file.seek(0)
                        detail_file.truncate() 
                        detail_Head =['SNP', 'SNP_CHR', 'SNP_POS', 'SNP_DHS_ID',\
                                          'SNP_DHS_START', 'SNP_DHS_END', \
                                          'SNP_DHS_SIGNALS', 'SNP_CLUSTER_ID', \
                                          'SNP_CpG%', 'SNP_PROMOTER%', \
                                          'SNP_CONSERVED%', 'SNP_CLUSTER_TFs',\
                                          'SNP_CLUSTER_DHS_COUNT', \
                                          'SNP_OPEN_SAMPLE_COUNT', \
                                          'SNP_OPEN_CELLTYPE_COUNT', \
                                          'SNP_OPEN_TISSUE_COUNT', \
                                          'SNP_MAX_SAMPLE', 'SNP_MAX_CELLTYPE',\
                                          'SNP_MAX_TISSUE', 'GENE', 'GENE_CHR',\
                                          'GENE_START', 'GENE_END', \
                                          'GENE_DHS_ID', 'GENE_DHS_START', \
                                          'GENE_DHS_END', 'DHS_PVAL', \
                                          'GENE_DHS_SIGNALS', 'GENE_CLUSTER_ID',\
                                          'GENE_CpG%', 'GENE_PROMOTER%', \
                                          'GENE_CONSERVED%', 'GENE_CLUSTER_TFs',\
                                          'GENE_CLUSTER_DHS_COUNT',\
                                          'GENE_OPEN_SAMPLE_COUNT', \
                                          'GENE_OPEN_CELLTYPE_COUNT', \
                                          'GENE_OPEN_TISSUE_COUNT', \
                                          'GENE_MAX_SAMPLE', 'GENE_MAX_CELLTYPE',\
                                          'GENE_MAX_TISSUE', 'P_VALUE', 'Q_VALUE',\
                                          'CIS?', 'SNP_GENE_DISTANCE']
                        detail.writerow(detail_Head)
                        detail_Header = True
                    detail.writerow([snp, snp_chr, snp_pos, this_snp[3], \
                                         this_snp[4], this_snp[5], this_snp[6],\
                                         this_snp[7], this_snp[8], this_snp[9],\
                                         this_snp[10], this_snp[11], \
                                         this_snp[12], this_snp[13], \
                                         this_snp[14], this_snp[15], \
                                         this_snp[16], this_snp[17], \
                                         this_snp[18], dhs[0], gene_chr,\
                                         gene_start, gene_end, dhs[1], dhs[2],\
                                         dhs[3], dhs[4], dhs[5], dhs[6],dhs[7],\
                                         dhs[8], dhs[9], dhs[10], dhs[11], \
                                         dhs[12], dhs[13], dhs[14], dhs[15], \
                                         dhs[16], dhs[17], pval, qval, cis, \
                                         snp_gene_distance])


            row = row+1

    with open(filepath + "/dhs_summary.txt", 'wb') as summary_file:
        summary = csv.writer(summary_file, delimiter = '\t')
        summary.writerow(['SNP', 'SNP_CHR', 'SNP_POS', 'SNP_DHS_ID',\
                              'SNP_DHS_SIGNALS', 'SNP_CLUSTER_ID', 'SNP_TFs',\
                              'SNP_MAX_TISSUE', 'GENE','GENE_CHR','GENE_START',\
                              'GENE_END', 'GENE_DHS_ID', 'GENE_DHS_START',\
                              'GENE_DHS_END', 'DHS_PVAL', 'GENE_DHS_SIGNALS',\
                              'GENE_CLUSTER_ID', 'GENE_CLUSTER_TFs', \
                              'MAX_CELLTYPE', 'MAX_TISSUE'])
        summary.writerows(summary_data)
    summary_file.close()

def match_tissues(code3d_data, dhs_filepath, dhsDB):
    placeholder = []
    snp_gene_data = []    # Unique SNP-gene pairs and eQTL tissues
    dhs_raw = []
    dhs_processed = {}
    conn = sqlite3.connect(dhsDB)
    conn.text_factory = str
    cur = conn.cursor()
    for line in code3d_data:
        snp = line[0]
        gene = line[3]
        gtex_tissue = line[7]
        hic_cell_line = line[10]
        gtex_max_tissue = line[15]
        snp_chr = line[1]
        snp_pos = line[2]
        gene_chr = line[4]
        pval = line[8]
        qval = line[9]
        snp_gene_distance = line[13]
        tester = snp + '-' + gene
        if tester not in placeholder:
            snp_gene_data.append([snp, snp_chr, snp_pos, gene, gene_chr, pval, \
                                      qval, snp_gene_distance,hic_cell_line, \
                                      gtex_tissue, gtex_max_tissue])
            placeholder.append(tester)
    for pair in snp_gene_data:
        tissues = []
        for hic_rec in code3d_data:
            if (hic_rec[0] == pair[0] and hic_rec[3] == pair[3]):
                tissues.append(hic_rec[7])
        pair.append(tissues)


    gene_dhs = {}  # genes, their DHS and their open cells        
    gene_celltypes = {}    # genes, celltypes and the number of DHS they have
    snp_dhs = {}   # SNPs and their celltypes if in a DHS
    snp_info = []
    snps_done = []
    print ('\t Comparing open chromatin states of HiC and GTEx tissues... ')
    with open(dhs_filepath, 'r') as dhsfile:
        dhs = csv.reader(dhsfile, delimiter = '\t')
        for code3d in snp_gene_data:
            c_snp = code3d[0]
            c_gene = code3d[3]
            # print ('\t\t' + c_snp  + '-' + c_gene + ' pair')
            dhs_celltypes = {}    # Gene DHS and their open cells
            snp_dhs_celltypes = {} # SNP DHS and their open cells
            for dhs_rec in dhs:
                #dhs_rec = dhs_rec.split('\t')
                snp = dhs_rec[0]
                snp_dhs_id = dhs_rec[3]
                snp_dhs_celltypes = dhs_rec[6]
                snp_max_celltype = dhs_rec[17]
                snp_max_tissue = dhs_rec[18]
                gene = dhs_rec[19]
                gene_dhs_id = dhs_rec[23]
                gene_dhs_celltypes = dhs_rec[27]
                gene_dhs_celltype = dhs_rec[38]
                gene_max_tissue = dhs_rec[39]

                if (snp == c_snp and gene == c_gene):
                    #if c_snp == b:
                    #    print b, 'yes'
                    for cell in gene_dhs_celltypes.split():
                        start = 1
                        celltypes = []
                        for i in xrange(start, len(cell)-1):
                            if cell[i] == ',' and len(cell) -i >=3:
                                cell_dhs_signal = cell[start:i].replace("'", "")
                                cell_dhs_signal = cell_dhs_signal.strip()\
                                    .split(':')
                                celltype = cell_dhs_signal[0] 
                                celltypes.append(celltype)
                                start = i + 1
                            elif len(cell)-i <3:
                                cell_dhs_signal = cell[start:len(cell)-1].\
                                    replace("'", "")
                                cell_dhs_signal = cell_dhs_signal.strip().\
                                    split(':')
                                celltype = cell_dhs_signal[0] 
                                celltypes.append(celltype)
                        dhs_celltypes[gene_dhs_id] = celltypes
        
                        
                snp_celltypes = []
                if snp_dhs_celltypes != 'NA' and snp_dhs_celltypes != "['NA']":
                    start = 1
                    for i in xrange(start, len(snp_dhs_celltypes)-1):
                        if snp_dhs_celltypes[i] == ',' and \
                                len(snp_dhs_celltypes) -i >=3:
                            cell_dhs_signal = snp_dhs_celltypes[start:i].\
                                replace("'", "")
                            cell_dhs_signal = cell_dhs_signal.strip().split(':')
                            celltype = cell_dhs_signal[0] 
                            dhs_signal = float(cell_dhs_signal[1])
                            snp_celltypes.append([celltype, dhs_signal])
                            start = i + 1
                        elif len(snp_dhs_celltypes)-i <3:
                            cell_dhs_signal = snp_dhs_celltypes\
                                [start:len(snp_dhs_celltypes)-1].replace("'", "")
                            cell_dhs_signal = cell_dhs_signal.strip().split(':')
                            celltype = cell_dhs_signal[0]
                            if len(cell_dhs_signal) == 2:
                                dhs_signal = float(cell_dhs_signal[1])
                                snp_celltypes.append([celltype, dhs_signal])
                else:
                    snp_celltypes.append('NA')
                snp_info.append([snp, snp_dhs_id, snp_celltypes])

            gene_dhs[c_gene] = dhs_celltypes
    
        # Get tissues for SNP open celltypes 
    for row in snp_info:
        cells = row[2]
        if cells != ['NA']:
            new_cells = []
            for i in xrange(0, len(cells)):
                cell = cells[i][0]
                cur.execute("SELECT tissue from concordance where canonical \
                                LIKE ? ;", (cell,))
                tissue = cur.fetchone()
                if tissue:
                    new_cells.append([cells[i][0], cells[i][1], tissue[0]])
                else:
                    new_cells.append([cells[i][0], cells[i][1], 'NA'])
            row[2] = new_cells

    
        # get number (and percentage) of a gene's DHSs in cells
    dhs_info = []
    for gene in gene_dhs.keys():
        cell_count = {}
        cell_dhs = {}
        dhs_num = 0
        for dhs in gene_dhs[gene].keys():
            for cell in gene_dhs[gene][dhs]:
                if cell not in cell_count:
                    cell_count[cell] = 1 
                else:
                    cell_count[cell] += 1 
            dhs_num += 1
        for celltype in cell_count:
            cell_dhs_num = cell_count[celltype]
            cell_dhs_percent = float(cell_dhs_num) / dhs_num 
            cell_dhs[celltype] = [dhs_num, cell_dhs_num, \
                                      round(cell_dhs_percent, 2)]
        gene_celltypes[gene]= cell_dhs 

        #   print cell, dhs_percent
        if cell_dhs:
            for hic in snp_gene_data:
                snp = hic[0]
                snp_chr = hic[1]
                snp_pos = hic[2]
                gene_chr = hic[4]
                pval = hic[5]
                qval = hic[6]
                hic_celltype = hic[8]
                snp_gene_distance = hic[7]
                eqtl_tissue = hic[9]
                gtex_max_tissue = hic[10]
                cis = ''
                if not snp_gene_distance == 'NA':
                    if (snp_chr == gene_chr and int(snp_gene_distance) \
                            <= 1000000):
                        cis = 'Cis'
                    else:
                        cis = 'Trans'
                else:
                    cis = 'Trans'
                if gene == hic[3]:
                    hic_tissues = []
                    hic_percent = []
                    for hic_cell in hic_celltype.split(','):
                        
                        cur.execute("SELECT tissue FROM concordance "  
                                       "WHERE canonical LIKE ?;", (hic_cell,))
                        hic_tissue = cur.fetchone()
                        if hic_tissue:
                            hic_tissues.append(hic_tissue[0])
                        else:
                            hic_tissues.append(hic_cell)
                        if hic_cell in cell_dhs:
                            hic_percent.append(cell_dhs[hic_cell][2])
                        else:
                            hic_percent.append('NA')

                    dhs_percent  = cell_dhs[cell][2]
                    max_cell = max(cell_dhs.keys(), \
                                       key = (lambda dhs_percent: \
                                                  cell_dhs[dhs_percent]))
                    max_percent = cell_dhs[max_cell][2]
                        # Correct H1-heSC, LNCaP_andra in concordance table
                    cur.execute("SELECT tissue  FROM concordance WHERE canonical \
                                    LIKE ?;", (max_cell,))
                    max_tissue = cur.fetchone()
                    if max_tissue is None:
                        max_tissue = max_cell
                    else:
                        max_tissue = max_tissue[0]
                    dhs_info.append([snp, snp_chr, gene, gene_chr, cis, \
                                         hic_celltype, hic_percent, hic_tissues,\
                                         max_cell, max_percent, max_tissue, \
                                         eqtl_tissue, gtex_max_tissue])
        else:
            for hic in snp_gene_data:
                snp = hic[0]
                snp_chr = hic[1]
                snp_pos = hic[2]
                gene_chr = hic[4]
                pval = hic[5]
                qval = hic[6]
                hic_celltype = hic[8]
                snp_gene_distance = hic[7]
                eqtl_tissue = hic[9]
                gtex_max_tissue = hic[10]
                cis = ''
                if not snp_gene_distance == 'NA':
                    if (snp_chr == gene_chr and int(snp_gene_distance) \
                            <= 1000000):
                        cis = 'Cis'
                    else:
                        cis = 'Trans'
                else:
                    cis = 'Trans'

                if gene == hic[3]:
                    hic_tissues = []
                    for hic_cell in hic_celltype.split(','):
                        cur.execute("SELECT tissue FROM concordance \
                                      WHERE canonical LIKE ?;", (hic_celltype,))
                        hic_tissue = cur.fetchone()
                        if hic_tissue:
                            hic_tissues.append(hic_tissue[0])
                        else:
                            hic_tissues.append(hic_cell)
                    dhs_info.append([snp, snp_chr, gene, gene_chr, cis, \
                                         hic_celltype, 'NA', hic_tissues, \
                                         'NA', 'NA', 'NA', eqtl_tissue, \
                                         gtex_max_tissue])

    with open(filepath + '/match.txt', 'wb') as match:
        w = csv.writer(match, delimiter = '\t')
        w.writerow(['SNP', 'SNP_CHR', 'SNP_DHS_ID', 'SNP_OPEN_CELLTYPES', \
                        'GENE', 'GENE_CHR', 'CIS','HiC_CELLTYPE', 'HiC_DHS%',\
                        'HiC_tissue', 'MAX_DHS_CELL', 'MAX_DHS_CELL%', \
                        'MAX_DHS_TISSUE', 'eQTL_TISSUE', 'GTeX_MAX_TISSUE'])
        written = []
        for gene in dhs_info:
            for snp in snp_info:
                test = snp+gene
                if snp[0] == gene[0] and test not in written:
                    written.append(test)
                    p  = snp[0], gene[1], snp[1], snp[2], gene[2], gene[3], \
                        gene[4], gene[5], gene[6], gene[7], gene[8], gene[9], \
                        gene[10], gene[11], gene[12]
                    w.writerow(p)
    match.close()
                


if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, \
                            help = "The \'summary.txt\' output from the hiC \
                            query developed by Cam.")
    parser.add_argument("-d", "--dhsDB", default = "dhs.db", \
                            help = "The DNA Regulatory Elements database.")
    args = parser.parse_args()
    filepath = set_filepath(args.input)
    snp_genes = sig_SNP_genes(args.input)
    
    get_SNPGeneDHS(snp_genes, args.dhsDB)
    dhs_filepath = filepath + '/dhs_details.txt'
    match_tissues(snp_genes, dhs_filepath, args.dhsDB)

#!usr/bin/env python
import sys, sqlite3, csv, os, argparse

    

def prep_input(input_fp):
    snp_eqtls = []
    if os.path.isfile(input_fp):

        with open(input_fp, 'rb') as infile:
            for line in infile:
                line = line.strip().split()
                row = line[0:7]
                if row not in snp_eqtls:
                    snp_eqtls.append(row)
    return snp_eqtls





def get_openCellTypes(cluster_id, dhsDB):
    conn = sqlite3.connect(dhsDB)
    conn.text_factory =str
    cur = conn.cursor()
    cellTypes = []
    
    output = [cluster_id, 'NA', 'NA', 'NA']
    if isinstance(cluster_id, int):
        cur.execute("SELECT * FROM openCellTypes WHERE cluster_id = ?;", (cluster_id,))
        data = cur.fetchone()
        if isinstance(data, tuple):
            cellTypes =data
        else:
            cellTypes =output
    else:
        cellTypes =output
 
    return cellTypes
    cur.close()
    conn.close()
 

def get_openSamples(cluster_id, dhsDB):
    conn = sqlite3.connect(dhsDB)
    conn.text_factory =str
    cur = conn.cursor()
    openSamples = []
    #print "\t\t Gathering info on the open cell types and tissues in which the SNPs are found..."

    output = [cluster_id, 'NA', 'NA', 'NA','NA', 'NA', 'NA', 'NA','NA', 'NA', 'NA', 'NA']
    if isinstance(cluster_id, int):
        cur.execute("SELECT * FROM openSamples WHERE cluster_id = ?;", (cluster_id,))
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
    #print "\t\t Verifying if SNP lies in a CpG or promoter region... "

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
    conn.text_factory =str
    cur = conn.cursor()
    motifs = []
    #print "\t\t Verifying if SNP lies in motif..."

    # TODO: Modify to retrieve all significant (e < 10-6) TFs?

    if isinstance(cluster_id, int):
        cur.execute("SELECT * FROM motifJaspar WHERE cluster_id = ?;", (cluster_id,))
        data = cur.fetchall()
        tf_list = []
        if len(data) is 0:
            tf_list = 'NA'
        else:
            for row in data:
                tf_list.append((row[4], row[5]))
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

    #print ("\t\t Scanning for cell types with significant DHS signals...")

    if not isinstance(dhs_id, int):
        sig_Samples = [dhs_id, ['NA']]
            
    else:
        cur.execute("SELECT * FROM dhs112 WHERE rowid =?;", (dhs_id+1,))
        dhs_signal = cur.fetchone()
        head = [description[0] for description in cur.description]
        sig_dhs =[]
           
        for i in xrange(3, len(dhs_signal)):
            if (dhs_signal[i] >= 0.1):
                sig_dhs.append((head[i],dhs_signal[i]))
            i = i+1
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

    
    snp_chr = "chr"+snp_chr

    cur.execute("SELECT rowid, chr, start, end FROM dhs112 WHERE chr = ? AND start <=?  AND end >=?;",(snp_chr,snp_pos, snp_pos))
    data = cur.fetchone()
    if (data == None):
        dhs_data = [snp, snp_chr, snp_pos, 'NA', 'NA', 'NA']
        #dhs_data.append(to_data)
                
    else:
        dhs_data = [snp, snp_chr, snp_pos, data[0]-1, data[2], data[3]]
        #dhs_data.append(to_data)

    
    dhs_id = dhs_data[3]
    if isinstance(dhs_id, int):
        cur.execute("SELECT chr, start, end, refined_cluster FROM dhsCluster LIMIT 1 OFFSET "+ str(dhs_id-1)+";")
        data = cur.fetchone()
        dhs_data.append(data[3])

        #print dhs_data[i]
    else:
        dhs_data.append('NA')


    cluster_id = dhs_data[6]
    
    open_Samples = get_openSamples(cluster_id, dhsDB)
    overlaps = get_overlaps(cluster_id, dhsDB)
    motifs = get_motifs(cluster_id, dhsDB)
    sig_dhsSamples = get_sampleDHS_signal(dhs_id, dhsDB)
    
    snp_dhs = (snp, snp_chr, snp_pos, dhs_id, dhs_data[4], dhs_data[5], sig_dhsSamples[1], cluster_id, overlaps[1], overlaps[2], overlaps[3], motifs[1], open_Samples[1], open_Samples[6], open_Samples[7], open_Samples[8], open_Samples[9], open_Samples[10], open_Samples[11])

    return snp_dhs
    

    cur.close()
    conn.close()

def get_geneDHS(gene, dhsDB):
    conn = sqlite3.connect(dhsDB)
    conn.text_factory = str
    cur = conn.cursor()
    gene_dhs = []
    dhs_data = []
    open_Samples = []
    overlaps = []
    motifs = []
    sig_dhsSamples = []
    
    print("\t\t Finding DHSs in "+gene)
    cur.execute("SELECT gene_name, dhs_id, dhs_start, dhs_end, pval from gene_correlations_p05 WHERE pval >= 0.95 AND gene_name = ?;", (gene,))
    
    cor_query = cur.fetchall()
    if cor_query is None:
        dhs_data.append((gene, 'NA', 'NA', 'NA', 'NA', 'NA'))
    else:
        for item in cor_query:
            dhs_id = item[1]
            dhs_start = int(item[2])
            dhs_end = int(item[3])
            pval = 1-item[4]
            cluster_id = 0

            cur.execute("SELECT chr, start, end, refined_cluster FROM dhsCluster LIMIT 1 OFFSET "+ str(dhs_id-1)+";")
            cluster_query = cur.fetchone()
            if cluster_query is None:
                #cluster_id = 'NA'
                dhs_data.append((gene, dhs_id, dhs_start, dhs_end, pval,'NA'))
            else:
                #cluster_id = cluster_query[3]
                dhs_data.append((gene, dhs_id, dhs_start, dhs_end, pval, cluster_query[3]))
            
    for dhs in dhs_data:
        cluster_id = dhs[5]
        dhsID = dhs[1]
        open_Samples = get_openSamples(cluster_id, dhsDB)
        overlaps = get_overlaps(cluster_id, dhsDB)
        motifs = get_motifs(cluster_id, dhsDB)
        sig_dhsSamples = get_sampleDHS_signal(dhsID, dhsDB)


        output = (dhs[0], dhs[1], dhs[2], dhs[3], dhs[4], sig_dhsSamples[1], dhs[5], overlaps[1], overlaps[2], overlaps[3], motifs[1], open_Samples[1], open_Samples[6], open_Samples[7], open_Samples[8], open_Samples[9], open_Samples[10],open_Samples[11]) 
        gene_dhs.append(output)
    return gene_dhs


def get_SNPGeneDHS(snp_eqtls, dhsDB):
    conn = sqlite3.connect(dhsDB)
    conn.text_factory = str
    cur = conn.cursor()
     
    gene_cor = []
    open_samples = []
    dhs_result = []
    snps_done =[]
    snps_dhs_data = []
    genes_done = []
    detail_Header = False

    with open(output_dir+"dhs_summary.txt", 'wb') as summary_file:
        summary = csv.writer(summary_file, delimiter = '\t')
        summary.writerow(['SNP', 'SNP_CHR', 'SNP_POS', 'SNP_DHS_ID', 'SNP_DHS_SIGNALS', 'SNP_CLUSTER_ID', 'SNP_TFs', 'SNP_MAX_TISSUE', 'GENE', 'GENE_CHR', 'GENE_START', 'GENE_END', 'GENE_DHS_ID', 'GENE_DHS_START', 'GENE_DHS_END', 'DHS_PVAL', 'GENE_DHS_SIGNALS', 'GENE_CLUSTER_ID', 'GENE_CLUSTER_TFs', 'MAX_CELLTYPE', 'MAX_TISSUE' ] )
        for row in xrange(1,len(snp_eqtls)):
            snp = snp_eqtls[row][0]
            snp_chr = snp_eqtls[row][1]
            snp_pos = snp_eqtls[row][2]
            gene = snp_eqtls[row][3]
            gene_chr = snp_eqtls[row][4]
            gene_start = snp_eqtls[row][5]
            gene_end = snp_eqtls[row][6]
            genes_data = []
        
            if snp not in snps_done:
                print ("\t Checking if "+snp+" lies in a DHS")
                snp_dhs = get_snpDHS(snp, snp_chr, snp_pos, dhsDB)
                snps_dhs_data.append(snp_dhs)
                snps_done.append(snp)
                if gene not in genes_done:
                    gene_dhs = get_geneDHS(gene, dhsDB)
                    genes_data.append(gene_dhs)
                    genes_done.append(gene)

            else:
                if gene not in genes_done:
                    gene_dhs = get_geneDHS(gene, dhsDB)
                    genes_data.append(gene_dhs)
                    genes_done.append(gene)
            for gene_rows in genes_data:
                for  dhs in gene_rows:
                    summary.writerow([snp, snp_chr, snp_pos, snp_dhs[3], snp_dhs[6], snp_dhs[7], snp_dhs[11], snp_dhs[18], dhs[0], gene_chr, gene_start, gene_end, dhs[1], dhs[2], dhs[3], dhs[4], dhs[5], dhs[6], dhs[10], dhs[16], dhs[17]])
                    
                    
                    with open(output_dir+"/dhs_details.txt", 'a') as detail_file:
                        detail = csv.writer(detail_file, delimiter = '\t')
                        if detail_Header == False:
                            detail_file.seek(0)
                            detail_file.truncate() 
                            detail_Head =['SNP', 'SNP_CHR', 'SNP_POS', 'SNP_DHS_ID', 'SNP_DHS_START', 'SNP_DHS_END', 'SNP_DHS_SIGNALS', 'SNP_CLUSTER_ID', 'SNP_CpG%', 'SNP_PROMOTER%', 'SNP_CONSERVED%', 'SNP_CLUSTER_TFs', 'SNP_CLUSTER_DHS_COUNT', 'SNP_OPEN_SAMPLE_COUNT', 'SNP_OPEN_CELLTYPE_COUNT', 'SNP_OPEN_TISSUE_COUNT', 'SNP_MAX_SAMPLE', 'SNP_MAX_CELLTYPE',  'SNP_MAX_TISSUE', 'GENE', 'GENE_CHR', 'GENE_START', 'GENE_END', 'GENE_DHS_ID', 'GENE_DHS_START', 'GENE_DHS_END', 'DHS_PVAL', 'GENE_DHS_SIGNALS', 'GENE_CLUSTER_ID','GENE_CpG%', 'GENE_PROMOTER%', 'GENE_CONSERVED%', 'GENE_CLUSTER_TFs', 'GENE_CLUSTER_DHS_COUNT', 'GENE_OPEN_SAMPLE_COUNT', 'GENE_OPEN_CELLTYPE_COUNT', 'GENE_OPEN_TISSUE_COUNT', 'GENE_MAX_SAMPLE', 'GENE_MAX_CELLTYPE',  'GENE_MAX_TISSUE' ] 
                            detail.writerow(detail_Head)
                            detail_Header = True
                        detail.writerow([snp, snp_chr, snp_pos, snp_dhs[3], snp_dhs[4], snp_dhs[5], snp_dhs[6], snp_dhs[7], snp_dhs[8], snp_dhs[9], snp_dhs[10], snp_dhs[11], snp_dhs[12], snp_dhs[13], snp_dhs[14], snp_dhs[15], snp_dhs[16], snp_dhs[17], snp_dhs[18], dhs[0], gene_chr, gene_start, gene_end, dhs[1], dhs[2], dhs[3], dhs[4], dhs[5], dhs[6],dhs[7], dhs[8], dhs[9], dhs[10], dhs[11], dhs[12], dhs[13], dhs[14], dhs[15], dhs[16], dhs[17]])
                        
                    detail_file.close()
    
        row = row+1
    summary_file.close()


def resolve_output_fp(input_fp):
    output_fp = input_fp
    if ('/' in output_fp):
        while not (output_fp.endswith('/')):
            output_fp = output_fp[:len(output_fp)-1]
    else:
        output_fp = ''
    print output_fp+"dhs_summary.txt"
    return output_fp


 
 

    

if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required = True, help = "The \'summary.txt\' output from the hiC query developed by Cam.")
    parser.add_argument("-d", "--dhsDB", default = "/mnt/3dgenome/projects/tfad334/dhs/dhs-mining/dhs.db", help = "The DNA Regulatory Elements database.")
    args = parser.parse_args()

    output_dir = resolve_output_fp(args.input)

    snp_genes = prep_input(args.input)
    
    get_SNPGeneDHS(snp_genes, args.dhsDB)


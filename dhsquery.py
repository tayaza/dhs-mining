#usr/bin/env python
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
    



def create_corrTB(corr_file):
    corr_conn = sqlite3.connect('dhs.db')
    corr_conn.text_factory = str
    corr_cur = corr_conn.cursor()

    corr_cur.execute("DROP TABLE IF EXISTS gene_correlations_p05")
    corr_cur.execute("CREATE TABLE gene_correlations_p05 (dhs_chr TEXT, dhs_start INTEGER, dhs_stop INTEGER,dhs_id INTEGER, gene_chr TEXT, gene_start INTEGER, gene_stop INTEGER, gene_name TEXT, metaprobeset_id REAL, ensemblID TEXT, cor REAL, pval REAL )" )

    with open(corr_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        corr_cur.executemany("INSERT INTO gene_correlations_p05 VALUES (?,?,?,?,?,?,?,?,?,?,?,?);", (row,))
    corr_conn.commit()
    corr_conn.close()


def create_overlapTB(overlap_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    cur.execute("DROP TABLE IF EXISTS overlap")
    cur.execute("CREATE TABLE overlap (clusterID integer,cpg INTEGER,promoter INTEGER,phastVert)" )

    with open(overlap_file, 'rb') as infile:
        data = [row for row in csv.reader(infile.read().splitlines())]
    for row in data:
        row = row[0].split()
        cur.executemany("INSERT INTO overlap VALUES (?,?,?,?);", (row,))
    conn.commit()
    conn.close()


def create_concordanceTB(concordance_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    cur.execute("DROP TABLE IF EXISTS concordance")
    cur.execute("CREATE TABLE concordance (canonical TEXT,celltype TEXT,tissue TEXT,malignant TEXT, sex TEXT, description TEXT)" )

    with open(concordance_file, 'rb') as infile:
        data = [row for row in csv.reader(infile.read().splitlines())]

    for row in data:
        print i, row[5]
        cur.executemany("INSERT INTO concordance VALUES (?,?,?,?,?,?);", (row,))
    conn.commit()
    conn.close()


def create_openCellTypesTB(openCellTypes_file):
    corr_conn = sqlite3.connect('dhs.db')
    corr_conn.text_factory = str
    corr_cur = corr_conn.cursor()

    corr_cur.execute("DROP TABLE IF EXISTS openCellTypes")
    corr_cur.execute("CREATE TABLE openCellTypes (clusterID INTEGER,openSample TEXT, opoenCellType TEXT, openTissue TEXT)" )

    with open(openCellTypes_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        corr_cur.executemany("INSERT INTO openCellTypes VALUES (?,?,?,?);", (row,))
    corr_conn.commit()
    corr_conn.close()


#TODO Learn how to work with BED files
def create_dhsPredictorsTB(dhsPredictors_file):
    corr_conn = sqlite3.connect('dhs.db')
    corr_conn.text_factory = str
    corr_cur = corr_conn.cursor()

    corr_cur.execute("DROP TABLE IF EXISTS dhsPredictors")
    corr_cur.execute("CREATE TABLE dhsPredictors (chr TEXT, start INTEGER, stop INTEGER, dhs_id INTEGER, cluster_id INTEGER, tissue TEXT)" )

    with open(dhsPredictors_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        corr_cur.executemany("INSERT INTO dhsPredictors VALUES (?,?,?,?,?,?);", (row,))
    corr_conn.commit()
    corr_conn.close()


def create_openSamplesTB(openSamples_file):
    corr_conn = sqlite3.connect('dhs.db')
    corr_conn.text_factory = str
    corr_cur = corr_conn.cursor()

    corr_cur.execute("DROP TABLE IF EXISTS openSamples")
    corr_cur.execute("CREATE TABLE openSamples (clusterID INTEGER, dhs_count INTEGER, chr TEXT,  start INTEGER, stop INTEGER, dhs_id INTEGER, openSampleCount INTEGER, openCelltypeCount INTEGER, openTissueCount INTEGER, maxSample TEXT, maxCellType TEXT, maxTissue TEXT)" )

    with open(openSamples_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        corr_cur.executemany("INSERT INTO openSamples VALUES (?,?,?,?,?,?,?,?,?,?,?,?);", (row,))
    corr_conn.commit()
    corr_conn.close()

def create_motifJasparTB(motifJaspar_file):
    corr_conn = sqlite3.connect('dhs.db')
    corr_conn.text_factory = str
    corr_cur = corr_conn.cursor()

    corr_cur.execute("DROP TABLE IF EXISTS motifJaspar")
    corr_cur.execute("CREATE TABLE motifJaspar (clusterID INTEGER, iMotif INTEGER, siteCount INTEGER,  eDiscovery REAL, TF TEXT, eMatch REAL,database TEXT)" )

    with open(motifJaspar_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        corr_cur.executemany("INSERT INTO motifJaspar VALUES (?,?,?,?,?,?,?);", (row,))
    corr_conn.commit()
    corr_conn.close()



def create_malignantTissuesTB(malignantTissues_file):
    corr_conn = sqlite3.connect('dhs.db')
    corr_conn.text_factory = str
    corr_cur = corr_conn.cursor()

    corr_cur.execute("DROP TABLE IF EXISTS malignantTissues")
    corr_cur.execute("CREATE TABLE malignantTissues (canonical TEXT,isMalig TEXT, trueTissue TEXT, myGuess TEXT,isCorrect TEXT,brain REAL, endothelial REAL, epithelial REAL, fibroblast REAL, hematopoietic REAL, muscle REAL, stem REAL)" )

    with open(malignantTissues_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    i=0
    for row in data:
        row = row[0].split()
        #i = i+1
        #print i, row[11]
        corr_cur.executemany("INSERT INTO malignantTissues VALUES (?,?,?,?,?,?,?,?,?,?,?,?);", (row,))
    corr_conn.commit()
    corr_conn.close()



def create_sexTB(sex_file):
    corr_conn = sqlite3.connect('dhs.db')
    corr_conn.text_factory = str
    corr_cur = corr_conn.cursor()

    corr_cur.execute("DROP TABLE IF EXISTS sex")
    corr_cur.execute("CREATE TABLE sex (canonical TEXT,tissue TEXT, trueSex TEXT, myGuess TEXT,isCorrect TEXT,F REAL, M REAL)" )

    with open(sex_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        corr_cur.executemany("INSERT INTO sex VALUES (?,?,?,?,?,?,?);", (row,))
    corr_conn.commit()
    corr_conn.close()

### END OF DATABASE PROCESSING

###
dhsPredictors = 'dhs_data/TableS07-dhsPredictors_v2.bed' 

def create_dhsPredTB(dhsPredictors):
    content = []
    with open(dhsPredictors, 'rb') as f:
        for line in f:
            content.append(line.strip().split())
        return content

dhsPred = create_dhsPredTB(dhsPredictors)
#for row in dhsPred:
#    print row
####

def get_dhs(genes):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()
    
    gene_cor = []
    open_samples = []
    dhs_result = []
    
   

    for gene in genes:
        cur.execute("SELECT DISTINCT gene_name, gene_chr, gene_start, gene_stop, dhs_id, dhs_start, dhs_stop, pval from gene_correlations_p05  WHERE gene_name = ? AND pval >= 0.95;", (gene,))
        data = cur.fetchall()
        for item in data:
            rec = [item[0], item[1], item[2], item[3], item[4], item[5], item[6], 1-item[7]]
            gene_cor.append(rec)
   

   
    cur.execute("SELECT * FROM openSamples;")
    os_result = cur.fetchall()
    print "\t\t Fetching open samples for DHS."
    for row in os_result:
        open_samples.append(row)
   
    for gene in gene_cor:
        if ((gene[2] == 'NA') or (gene[3] == 'NA')):
            data=[gene[0], gene[4], gene[7], gene[1], gene[2], gene[3], 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']
            dhs_result.append(data)
        else:
            gStart = int(gene[5])
            gStop = int(gene[6])
            for os_row in open_samples:
                cStart = int(os_row[3])
                cStop = int(os_row[4])
                if ((gene[1] == os_row[2]) and ((abs(gStart - cStart) <= 100000) or (abs(gStart - cStop) <= 100000) or (abs(gStop - cStart) <= 100000) or (abs(gStop - cStop) <= 100000))):
                    data = [gene[0], gene[4], gene[7], gene[1], gene[2], gene[3], os_row[0], os_row[1], os_row[6],os_row[7],os_row[8], os_row[9], os_row[10], os_row[11]]
                    dhs_result.append(data)
    return dhs_result

def output_dhs(dhs_result):
    with open('dhs_outfile.txt2', 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        header =['gene_name', 'dhs_id', 'pval', 'chr', 'dhs_start', 'dhs_end', 'cluster_id', 'dhs_count',  'openSampleCount', 'openCellTypeCount', 'openTissueCount', 'maxSample', 'maxCellType', 'maxTissue']
        writer.writerow(header)
        writer.writerows(dhs_result)
 
 

    

    
    
input_file = 'MWC/obesity/obesity_genes.txt'

#correlation_input = 'dhs_data/allGeneCorrelations100000.p05_v3.txt'
#overlap_input = 'dhs_data/TableS05-overlapSummary.txt'
#concordance_input = 'dhs_data/TableS01-concordanceSummary.csv'
#openCellTypes_input = 'dhs_data/TableS04-cluster-to-openCellTypes_v3.txt'
#openSamples_input = 'dhs_data/TableS02-openSamples_v3.txt'
#motifJaspar_input = 'dhs_data/TableS09-motifJaspar.txt'
#malignantTissues_input = 'dhs_data/TableS06-results_tissueMalig_v2.txt'
#sex_input = 'dhs_data/TableS08-results_sex_v3.txt'

#dhsPredictors = 'dhs_data/TableS07-dhsPredictors_v2.bed' #TODO learn how to handle BED files
dhs112_file = 'dhs_data/dhs112_v3.bed'

genes = prep_gene_list(input_file)

#create_corrTB(correlation_input)
#create_overlapTB(overlap_input)
#create_concordanceTB(concordance_input)
#create_openCellTypesTB(openCellTypes_input)
#create_dhsPredictorsTB(dhsPredictors_input)
#create_openSamplesTB(openSamples_input)
#create_motifJasparTB(motifJaspar_input)
#create_malignantTissuesTB(malignantTissues_input)
#create_sexTB(sex_input)

#dhsPred = create_dhsPredTB(dhsPredictors)
query = get_dhs(genes)
output_dhs(query)

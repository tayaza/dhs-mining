#usr/bin/env python
import sys, sqlite3, csv



def create_corrTB(corr_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()
    
    print("\t\t Creating Gene Correlations table...")
    cur.execute("DROP TABLE IF EXISTS gene_correlations_p05")
    cur.execute("CREATE TABLE gene_correlations_p05 (dhs_chr TEXT, dhs_start INTEGER, dhs_end INTEGER,dhs_id INTEGER, gene_chr TEXT, gene_start INTEGER, gene_end INTEGER, gene_name TEXT, metaprobeset_id REAL, ensemblID TEXT, cor REAL, pval REAL )" )

    with open(corr_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        cur.executemany("INSERT INTO gene_correlations_p05 VALUES (?,?,?,?,?,?,?,?,?,?,?,?);", (row,))
    conn.commit()
    conn.close()


def create_overlapTB(overlap_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    print("\t\t Creating Overlap table...")

    cur.execute("DROP TABLE IF EXISTS overlap")
    cur.execute("CREATE TABLE overlap (cluster_id integer,cpg INTEGER,promoter INTEGER,phastVert)" )

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

    print("\t\t Creating Concordance table...")

    cur.execute("DROP TABLE IF EXISTS concordance")
    cur.execute("CREATE TABLE concordance (canonical TEXT,celltype TEXT,tissue TEXT,malignant TEXT, sex TEXT, description TEXT)" )

    with open(concordance_file, 'rb') as infile:
        data = [row for row in csv.reader(infile.read().splitlines())]

    for row in data:
        cur.executemany("INSERT INTO concordance VALUES (?,?,?,?,?,?);", (row,))
    
    conn.commit()
    conn.close()


def create_openCellTypesTB(openCellTypes_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    cur.execute("DROP TABLE IF EXISTS openCellTypes")
    cur.execute("CREATE TABLE openCellTypes (cluster_id INTEGER,openSample TEXT, opoenCellType TEXT, openTissue TEXT)" )

    with open(openCellTypes_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    
    for row in data:
        row = row[0].split()
        cur.executemany("INSERT INTO openCellTypes VALUES (?,?,?,?);", (row,))
    
    conn.commit()
    conn.close()



def create_dhsPredictorsTB(dhsPredictors_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    print("\t\t Creating Gene DHS Predictors table...")

    cur.execute("DROP TABLE IF EXISTS dhsPredictors")
    cur.execute("CREATE TABLE dhsPredictors (chr TEXT, start INTEGER, end INTEGER, dhs_id INTEGER, cluster_id INTEGER, tissue TEXT)" )

    with open(dhsPredictors_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        cur.executemany("INSERT INTO dhsPredictors VALUES (?,?,?,?,?,?);", (row,))
    conn.commit()
    conn.close()


def create_openSamplesTB(openSamples_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    print("\t\t Creating Gene Open Samples table...")

    cur.execute("DROP TABLE IF EXISTS openSamples")
    cur.execute("CREATE TABLE openSamples (cluster_id INTEGER, dhs_count INTEGER, chr TEXT,  start INTEGER, end INTEGER, dhs_id INTEGER, openSampleCount INTEGER, openCelltypeCount INTEGER, openTissueCount INTEGER, maxSample TEXT, maxCellType TEXT, maxTissue TEXT)" )

    with open(openSamples_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        cur.executemany("INSERT INTO openSamples VALUES (?,?,?,?,?,?,?,?,?,?,?,?);", (row,))
    conn.commit()
    conn.close()

def create_motifJasparTB(motifJaspar_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    print("\t\t Creating Motif table...")

    cur.execute("DROP TABLE IF EXISTS motifJaspar")
    cur.execute("CREATE TABLE motifJaspar (cluster_id INTEGER, iMotif INTEGER, siteCount INTEGER,  eDiscovery REAL, TF TEXT, eMatch REAL,database TEXT)" )

    with open(motifJaspar_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        cur.executemany("INSERT INTO motifJaspar VALUES (?,?,?,?,?,?,?);", (row,))
    conn.commit()
    conn.close()



def create_malignantTissuesTB(malignantTissues_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    print("\t\t Creating Malignant Tissues table...")

    cur.execute("DROP TABLE IF EXISTS malignantTissues")
    cur.execute("CREATE TABLE malignantTissues (canonical TEXT,isMalig TEXT, trueTissue TEXT, myGuess TEXT,isCorrect TEXT,brain REAL, endothelial REAL, epithelial REAL, fibroblast REAL, hematopoietic REAL, muscle REAL, stem REAL)" )

    with open(malignantTissues_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    i=0
    for row in data:
        row = row[0].split()
        #i = i+1
        #print i, row[11]
        cur.executemany("INSERT INTO malignantTissues VALUES (?,?,?,?,?,?,?,?,?,?,?,?);", (row,))
    conn.commit()
    conn.close()



def create_sexTB(sex_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    print("\t\t Creating Sex table...")

    cur.execute("DROP TABLE IF EXISTS sex")
    cur.execute("CREATE TABLE sex (canonical TEXT,tissue TEXT, trueSex TEXT, myGuess TEXT,isCorrect TEXT,F REAL, M REAL)" )

    with open(sex_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        cur.executemany("INSERT INTO sex VALUES (?,?,?,?,?,?,?);", (row,))
    conn.commit()
    conn.close()


def create_dhs112TB(dhs_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    print("\t\t Creating DHS 112 table...")

    cur.execute("DROP TABLE IF EXISTS dhs112")
    cur.execute("CREATE TABLE dhs112 ('chr' TEXT, 'start' INTEGER, 'end' INTEGER, 'A549' REAL, 'AG04449' REAL, 'AG04450' REAL, 'AG09309' REAL, 'AG09319' REAL, 'AG10803' REAL, 'AoAF' REAL, 'AoSMC_SFM' REAL, 'BE2_C' REAL, 'BJ' REAL, 'CD14' REAL, 'Chorion' REAL, 'CLL' REAL, 'CMK' REAL, 'Colo829' REAL, 'Fibrobl' REAL, 'FibroP' REAL, 'Gliobla' REAL, 'GM06990' REAL, 'GM12864' REAL, 'GM12865' REAL, 'GM12878' REAL, 'GM12891' REAL, 'GM12892' REAL, 'GM18507' REAL, 'GM19238' REAL, 'GM19239' REAL, 'GM19240' REAL, 'H1-hE5C' REAL, 'H7-hE5C' REAL, 'H9-hE5C' REAL, 'HA-c' REAL, 'HA-sp' REAL, 'HAEpiC' REAL, 'HAh' REAL, 'HBMEC' REAL, 'HCF' REAL, 'HCFaa' REAL, 'HCM' REAL, 'HConF' REAL, 'HCPEpiC' REAL, 'HCT-116' REAL, 'HEEpiC' REAL, 'HeLa-S3' REAL, 'HeLa-S3_IFNA' REAL, 'Hepatocytes' REAL, 'HepG2' REAL, 'HFF' REAL, 'HFF_Myc' REAL, 'HGF' REAL, 'HIPEpiC' REAL, 'HL-60' REAL, 'HMEC' REAL, 'HMF' REAL, 'HMVEC-dBl-Ad' REAL, 'HMVEC-dBl-Neo'  REAL, 'HMVEC-dLy-Ad' REAL, 'HMVEC-dLy-Neo' REAL, 'HMVEC-dNeo' REAL, 'HMVEC-LBl' REAL, 'HMVEC-Lly' REAL, 'HMVECdAd' REAL, 'HNPCEpiC' REAL, 'HPAEC' REAL, 'HPAF' REAL, 'HPDE6-E6E7' REAL, 'HPdLF' REAL, 'HPF' REAL, 'HRCE' REAL, 'HRE' REAL, 'HRGEC' REAL, 'HRPEpiC' REAL, 'HSMM' REAL, 'HSMMtube' REAL, 'Htr8' REAL, 'Huh-7' REAL, 'Huh-75' REAL, 'HUVEC' REAL, 'HVMF' REAL, 'iPS' REAL, 'Jurkat' REAL, 'K562' REAL, 'LNCaP' REAL, 'LNCap_andro' REAL, 'MCF-7' REAL, 'MCF-7_hyp_lac' REAL, 'Medullo' REAL, 'Melano' REAL, 'Myometr' REAL, 'NB4' REAL, 'NHA' REAL, 'NHDF-Ad' REAL, 'NHDF-neo' REAL, 'NHEK' REAL, 'NHLF' REAL, 'Ntera2' REAL, 'Osteobl' REAL, 'PA-TU-8988T' REAL, 'PANC-1' REAL, 'PrEC' REAL, 'ProgFib' REAL, 'RPTEC' REAL, 'SAEC' REAL, 'SK-N-SH_RA' REAL, 'SKMC' REAL, 'SKNMC' REAL, 'Stellate' REAL, 'Th1' REAL, 'Th2' REAL, 'Urothelia_UT189' REAL, 'WI-38' REAL, 'WI-38-TAM' REAL)" )

    with open(dhs_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        
        #NOTE: dhs_id is (row_id-1)
        cur.executemany("INSERT INTO dhs112 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);", (row,))
    conn.commit()
    conn.close()


def create_exp112TB(exp_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    print("\t\t Creating Expression table...")

    cur.execute("DROP TABLE IF EXISTS exp112")
    cur.execute("CREATE TABLE exp112 ('chr' TEXT, 'start' INTEGER, 'end' INTEGER, 'geneSymbol' TEXT, 'probeset_id' INTEGER, 'A549' REAL, 'AG04449' REAL, 'AG04450' REAL, 'AG09309' REAL, 'AG09319' REAL, 'AG10803' REAL, 'AoAF' REAL, 'AoSMC_SFM' REAL, 'BE2_C' REAL, 'BJ' REAL, 'CD14' REAL, 'Chorion' REAL, 'CLL' REAL, 'CMK' REAL, 'Colo829' REAL, 'Fibrobl' REAL, 'FibroP' REAL, 'Gliobla' REAL, 'GM06990' REAL, 'GM12864' REAL, 'GM12865' REAL, 'GM12878' REAL, 'GM12891' REAL, 'GM12892' REAL, 'GM18507' REAL, 'GM19238' REAL, 'GM19239' REAL, 'GM19240' REAL, 'H1-hE5C' REAL, 'H7-hE5C' REAL, 'H9-hE5C' REAL, 'HA-c' REAL, 'HA-sp' REAL, 'HAEpiC' REAL, 'HAh' REAL, 'HBMEC' REAL, 'HCF' REAL, 'HCFaa' REAL, 'HCM' REAL, 'HConF' REAL, 'HCPEpiC' REAL, 'HCT-116' REAL, 'HEEpiC' REAL, 'HeLa-S3' REAL, 'HeLa-S3_IFNA' REAL, 'Hepatocytes' REAL, 'HepG2' REAL, 'HFF' REAL, 'HFF_Myc' REAL, 'HGF' REAL, 'HIPEpiC' REAL, 'HL-60' REAL, 'HMEC' REAL, 'HMF' REAL, 'HMVEC-dBl-Ad' REAL, 'HMVEC-dBl-Neo'  REAL, 'HMVEC-dLy-Ad' REAL, 'HMVEC-dLy-Neo' REAL, 'HMVEC-dNeo' REAL, 'HMVEC-LBl' REAL, 'HMVEC-Lly' REAL, 'HMVECdAd' REAL, 'HNPCEpiC' REAL, 'HPAEC' REAL, 'HPAF' REAL, 'HPDE6-E6E7' REAL, 'HPdLF' REAL, 'HPF' REAL, 'HRCE' REAL, 'HRE' REAL, 'HRGEC' REAL, 'HRPEpiC' REAL, 'HSMM' REAL, 'HSMMtube' REAL, 'Htr8' REAL, 'Huh-7' REAL, 'Huh-75' REAL, 'HUVEC' REAL, 'HVMF' REAL, 'iPS' REAL, 'Jurkat' REAL, 'K562' REAL, 'LNCaP' REAL, 'LNCap_andro' REAL, 'MCF-7' REAL, 'MCF-7_hyp_lac' REAL, 'Medullo' REAL, 'Melano' REAL, 'Myometr' REAL, 'NB4' REAL, 'NHA' REAL, 'NHDF-Ad' REAL, 'NHDF-neo' REAL, 'NHEK' REAL, 'NHLF' REAL, 'Ntera2' REAL, 'Osteobl' REAL, 'PA-TU-8988T' REAL, 'PANC-1' REAL, 'PrEC' REAL, 'ProgFib' REAL, 'RPTEC' REAL, 'SAEC' REAL, 'SK-N-SH_RA' REAL, 'SKMC' REAL, 'SKNMC' REAL, 'Stellate' REAL, 'Th1' REAL, 'Th2' REAL, 'Urothelia_UT189' REAL, 'WI-38' REAL, 'WI-38-TAM' REAL)" )

    with open(exp_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        cur.executemany("INSERT INTO exp112 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);", (row,))
    conn.commit()
    conn.close()

def create_dhsClusterTB(dhsCluster_file):
    conn = sqlite3.connect('dhs.db')
    conn.text_factory = str
    cur = conn.cursor()

    print("\t\t Creating DHS-to-Cluster table...")

    cur.execute("DROP TABLE IF EXISTS dhsCluster")
    cur.execute("CREATE TABLE dhsCluster ('chr' TEXT, 'start' INTEGER, 'end' INTEGER, 'original_cluster' INTEGER, 'refined_cluster' INTEGER, 'original_distance' REAL, 'refined_distance' REAL)")

    with open(dhsCluster_file, 'rb') as cfile:
        data = [row for row in csv.reader(cfile.read().splitlines())]
    for row in data:
        row = row[0].split()
        
        #NOTE: dhs_id is (rowID-1)
        cur.executemany("INSERT INTO dhsCluster VALUES (?,?,?,?,?,?,?);", (row,))
    conn.commit()
    conn.close()

    



def create_dhsDB():
    
    create_corrTB(correlation_input)
    create_overlapTB(overlap_input)
    create_concordanceTB(concordance_input)
    create_openCellTypesTB(openCellTypes_input)
    create_dhsPredictorsTB(dhsPredictors_input)
    create_openSamplesTB(openSamples_input)
    create_motifJasparTB(motifJaspar_input)
    create_malignantTissuesTB(malignantTissues_input)
    create_sexTB(sex_input)
    create_dhsPredictorsTB(dhsPredictors_input)
    create_dhs112TB(dhs112_input)
    create_exp112TB(exp112_input)
    create_dhsClusterTB(dhsCluster_input)


#------------------------------ END OF DATABASE POPULATION -------------------------------------------

if __name__ == "__main__":
    correlation_input = '../dhs_data/allGeneCorrelations100000.p05_v3.txt'
    overlap_input = '../dhs_data/TableS05-overlapSummary.txt'
    concordance_input = '../dhs_data/TableS01-concordanceSummary.csv'
    openCellTypes_input = '../dhs_data/TableS04-cluster-to-openCellTypes_v3.txt'
    openSamples_input = '../dhs_data/TableS02-openSamples_v3.txt'
    motifJaspar_input = '../dhs_data/TableS09-motifJaspar.txt'
    malignantTissues_input = '../dhs_data/TableS06-results_tissueMalig_v2.txt'
    sex_input = '../dhs_data/TableS08-results_sex_v3.txt'
    dhsPredictors_input = '../dhs_data/TableS07-dhsPredictors_v2.bed' 
    dhs112_input = '../dhs_data/dhs112_v3.bed'
    exp112_input = '../dhs_data/exp112.bed'
    dhsCluster_input = '../dhs_data/TableS03-dhs-to-cluster.txt'

    create_dhsDB()
    
    # TODO: Build index for tables.

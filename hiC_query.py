#!/usr/bin/python
import argparse
import os
import sqlite3
import pybedtools
from sets import Set

def build_snp_index(snp_database_fp,snp_dir):
	print "Building SNP index..."
	snp_index_db = sqlite3.connect(snp_database_fp)
	snp_index = snp_index_db.cursor()
	snp_index.execute("CREATE TABLE snps (rsID text, chr text, locus integer)")
	snp_index.execute("CREATE INDEX id ON snps (rsID,chr,locus)")
	for bed_fp in os.listdir(snp_dir):
		print "\tProcessing " + bed_fp
		bed = open(snp_dir + '/' + bed_fp,'r')
		for line in bed:
			if line.startswith("track name"):
				continue
			snp = line.strip().split('\t')
			snp_index.execute("INSERT INTO snps VALUES(?,?,?)",[snp[3],snp[0][snp[0].find("chr")+3:],int(snp[1])])
		bed.close()
	print "\tWriting SNP index to file..."
	snp_index_db.commit()
	print "Done building SNP index."

def process_inputs(inputs,snp_database_fp,snp_dir):
	snps = {}
	if not os.path.isfile(snp_database_fp): #Need to build SNP index
		build_snp_index(snp_database_fp,snp_dir)
	print "Loading SNP index..."
	snp_index_db = sqlite3.connect(snp_database_fp)
	snp_index_db.text_factory = str
	snp_index = snp_index_db.cursor()
	print "Processing input..."
	for input in inputs:
		snp = None
		if input.startswith("rs"):
			snp_index.execute("SELECT * FROM snps WHERE rsID=?",(input,))
		else:
			chr = input[input.find("chr")+3:input.find(':')]
			locus = int(input[input.find(':')+1:])
			snp_index.execute("SELECT * FROM snps WHERE chr=? and locus=?",[chr,locus])
		snp = snp_index.fetchone()
		if snp == None:
			print "Warning: %s does not exist in SNP database." % input
		else:
			snps[snp[0]]=(snp[1],snp[2])
	return snps

def find_interactions(snps,hic_data_dir,distance,include,exclude):
	print "Finding interactions..."
	if include:
		include_set = Set(include)
	if exclude:
		exclude_set = Set(exclude)
	interactions = {}
	for snp in snps.keys():
		interactions[snp] = []
	for cell_line in os.listdir(hic_data_dir):
		if (include and not cell_line in include) or (exclude and cell_line in exclude):
			continue
		if os.path.isdir(hic_data_dir + '/' + cell_line):
			print "\tSearching cell line " + cell_line
			for replicate in os.listdir(hic_data_dir + '/' + cell_line):
				if replicate.endswith(".db"):
					db = hic_data_dir + '/' + cell_line + '/' + replicate
					print "\t\tSearching replicate " + replicate
					table_index_db = sqlite3.connect(db)
					table_index_db.text_factory = str
					table_index = table_index_db.cursor()
					for snp in snps.keys():
						print "\t\t\tFinding interactions for " + str(snp)
						snp_chr = snps[snp][0]
						snp_locus = snps[snp][1]
						for interaction in table_index.execute("SELECT chr1, fragment1 FROM chr%s_interactions WHERE chr2=? and (locus2 >= ? and locus2 <= ?)" % snp_chr,[snp_chr,snp_locus-distance,snp_locus+distance]):
							interactions[snp].append(interaction)
						for interaction in table_index.execute("SELECT chr2, fragment2 FROM chr%s_interactions WHERE chr1=? and (locus1 >= ? and locus1 <= ?)" % snp_chr,[snp_chr,snp_locus-distance,snp_locus+distance]):
							interactions[snp].append(interaction)
	return interactions

def find_genes(interactions,fragment_bed_fp,gene_bed_fp):
	#For each SNP
		#Generate BED file of all fragments interacting with SNP-containing fragment
		#Get intersection of this BED file with BED file detailing gene locations
		#Return a list of genes with which SNP is interacting

def find_eqtls(interactions,eqtl_data_dir):
	print "Identifying eQTLs of interest..."
	eqtls = {} #A mapping of SNPs to different tissue types, which are in turn mapped to any eQTLs that coincide with a spatial interaction with said SNPs
	for snp in interactions.keys():
		eqtls[snp] = {}
	for db in os.listdir(eqtl_data_dir): #Iterate through databases of eQTLs by tissue type
		print "\tQuerying " + db[:db.rfind('.')]
		eqtl_index_db = sqlite3.connect(eqtl_data_dir + '/' + db)
		eqtl_index_db.text_factory = str
		eqtl_index = eqtl_index_db.cursor()
		tissue = db[:db.rfind('.')]
		for snp in interactions.keys():
			eqtls[snp][tissue] = Set([]) #Create a list of relevant eQTLs for this tissue type
			for eqtl in eqtl_index.execute("SELECT gene_name, gene_chr, gene_start, gene_end FROM eqtls WHERE rsID=?",(snp,)): #Pull down all eQTLs related to a given SNP to test for relevance:
				for interaction in interactions[snp]:
					if eqtl[1] == interaction[0] and (eqtl[2] <= interaction[1] and eqtl[3] >= interaction[1]): #If an eQTL coincides with a spatial interaction, save it
						eqtls[snp][tissue].add(eqtl)
	return eqtls


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--inputs",nargs='+',required=True,help="The the dbSNP IDs or loci of SNPs of interest in the format \"chr<x>:<locus>\"")
	parser.add_argument("-n","--include_cell_lines",nargs='+',help="Space-separated list of cell lines to include (others will be ignored). NOTE: Mutually exclusive with EXCLUDE_CELL_LINES.")
	parser.add_argument("-x","--exclude_cell_lines",nargs='+',help="Space-separated list of cell lines to exclude (others will be included). NOTE: Mutually exclusive with INCLUDE_CELL_LINES.")
	parser.add_argument("-d","--distance",type=int,default=500,help="The allowed distance from the locus of interest for a fragment to be considered.")
	parser.add_argument("-s","--snp_database_fp",default="snpIndex.db",help="The database of SNPs to search for details on input dbSNP IDs (Optional - will try to build from files in SNP_DIR if non-existent.")
	parser.add_argument("-t","--snp_dir",default="snps",help="The directory containing (BED) files of SNPs to use as a reference (the dbSNP Build 144 (GRCh37.p13) BED files by default).")
	parser.add_argument("-c","--hic_data_dir",default="hic_data",help="The directory containing the directories representing cell lines of HiC experiment tables.")
	parser.add_argument("-f","--fragment_bed_fp",default="Homo_sapiens.ensembl.release-74.MboI.fragments.bed",help="Bed file detailing the start and end points of HiC experiment fragments (hg19 fragmented by MboI by default.)")
	parser.add_argument("-g","--gene_bed_fp",default="hg19_genes.bed",help="Bed file detailing locations of genes in genome (the UCSC known gene list for hg19 by default).")
	parser.add_argument("-e","--eqtl_data_dir",default="eQTLs",help="The directory containing databases of eQTL data from various tissues (the GTEx Analysis V6 by default).")
	args = parser.parse_args()
	snps = process_inputs(args.inputs,args.snp_database_fp,args.snp_dir)
	interactions = find_interactions(snps,args.hic_data_dir,args.distance,args.include_cell_lines,args.exclude_cell_lines)
	genes = find_genes(interactions,args.fragment_bed_fp,args.gene_bed_fp)
	eqtls = find_eqtls(interactions,args.eqtl_data_dir)
	"""for snp in interactions.keys():
		print snp + ':'
		for interaction in interactions[snp]:
			print "\t" + str(interaction)"""
	for snp in interactions.keys():
		print snp + ':'
		for tissue in eqtls[snp].keys():
			print '\t' + tissue + ':'
			for eqtl in eqtls[snp][tissue]:
				print '\t\t' + str(eqtl)

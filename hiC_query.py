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
		if os.path.isfile(input):
			with open(input,'r') as infile:
				for line in infile:
					id = line.strip().split(' ')[0]
					snp = None
					if id.startswith("rs"):
						snp_index.execute("SELECT * FROM snps WHERE rsID=?",(id,))
					else:
						chr = id[id.find("chr")+3:id.find(':')]
						locus = int(id[id.find(':')+1:])
						snp_index.execute("SELECT * FROM snps WHERE chr=? and locus=?",[chr,locus])
					snp = snp_index.fetchone()
					if snp == None:
						print "Warning: %s does not exist in SNP database." % id
					else:
						snps[snp[0]]=(snp[1],snp[2])
		else:
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

def find_interactions(snps,fragment_database_fp,hic_data_dir,distance,include,exclude,cis_interactions_only):
	print "Finding interactions..."
	#Query fragmentIndex to find to which fragment the SNP belongs
	fragment_index_db = sqlite3.connect(fragment_database_fp)
	fragment_index_db.text_factory = str
	fragment_index = fragment_index_db.cursor()
	#Look for all interactions involving said fragment in the HiC databases
	if include:
		include_set = Set(include)
	if exclude:
		exclude_set = Set(exclude)
	interactions = {} #A mapping of each SNP to the fragments with which the fragment it is on interacts
	for snp in snps.keys():
		interactions[snp] = Set([])
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
						fragment_index.execute("SELECT fragment FROM fragments WHERE chr=? AND start<=? AND end>=?",["chr" + snps[snp][0],snps[snp][1],snps[snp][1]])
						snp_fragment_result = fragment_index.fetchone()
						if snp_fragment_result == None:
							print "Warning: error retrieving SNP fragment for SNP " + snp
							continue
						snp_fragment = snp_fragment_result[0]
						snp_chr = snps[snp][0]
						if(cis_interactions_only):
							for interaction in table_index.execute("SELECT chr1, fragment1 FROM chr%s_interactions WHERE (chr1=? and chr2=?) and fragment2=?" % snp_chr,[snp_chr,snp_chr,snp_fragment]):
								interactions[snp].add(interaction)
							for interaction in table_index.execute("SELECT chr2, fragment2 FROM chr%s_interactions WHERE (chr1=? and chr2=?) and fragment1=?" % snp_chr,[snp_chr,snp_chr,snp_fragment]):
								interactions[snp].add(interaction)
						else:
							for interaction in table_index.execute("SELECT chr1, fragment1 FROM chr%s_interactions WHERE chr2=? and fragment2=?" % snp_chr,[snp_chr,snp_fragment]):
								interactions[snp].add(interaction)
							for interaction in table_index.execute("SELECT chr2, fragment2 FROM chr%s_interactions WHERE chr1=? and fragment1=?" % snp_chr,[snp_chr,snp_fragment]):
								interactions[snp].add(interaction)
	return interactions

def find_genes(interactions,fragment_database_fp,gene_bed_fp):
	fragment_index_db = sqlite3.connect(fragment_database_fp)
	fragment_index_db.text_factory = str
	fragment_index = fragment_index_db.cursor()
	hs_gene_bed = pybedtools.BedTool(gene_bed_fp)
	genes = {}
	for snp in interactions.keys():
		#Generate BED file of all fragments interacting with SNP-containing fragment
		genes[snp] = Set([])
		temp_snp_bed = open("temp_snp_bed.bed",'w')
		for interaction in interactions[snp]:
			fragment_index.execute("SELECT start, end FROM fragments WHERE chr=? and fragment=?",["chr" + interaction[0],interaction[1]])
			fragment_pos = fragment_index.fetchone()
			if fragment_pos == None:
				print "Warning: error retrieving fragment %s on chromosome %s" % (interaction[1],interaction[0])
			temp_snp_bed.write("%s\t%s\t%s\n" % ("chr" + interaction[0],fragment_pos[0],fragment_pos[1]))
		temp_snp_bed.close()
		int_bed = pybedtools.BedTool("temp_snp_bed.bed")
		#Get intersection of this BED file with BED file detailing gene locations
		gene_bed = hs_gene_bed.intersect(int_bed,u=True)
		#Return a list of genes with which SNP is interacting
		for feat in gene_bed:
			genes[snp].add(str(feat.name))
	return genes

def find_eqtls(genes,eqtl_data_dir):
	print "Identifying eQTLs of interest..."
	eqtls = {} #A mapping of SNPs to genes with which they have an eQTL relationship, which in turn maps to a list of tissues in which this eQTL occurs
	for snp in genes.keys():
		eqtls[snp] = {} #Stores eQTLs relevant to this SNP, and, for each eQTL, a list of the tissues in which it is found
	for db in os.listdir(eqtl_data_dir): #Iterate through databases of eQTLs by tissue type
		tissue = db[:db.rfind('.')]
		print "\tQuerying " + tissue
		eqtl_index_db = sqlite3.connect(eqtl_data_dir + '/' + db)
		eqtl_index_db.text_factory = str
		eqtl_index = eqtl_index_db.cursor()
		for snp in genes.keys():
			for gene in genes[snp]:
				for eqtl in eqtl_index.execute("SELECT gene_name, ens_id, gene_chr, gene_start, gene_end FROM eqtls WHERE rsID=? AND ens_id=?",(snp,gene)): #Pull down all eQTLs related to a given SNP to test for relevance:
					if not eqtls[snp].has_key(eqtl):
						eqtls[snp][eqtl] = Set([])
					eqtls[snp][eqtl].add(tissue)
	for snp in eqtls.keys():
		print snp + ':'
		for eqtl in eqtls[snp].keys():
			print '\t' + str(eqtl)
			print '\t\t' + str(eqtls[snp][eqtl])
	return eqtls

def produce_output(snps,interactions,eqtls,output_dir):
	print "Producing output..."
	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)
	summary_table = open(output_dir + "/summary_table.txt",'w')
	summary_table.write("SNP\tchr\tpos\tnum_contacts\tnum_eQTLs\n")
	for snp in snps.keys():
		summary_table.write(snp + '\t' + snps[snp][0] + '\t' + str(snps[snp][1]) + '\t' + str(len(interactions[snp])) + '\t' + str(len(eqtls[snp])) + '\n')
		snp_summary = open(output_dir + '/' + snp + ".txt",'w')
		snp_summary.write("gene_name\tensembl_id\tchr\tstart_pos\tend_pos\tnum_tissues\ttissues\tdistance_from_snp\n")
		for eqtl in eqtls[snp].keys():
			distance_from_snp = 0
			if(not snps[snp][0] == eqtl[2]):
				distance_from_snp = None #Not applicable to trans interactions
			elif(snps[snp][1] < eqtl[3]):
				distance_from_snp = eqtl[3] - snps[snp][1]
			elif(snps[snp][1] > eqtl[4]):
				distance_from_snp = snps[snp][1] - eqtl[4]
			snp_summary.write(eqtl[1] + '\t' + eqtl[0] + '\t' +eqtl[2] + '\t' + str(eqtl[3]) + '\t' + str(eqtl[4]) + '\t' + str(len(eqtls[snp][eqtl])) + '\t') 
			for tissue in eqtls[snp][eqtl]:
				snp_summary.write(tissue + ',') 
			snp_summary.write('\t' + str(distance_from_snp) + '\n')
		snp_summary.close()
	summary_table.close()

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
	parser.add_argument("-b","--fragment_database_fp",default="fragmentIndex.db",help="The database of fragments to search when assigning SNPs to fragments (GRCh37 digested by MboI by default.)")
	parser.add_argument("-g","--gene_bed_fp",default="hg19_genes.bed",help="Bed file detailing locations of genes in genome (the UCSC known gene list for hg19 by default).")
	parser.add_argument("-e","--eqtl_data_dir",default="eQTLs",help="The directory containing databases of eQTL data from various tissues (the GTEx Analysis V6 by default).")
	parser.add_argument("-o","--output_dir",default="hiCquery_output",help="The directory in which to output results (\"hiCquery_output\" by default).")
	parser.add_argument("-j","--cis_interactions_only",action="store_true",default=False,help="Consider only cis interactions, i.e. interactions occurring on the same chromosome.")
	args = parser.parse_args()
	snps = process_inputs(args.inputs,args.snp_database_fp,args.snp_dir)
	interactions = find_interactions(snps,args.fragment_database_fp,args.hic_data_dir,args.distance,args.include_cell_lines,args.exclude_cell_lines,args.cis_interactions_only)
	genes = find_genes(interactions,args.fragment_database_fp,args.gene_bed_fp)
	with open("test_gene_output.txt",'w') as test_gene_output:
		test_gene_output.write(str(genes))
	eqtls = find_eqtls(genes,args.eqtl_data_dir)
	produce_output(snps,interactions,eqtls,args.output_dir)

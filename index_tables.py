#!/usr/bin/python
import argparse
import sqlite3
import os
from sets import Set

def build_hic_table_index(table_fp,mapq_cutoff):
	chr_set = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"]
	if os.path.isfile(table_fp[:table_fp.rfind('.')] + ".db"):
		os.remove(table_fp[:table_fp.rfind('.')] + ".db")
	
	table_index_db = sqlite3.connect(table_fp[:table_fp.rfind('.')] + ".db")
	table_index = table_index_db.cursor()
	for chr in chr_set:
		table_index.execute("CREATE TABLE chr%s_interactions (chr1 text, locus1 integer, fragment1 integer, chr2 text, locus2 integer, fragment2 integer)" % chr)
		table_index.execute("CREATE INDEX chr%s_loci ON chr%s_interactions (chr1,locus1,chr2,locus2)" % (chr,chr))
	##Do line count for progress meter
	do_linecount = True
	print "Determining table size..."
	eqtl_table = open(table_fp,'r')
	lines = 0
	for i in eqtl_table:
		lines += 1
	eqtl_table.close()
	lines = lines//100*100 #Get an approximation
	do_linecount = not lines == 0
	
	open_dbs = []
	with open(table_fp,'r') as rao_table:
		print "Indexing HiC interaction table..."
		for i,line in enumerate(rao_table):
			if do_linecount:
				if i % (lines/100) == 0:
					print "\tProcessed %d%%..." % ((float(i)/float(lines))*100)
			interaction = line.strip().split(' ')
			chr1 = interaction[2]
			chr2 = interaction[6]
			if int(interaction[9]) >= mapq_cutoff and int(interaction[10]) >= mapq_cutoff:
				try:
					table_index.execute("INSERT INTO chr%s_interactions VALUES (?,?,?,?,?,?)" % chr1,[chr1,int(interaction[3]),int(interaction[4]),chr2,int(interaction[7]),int(interaction[8])])
					if not chr1 == chr2:
						table_index.execute("INSERT INTO chr%s_interactions VALUES (?,?,?,?,?,?)" % chr2,[chr1,int(interaction[3]),int(interaction[4]),chr2,int(interaction[7]),int(interaction[8])])
				except sqlite3.OperationalError:
					print "Warning: error inserting line \"%s\" into database." % line.strip()
					continue
	table_index_db.commit()
	print "Done indexing HiC interaction table."

def build_eqtl_table_index(table_fp):
	if os.path.isfile(table_fp[:table_fp.rfind('.')] + ".db"):
		os.remove(table_fp[:table_fp.rfind('.')] + ".db")
	table_index_db = sqlite3.connect(table_fp[:table_fp.rfind('.')] + ".db")
	table_index = table_index_db.cursor()
	table_index.execute("CREATE TABLE eqtls (rsID text, gene_name text, ens_id text, gene_chr text, gene_start integer, gene_end integer)")
	table_index.execute("CREATE INDEX id ON eqtls (rsID)")
	##Do line count for progress meter
	do_linecount = True
	print "Determining table size..."
	eqtl_table = open(table_fp,'r')
	lines = 0
	for i in eqtl_table:
		lines += 1
	eqtl_table.close()
	lines = lines//100*100 #Get an approximation
	do_linecount = not lines == 0
	
	eqtl_table = open(table_fp,'r')
	for i,line in enumerate(eqtl_table):
		if do_linecount:
			if i % (lines/100) == 0:
				print "\tProcessed %d%%..." % ((float(i)/float(lines))*100)
		if i == 0:
			continue
		eqtl = line.strip().split('\t')
		table_index.execute("INSERT INTO eqtls VALUES (?,?,?,?,?,?)",[eqtl[22],eqtl[26],eqtl[1],eqtl[28],eqtl[29],eqtl[30]])
	eqtl_table.close()
	table_index_db.commit()
	print "Done indexing eQTL table."

def build_fragment_index(table_fp):
	if os.path.isfile(table_fp[:table_fp.rfind('.')] + ".db"):
		os.remove(table_fp[:table_fp.rfind('.')] + ".db")
	fragment_index_fp = sqlite3.connect(table_fp[:table_fp.rfind('.')] + ".db")
	fragment_index = fragment_index_fp.cursor()
	fragment_index.execute("CREATE TABLE fragments (chr text, start integer, end integer, fragment integer)")
	fragment_index.execute("CREATE INDEX f_index ON fragments (chr,fragment)")
	
	fragments_bed = open(table_fp,'r')
	for line in fragments_bed:
		fragment = line.strip().split('\t')
		fragment_index.execute("INSERT INTO fragments VALUES (?,?,?,?)", [fragment[0][fragment[0].find("chr")+3:],fragment[1],fragment[2],fragment[3]])
	fragment_index.commit()
	fragments_bed.close()

def build_snp_index(snp_dir):
	print "Building SNP index..."
	if not os.path.isdir(snp_dir):
		print "Error: argument to build SNP index must be a directory."
		break
	if os.path.isfile(snp_dir + "/snpIndex.db"):
		os.remove(snp_dir + "/snpIndex.db")
	snp_index_db = sqlite3.connect(snp_dir + "/snpIndex.db")
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


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Indexes HiC and SNP tables into searchable databases. NOTE: must supply an argument to INPUT_TYPE")
	parser.add_argument("-t","--table")
	parser.add_argument("-i","--input_type",required=True,help="The input type of the table (valid options are 'h', for HiC, 'e', for eqtl, 'f', for fragmented genome BED, or 's', for SNP. In the case of SNP, a directory is the expected argument).")
	parser.add_argument("-m","--mapq_cutoff",type=int,default=150,help="Mapq cutoff used to index HiC tables, not used for eQTL tables.")
	args = parser.parse_args()
	opts = Set(['i','e','f','s'])
	if not args.input_type in opts:
		raise InputError("A valid input type must be specified.")
	if args.input_type == 'h':
		build_hic_table_index(args.table,args.mapq_cutoff)
	elif args.input_type == 'e':
		build_eqtl_table_index(args.table)
	elif args.input_type == 'f':
		build_fragment_index(args.table)
	elif args.input_type == 's':
		build_snp_index(args.table)

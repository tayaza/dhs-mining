#!/usr/bin/python
import argparse
import os
import sqlite3
from sets import Set
import time

generate_circos_output = False

def build_snp_index():
	print "Building SNP index..."
	snp_index_db = sqlite3.connect("snpIndex.test.db")
	snp_index = snp_index_db.cursor()
	snp_index.execute("CREATE TABLE snps (rsID text, chr text, locus integer)")
	snp_index.execute("CREATE UNIQUE INDEX id ON snps (rsID)")
	for bed_fp in os.listdir("snps"):
		print "\tProcessing " + bed_fp
		bed = open("snps/" + bed_fp,'r')
		for line in bed:
			if line.startswith("track name"):
				continue
			snp = line.strip().split('\t')
			snp_index.execute("INSERT INTO snps VALUES(?,?,?)",[snp[3],snp[0][snp[0].find("chr")+3:],int(snp[1])])
		bed.close()
	print "\tWriting SNP index to file..."
	snp_index_db.commit()
	print "Done building SNP index."

def process_inputs(inputs):
	loci = []
	
	for input in inputs:
		if input.startswith("rs"):
			if not os.path.isfile("snpIndex.test.db"): #Need to build SNP index
				build_snp_index()
			print "Loading SNP index..."
			snp_index_db = sqlite3.connect("snpIndex.test.db")
			snp_index_db.text_factory = str
			snp_index = snp_index_db.cursor()
			break
	print "Processing input..."
	for input in inputs:
		#print "\tProcessing " + input
		if input.startswith("rs"):
			snp_index.execute("SELECT chr, locus FROM snps WHERE rsID=?",(input,))
			locus = snp_index.fetchone()
			if locus == None:
				print "Error: %s does not exist in database." % input
			else:
				loci.append(locus)
		else:
			loci.append((input[input.find("chr")+3:input.find(':')],int(input[input.find(':')+1:])))
	return loci

"""def build_table_index(table_fp,mapq_cutoff):
	print "Indexing HiC interaction tables..."
	table_index_db = sqlite3.connect("tableIndex.test.db")
	table_index = table_index_db.cursor()
	table_index.execute("CREATE TABLE interactions (readName text, chr1 text, locus1 integer, chr2 text, locus2 integer)")
	table_index.execute("CREATE INDEX loci ON interactions (chr1,locus1,chr2,locus2)")
	rao_table = open(table_fp,'r')
	for line in rao_table:
		interaction = line.strip().split(' ')
		if int(interaction[9]) >= mapq_cutoff and int(interaction[10]) >= mapq_cutoff:
			table_index.execute("INSERT INTO interactions VALUES (?,?,?,?,?)",[interaction[0],interaction[2],int(interaction[3]),interaction[6],int(interaction[7])])
	rao_table.close()
	print "Writing table index to file..."
	table_index_db.commit()
	print "Done indexing HiC interaction tables."""

"""def find_interactions_no_index(loci,table_fp,distance): ##Read database file
	interactions = {}
	#if not os.path.isfile("tableIndex.test.db"): #Need to index Rao table
		#build_table_index(table_fp,mapq_cutoff)
	table_index_db = sqlite3.connect("HIC003_mapq150.db")
	table_index_db.text_factory = str
	table_index = table_index_db.cursor()
	for locus in loci:
		print "\tFinding interactions for " + str(locus)
		interactions[locus] = []
		for interaction in table_index.execute("SELECT chr1, locus1 FROM interactions WHERE chr2=? and (locus2 >= ? and locus2 <= ?)",[locus[0],locus[1]-distance,locus[1]+distance]):
			interactions[locus].append(interaction)
		for interaction in table_index.execute("SELECT chr2, locus2 FROM interactions WHERE chr1=? and (locus1 >= ? and locus1 <= ?)",[locus[0],locus[1]-distance,locus[1]+distance]):
			interactions[locus].append(interaction)
	return interactions"""

def find_interactions(loci,db_dir,distance): ##Read database file
	interactions = {}
	open_dbs = []
	db_dir = db_dir.strip('/') #Strip trailing '/' if it is there
	#if not os.path.isfile("tableIndex.test.db"): #Need to index Rao table
		#build_table_index(table_fp,mapq_cutoff)
	table_index_db = sqlite3.connect(db_dir + "/dummy.db")
	table_index_db.text_factory = str
	table_index = table_index_db.cursor()
	for locus in loci:
		print "\tFinding interactions for " + str(locus)
		if locus[0] not in open_dbs:
			if len(open_dbs) >= 10:
				table_index.execute("DETACH db" + open_dbs.pop())
			table_index.execute("ATTACH '" + db_dir + '/' + db_dir[db_dir.rfind('/'):db_dir.rfind('_')] + '_' + locus[0] + ".db' AS db" + locus[0])
		else:
			open_dbs.remove(locus[0])
		open_dbs.insert(0,locus[0])
		interactions[locus] = []
		for interaction in table_index.execute("SELECT chr1, locus1 FROM db" + locus[0] + ".interactions WHERE chr2=? and (locus2 >= ? and locus2 <= ?)",[locus[0],locus[1]-distance,locus[1]+distance]):
			interactions[locus].append(interaction)
		for interaction in table_index.execute("SELECT chr2, locus2 FROM db" + locus[0] + ".interactions WHERE chr1=? and (locus1 >= ? and locus1 <= ?)",[locus[0],locus[1]-distance,locus[1]+distance]):
			interactions[locus].append(interaction)
	outfile = open("db_index_output_test.txt",'w')
	for locus in interactions.keys():
		for interaction in interactions[locus]:
			outfile.write(locus[0] + '\t' + str(locus[1]) + '\t' + interaction[0] + '\t' + str(interaction[1]) + '\n')
	os.remove(db_dir + "/dummy.db")
	return interactions

"""def find_interactions_raw(loci,table_fp,distance,mapq_cutoff): ##Read raw HiC file, no indexing
	interactions = {}
	for locus in loci:
		print "\tFinding interactions for " + str(locus)
		rao_table = open(table_fp,'r')
		chr = locus[0]
		coord = locus[1]
		interactions[locus] = []
		for line in rao_table:
			interaction = line.strip().split(' ')
			chr1 = interaction[2]
			chr2 = interaction[6]
			locus1 = int(interaction[3])
			locus2 = int(interaction[7])
			if ((chr1 == chr and (locus1-distance < coord and locus1+distance > coord)) or (chr2 == chr and (locus2-distance < coord and locus2+distance > coord))):
				if chr1 == chr and (locus1-distance < coord and locus1+distance > coord):
					interactions[locus].append((chr2,locus2,line.strip()))
				if chr2 == chr and (locus2-distance < coord and locus2+distance > coord):
					interactions[locus].append((chr1,locus1,line.strip()))
		rao_table.close()
	outfile = open("raw_output_test.txt",'w')
	for locus in interactions.keys():
		for interaction in interactions[locus]:
			outfile.write(locus[0] + '\t' + str(locus[1]) + '\t' + interaction[0] + '\t' + str(interaction[1]) + '\t' + interaction[2] + '\n')
	return interactions"""


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--inputs",nargs='+',help="The the dbSNP ID or loci of interest in the format \"chr<x>:<locus>\"")
	parser.add_argument("-t","--table_fp",default="HIC059_mapq150_chr21.test.txt",help="TEMP: The Rao table from which to pull connections.")
	parser.add_argument("-x","--db_dir",default="GSM1551552_HIC003_merged_nodups_db",help="TEMP: The directory containing the chromosome-indexed database to search.")
	parser.add_argument("-d","--distance",type=int,default=500,help="The allowed distance from the locus of interest for a fragment to be considered.")
	parser.add_argument("-m","--mapq_cutoff",type=int,default=150,help="The minimum mapq score allowed.")
	parser.add_argument("-c","--generate_circos_output",action="store_true",default=False,help="Produce Circos-compatible output files (default: False)")
	args = parser.parse_args()
	inputs = args.inputs
	table_fp = args.table_fp
	db_dir = args.db_dir
	distance = args.distance
	#mapq_cutoff = args.mapq_cutoff
	loci = process_inputs(inputs)
	interactions = find_interactions(loci,db_dir,distance)

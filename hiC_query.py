#!/usr/bin/python
import argparse
import os
import sqlite3
from sets import Set

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
	for input in inputs:
		print "Processing " + input
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
	
def find_interactions(loci,table_fp,distance,mapq_cutoff):
	interactions = {}
	for locus in loci:
		rao_table = open(table_fp,'r')
		chr = locus[0]
		coord = int(locus[1])
		interactions[locus] = []
		print "Processing " + chr + ' ' + str(coord)
		for line in rao_table:
			interaction = line.strip().split(' ')
			chr1 = interaction[2]
			chr2 = interaction[6]
			coord1 = int(interaction[3])
			coord2 = int(interaction[7])
			print '\t' + chr1 + ' ' + str(coord1) + '\n\t' + chr2 + ' ' + str(coord2)
			if ((chr1 == chr and (coord1-distance < coord and coord1+distance > coord)) or (chr2 == chr and (coord2-distance < coord and coord2+distance > coord))) and (int(interaction[9]) > mapq_cutoff and int(interaction[10]) > mapq_cutoff):
				if chr1 == chr and (coord1-distance < coord and coord1+distance > coord):
					interactions[locus].append((chr2,coord2))
				else:
					interactions[locus].append((chr1,coord1))
		rao_table.close()
	return interactions


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--inputs",nargs='+',help="The the dbSNP ID or loci of interest in the format \"chr<x>:<locus>\"")
	parser.add_argument("-t","--table_fp",default="HIC059_mapq150_chr21.test.txt",help="TEMP: The Rao table from which to pull connections...")
	parser.add_argument("-d","--distance",type=int,default=500,help="The allowed distance from the locus of interest for a fragment to be considered...")
	parser.add_argument("-m","--mapq_cutoff",type=int,default=150,help="The minimum mapq score allowed...")
	args = parser.parse_args()
	inputs = args.inputs
	table_fp = args.table_fp
	distance = args.distance
	mapq_cutoff = args.mapq_cutoff
	loci = process_inputs(inputs)
	print loci
	interactions = find_interactions(loci,table_fp,distance,mapq_cutoff)
	print interactions
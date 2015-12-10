#!/usr/bin/python
import argparse
import os
import sqlite3
from sets import Set

def build_snp_index():
	print "Building SNP index..."
	snp_index_db = sqlite3.connect("snpIndex.db")
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
			if not os.path.isfile("snpIndex.db"): #Need to build SNP index
				build_snp_index()
			print "Loading SNP index..."
			snp_index_db = sqlite3.connect("snpIndex.db")
			snp_index_db.text_factory = str
			snp_index = snp_index_db.cursor()
			break
	for input in inputs:
		print "Processing " + input
		if input.startswith("rs"):
			snp_index.execute("SELECT chr, locus FROM snps WHERE rsID=?",(input,))
			loci.append(snp_index.fetchone())
		else:
			loci.append((input[input.find("chr")+3:input.find(':')],int(input[input.find(':')+1:])))
	return loci
	
def find_interactions(args):
	chr = args.locus[args.locus.find("chr")+3:args.locus.find(':')]
	locus = int(args.locus[args.locus.find(':')+1:])
	rao_table = open(args.table_fp,'r')
	interactions = []
	for line in rao_table:
		interaction = line.strip().split(' ')
		chr1 = interaction[2]
		chr2 = interaction[6]
		locus1 = int(interaction[3])
		locus2 = int(interaction[7])
		if ((chr1 == chr and (locus1-args.distance < locus and locus1+args.distance > locus)) or (chr2 == chr and (locus2-args.distance < locus and locus2+args.distance > locus))) and (int(interaction[9]) > args.mapq_cutoff and int(interaction[10]) > args.mapq_cutoff):
			if chr1 == chr and (locus1-args.distance < locus and locus1+args.distance > locus):
				interactions.append(Interaction(chr2,locus2))
			else:
				interactions.append(Interaction(chr1,locus1))
		return interactions


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--inputs",nargs='+',help="The the dbSNP ID or loci of interest in the format \"chr<x>:<locus>\"")
	parser.add_argument("-t","--table_fp",default="HIC059_mapq150_chr21.txt",help="TEMP: The Rao table from which to pull connections...")
	parser.add_argument("-d","--distance",type=int,default=500,help="The allowed distance from the locus of interest for a fragment to be considered...")
	parser.add_argument("-m","--mapq_cutoff",type=int,default=150,help="The minimum mapq score allowed...")
	args = parser.parse_args()
	inputs = args.inputs
	table_fp = args.table_fp
	distance = args.distance
	mapq_cutoff = args.mapq_cutoff
	loci = process_inputs(inputs)
	print inputs
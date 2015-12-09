#!/usr/bin/python
import argparse
import os
import sqlite3
from sets import Set

def build_snp_index():
	print "Building SNP index..."
	snpIDs = Set([]) #Set of all dbSNP IDs, for quick validation
	snp_index = {} #Links all dbSNP IDSs to a chromosome and a locus
	for bed_fp in os.listdir("snps"):
		print "\tProcessing " + bed_fp
		bed = open("snps/" + bed_fp,'r')
		for line in bed:
			if line.startswith("track name"):
				continue
			snp = line.strip().split('\t')
			snpIDs.add(snp[3])
			snp_index[snp[3]] = (snp[0][snp[0].find("chr")+3:],int(snp[1]))
		bed.close()
	print "\tWriting ID set to file..."
	snpIDs_file = open("snpIDs.pkl",'wb')
	pickle.dump(snpIDs,snpIDs_file)
	snpIDs_file.close()
	print "\tWriting index to file..."
	snp_index_file = open("snpIndex.pkl",'wb')
	pickle.dump(snp_index,snp_index_file)
	snp_index_file.close()
	print "Done building SNP index."

def process_inputs(loci):
	inputs = []
	for locus in loci:
		if locus.startswith("rs"):
			if not os.path.isfile("snpIndex.pkl") or not os.path.isfile("snpIDs.pkl"): #Need to build SNP index
				build_snp_index()
			print "Loading SNP IDs..."
			snpIDs_file = open("snpIDs.pkl",'rb')
			snpIDs = pickle.load(snpIDs_file)
			snpIDs_file.close()
			print "Loading SNP index..."
			snp_index_file = open("snpIndex.pkl",'rb')
			snp_index = pickle.load(snp_index_file)
			snp_index_file.close()
			break
	for locus in loci:
		print "Processing " + locus
		if locus.startswith("rs"):
			if locus in snpIDs:
				inputs.append(snp_index[locus])
		else:
			inputs.append((locus[locus.find("chr")+3:locus.find(':')],int(locus[locus.find(':')+1:])))
	return inputs

	
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-l","--loci",nargs='+',help="The the dbSNP ID or loci of interest in the format \"chr<x>:<locus>\"")
	parser.add_argument("-t","--table_fp",default="HIC059_mapq150_chr21.txt",help="TEMP: The Rao table from which to pull connections...")
	parser.add_argument("-d","--distance",type=int,default=500,help="The allowed distance from the locus of interest for a fragment to be considered...")
	parser.add_argument("-m","--mapq_cutoff",type=int,default=150,help="The minimum mapq score allowed...")
	args = parser.parse_args()
	loci = args.loci
	table_fp = args.table_fp
	distance = args.distance
	mapq_cutoff = args.mapq_cutoff
	inputs = process_inputs(loci)
	print inputs
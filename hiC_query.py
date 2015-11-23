#!/usr/bin/python
import argparse

snps = {}

def process_input(args):
	chr = args.locus[args.locus.find("chr")+3:args.locus.find(':')]
	locus = int(args.locus[args.locus.find(':')+1:])
	rao_table = open(args.table_fp,'r')
	interactions = []
	for line in rao_table:
		interaction = line.strip().split(' ')
		locus1 = int(interaction[3])
		locus2 = int(interaction[7])
		if ((interaction[2] == chr and (locus1-args.distance < locus and locus1+args.distance > locus)) or (interaction[6] == chr and (locus2-args.distance < locus and locus2+args.distance > locus))) and (int(interaction[9]) > args.mapq_cutoff and int(interaction[10]) > args.mapq_cutoff):
			if interaction[2] == chr and (locus1-args.distance < locus and locus1+args.distance > locus):
				interactions.append((interaction[6],locus2))
			else:
				interactions.append((interaction[2],locus1))
		return interactions

def assign_snps(interactions):
	for interaction in interactions:
		chr = interaction[0]
		locus = interaction[1]
		if snps[chr] == None: #If snps are not yet loaded
			snps_to_load = open("snps/bed_chr_" + chr + ".bed",'r')
			for line in snps_to_load:
				snp = line.strip().split('\t')
				snps[chr].append((int(snp[1]),snp[3]))
			snps_to_load.close()
			snps[chr].sort()
		for snp in snps[chr]:
			snp[0] = snp_locus
			if snp_locus >= locus-args.distance and snp_locus <= locus+args.distance:
				#???

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-l","--locus",help="The locus of interest in the format \"chr<x>:<locus>\"")
	parser.add_argument("-t","--table_fp",help="TEMP: The Rao table from which to pull connections.")
	parser.add_argument("-d","--distance",type=int,default=5000,help="The allowed distance from the locus of interest for a fragment to be considered.")
	parser.add_argument("-m","--mapq_cutoff",type=int,default=150,help="The minimum mapq score allowed.")
	args = parser.parse_args()
	interactions = process_input(args)
	snp_list = assign_snps(interactions)
	print interactions
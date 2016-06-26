#!/usr/bin/python
import argparse
import os
import sqlite3
import pybedtools
import requests
import multiprocessing
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import style
from matplotlib.ticker import FuncFormatter
from itertools import cycle
import re
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

def find_interactions(snps,fragment_database_fp,hic_data_dir,distance,include,exclude):
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
		interactions[snp] = {}
	for cell_line in os.listdir(hic_data_dir):
		if (include and not cell_line in include) or (exclude and cell_line in exclude):
			continue
		if os.path.isdir(hic_data_dir + '/' + cell_line):
			interactions[snp][cell_line] = Set([])
			print "\tSearching cell line " + cell_line
			for replicate in os.listdir(hic_data_dir + '/' + cell_line):
				if replicate.endswith(".db"):
					db = hic_data_dir + '/' + cell_line + '/' + replicate
					print "\t\tSearching replicate " + replicate
					table_index_db = sqlite3.connect(db)
					table_index_db.text_factory = str
					table_index = table_index_db.cursor()
					for snp in snps.keys():
                                                interactions[snp][cell_line] = Set([])
						print "\t\t\tFinding interactions for " + str(snp)
						fragment_index.execute("SELECT fragment FROM fragments WHERE chr=? AND start<=? AND end>=?",["chr" + snps[snp][0],snps[snp][1],snps[snp][1]])
						snp_fragment_result = fragment_index.fetchone()
						if snp_fragment_result == None:
							print "Warning: error retrieving SNP fragment for SNP " + snp
							continue
						snp_fragment = snp_fragment_result[0]
						snp_chr = snps[snp][0]
						for interaction in table_index.execute("SELECT chr1, fragment1 FROM chr%s_interactions WHERE chr2=? and fragment2=?" % snp_chr,[snp_chr,snp_fragment]):
							interactions[snp][cell_line].add(interaction)
						for interaction in table_index.execute("SELECT chr2, fragment2 FROM chr%s_interactions WHERE chr1=? and fragment1=?" % snp_chr,[snp_chr,snp_fragment]):
							interactions[snp][cell_line].add(interaction)
	return interactions

def find_genes(interactions,fragment_database_fp,gene_bed_fp):
	print "Identifying interactions with genes..."
	fragment_index_db = sqlite3.connect(fragment_database_fp)
	fragment_index_db.text_factory = str
	fragment_index = fragment_index_db.cursor()
	hs_gene_bed = pybedtools.BedTool(gene_bed_fp)
	genes = {}
	for snp in interactions.keys():
		#Generate BED file of all fragments interacting with SNP-containing fragment
		genes[snp] = Set([])
		temp_snp_bed = open("temp_snp_bed.bed",'w')
		for cell_line in interactions[snp].keys():
			for interaction in interactions[snp][cell_line]:
				fragment_index.execute("SELECT start, end FROM fragments WHERE chr=? and fragment=?",["chr" + interaction[0],interaction[1]])
				fragment_pos = fragment_index.fetchone()
				if fragment_pos == None:
					print "Warning: error retrieving fragment %s on chromosome %s" % (interaction[1],interaction[0])
					continue
				temp_snp_bed.write("%s\t%s\t%s\n" % ("chr" + interaction[0],fragment_pos[0],fragment_pos[1]))
		temp_snp_bed.close()
		int_bed = pybedtools.BedTool("temp_snp_bed.bed")
		#Get intersection of this BED file with BED file detailing gene locations
		gene_bed = hs_gene_bed.intersect(int_bed,u=True)
		#Return a list of genes with which SNP is interacting
		for feat in gene_bed:
			genes[snp].add(str(feat.name))
	os.remove("temp_snp_bed.bed")
	return genes

def find_eqtls(snps,genes,eqtl_data_dir,gene_database_fp,local_databases_only,num_processes):
	print "Identifying eQTLs of interest..."
	eqtls = {} #A mapping of SNPs to genes with which they have an eQTL relationship, which in turn maps to a list of tissues in which this eQTL occurs
	for snp in genes.keys():
		eqtls[snp] = {} #Stores eQTLs relevant to this SNP, and, for each eQTL, a list of the tissues in which it is found
	gene_index_db = sqlite3.connect(gene_database_fp)
	gene_index_db.text_factory = str
	gene_index = gene_index_db.cursor()
	print "\tQuerying local databases."
	for db in os.listdir(eqtl_data_dir): #Iterate through databases of eQTLs by tissue type
		tissue = db[:db.rfind('.')]
		print "\t\tQuerying " + tissue
		eqtl_index_db = sqlite3.connect(eqtl_data_dir + '/' + db)
		eqtl_index_db.text_factory = str
		eqtl_index = eqtl_index_db.cursor()
		for snp in genes.keys():
			for gene in genes[snp]:
				for eqtl in eqtl_index.execute("SELECT rsID, gene_symbol FROM eqtls WHERE rsID=? AND gene_symbol=?",(snp,gene)): #Pull down all eQTLs related to a given SNP to test for relevance:
					gene_chr = None
					gene_start = None
					gene_end = None
					if not eqtls[snp].has_key(gene):
						try:
							max_length = 0
							for gene_stat in gene_index.execute("SELECT chr, start, end, p_thresh FROM genes WHERE symbol=?", [gene]):
								if gene_stat[2] - gene_stat[1] > max_length: #Consider "canonical" to be the longest record where multiple records are present
									if gene_stat[0].startswith("chr"):
										gene_chr = gene_stat[0][gene_stat[0].find("chr")+3:]
									else:
										gene_chr = gene_stat[0]
									gene_start = gene_stat[1]
									gene_end = gene_stat[2]
									p_thresh = gene_stat[3]
									max_length = gene_stat[2] - gene_stat[1]
						except TypeError:
							print "Warning: No entry in gene database for " + geneSymbol
							continue
						eqtls[snp][gene] = {}
						eqtls[snp][gene]["gene_chr"] = gene_chr
						eqtls[snp][gene]["gene_start"] = gene_start
						eqtls[snp][gene]["gene_end"] = gene_end
						eqtls[snp][gene]["p_thresh"] = p_thresh
						eqtls[snp][gene]["tissues"] = {}
						cis = gene_chr == snps[snp][0] and (snps[snp][1] > gene_start - 1000000 and snps[snp][1] < gene_start + 1000000) #eQTL is cis if the SNP is within 1Mbp of the gene
						if cis:
							eqtls[snp][gene]["cis?"] = True
						else:
							eqtls[snp][gene]["cis?"] = False
					try:
						eqtl_index.execute("SELECT pvalue FROM eqtls WHERE rsID=? AND gene_symbol=?",(snp,gene))
						p = eqtl_index.fetchone()
						eqtls[snp][gene]["tissues"][tissue] = {"pvalue": p}
					except sqlite3.OperationalError:
						eqtls[snp][gene]["tissues"][tissue] = {"pvalue": "NA"} #TODO
	if not local_databases_only:
		print "\tQuerying GTEx online database."
		get_GTEx_response(snps,genes,gene_database_fp,eqtls,num_processes)
	return eqtls

def get_GTEx_response(snps,genes,gene_database_fp,eqtls,num_processes):
	tissues = Set(["Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Artery_Coronary","Artery_Tibial","Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus","Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Breast_Mammary_Tissue","Cells_EBV-transformed_lymphocytes","Cells_Transformed_fibroblasts","Colon_Sigmoid","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Liver","Lung","Minor_Salivary_Gland","Muscle_Skeletal","Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood"])
	gene_index_db = sqlite3.connect(gene_database_fp)
	gene_index_db.text_factory = str
	gene_index = gene_index_db.cursor()
	manager = multiprocessing.Manager()
	reqLists = [[]]
	for snp in genes.keys():
		for gene in genes[snp]:
			if not eqtls[snp].has_key(gene): #We aren't interested in eQTLs already discovered from local databases
				for tissue in tissues:
					if len(reqLists[-1]) < 1000:
						reqLists[-1].append({"snpId":snp,"gencodeId":gene,"tissueName":tissue})
					else:
						reqLists.append([{"snpId":snp,"gencodeId":gene,"tissueName":tissue}])
	print "\t\tRequests to send: " + str(len(reqLists))
	gtexResponses = manager.list()
	procPool = multiprocessing.Pool(processes=num_processes)
	for reqList in reqLists:
		procPool.apply_async(send_GTEx_query, (reqList,gtexResponses))
	procPool.close()
	procPool.join()
	print "\t\tNumber of GTEx responses received: " + str(len(gtexResponses))
	results = []
	for response in gtexResponses:
		results += response.json()["result"]
	for result in results:
		geneSymbol = result["geneSymbol"]
		snpId = result["snpId"]
		try:
			max_length = 0
			for gene_stat in gene_index.execute("SELECT chr, start, end, p_thresh FROM genes WHERE symbol=?", [geneSymbol]):
				if gene_stat[2] - gene_stat[1] > max_length: #Consider "canonical" to be the longest record where multiple records are present
					if gene_stat[0].startswith("chr"):
						gene_chr = gene_stat[0][gene_stat[0].find("chr")+3:]
					else:
						gene_chr = gene_stat[0]
					gene_start = gene_stat[1]
					gene_end = gene_stat[2]
					max_length = gene_stat[2] - gene_stat[1]
		except TypeError:
			print "\t\tWarning: No entry in gene database for " + geneSymbol
			continue
		
		cis = gene_chr == snps[snp][0] and (snps[snp][1] > gene_start - 1000000 and snps[snp][1] < gene_start + 1000000) #eQTL is cis if the SNP is within 1Mbp of the gene
		if result["pvalue"] != "NA":
			if cis and float(result["pvalue"]) >= gene_stat[3]: #If a cis eQTL has a p-value higher than its p-value cutoff, ignore the result
				continue
		if not eqtls.has_key(snpId):
			eqtls[snpId] = {}
		if not eqtls[snpId].has_key(geneSymbol):
			eqtls[snpId][geneSymbol] = {}
			eqtls[snpId][geneSymbol]["ens_id"] = result["gencodeId"]
			eqtls[snpId][geneSymbol]["gene_chr"] = gene_chr
			eqtls[snpId][geneSymbol]["gene_start"] = gene_start
			eqtls[snpId][geneSymbol]["gene_end"] = gene_end
			eqtls[snpId][geneSymbol]["p_thresh"] = gene_stat[3]
			eqtls[snpId][geneSymbol]["tissues"] = {}
			if cis:
				eqtls[snpId][geneSymbol]["cis?"] = True
			else:
				eqtls[snpId][geneSymbol]["cis?"] = False
		eqtls[snpId][geneSymbol]["tissues"][result["tissueId"]] = {"pvalue": result["pvalue"]}

def send_GTEx_query(reqList,gtexResponses):
	print "\t\tSending request to GTEx API..."
	try:
		gtexResponses.append(requests.post("http://gtexportal.org/api/v6/dyneqtl?v=clversion", json=reqList))
	except requests.exceptions.ConnectionError:
		print "\t\tWarning: a connection error occurred with the GTEx service. Retrying..."
		gtexResponses.append(requests.post("http://gtexportal.org/api/v6/dyneqtl?v=clversion", json=reqList)) #Allow to crash if fails a second time
	print "\t\tResponse received."

def produce_output(snps,interactions,eqtls,snp_database_fp,output_dir,suppress_graphs):
	print "Producing output..."
	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)
	summary_table = open(output_dir + "/summary_table.txt",'w')
	summary_table.write("SNP\tchr\tpos\tnum_contacts\tcell_lines\tnum_genes\n")
	for snp in snps.keys():
		interaction_sum = 0
		for cell_line in interactions[snp].keys():
			interaction_sum += len(interactions[snp][cell_line])
		summary_table.write(snp + '\t' + snps[snp][0] + '\t' + str(snps[snp][1]) + '\t' + str(interaction_sum) + '\t')
		for cell_line in interactions[snp].keys():
			summary_table.write(cell_line + ',')
		summary_table.write('\t' + str(len(eqtls[snp])) + '\n')
		#Write two separate summaries, cis SNP summary and trans SNP summary?
		snp_summary = open(output_dir + '/' + snp + ".txt",'w')
		snp_summary.write("dbSNP_ID\tgene_name\ttissue\tcis_p_threshold\tp-value\tchr\tstart_pos\tend_pos\tsignificant_cis_interaction\tdistance_from_snp\n")
		for gene in eqtls[snp].keys():
			distance_from_snp = 0
			if(not snps[snp][0] == eqtls[snp][gene]["gene_chr"]):
				distance_from_snp = "NA" #Not applicable to trans interactions
			elif(snps[snp][1] < eqtls[snp][gene]["gene_start"]):
				distance_from_snp = eqtls[snp][gene]["gene_start"] - snps[snp][1]
			elif(snps[snp][1] > eqtls[snp][gene]["gene_end"]):
				distance_from_snp = snps[snp][1] - eqtls[snp][gene]["gene_end"]
			eqtl_tissue = []
			for tissue in eqtls[snp][gene]["tissues"].keys():
				eqtl_tissue.append((eqtls[snp][gene]["tissues"][tissue]["pvalue"],snp,gene,tissue,str(eqtls[snp][gene]["p_thresh"]),str(eqtls[snp][gene]["gene_chr"]),str(eqtls[snp][gene]["gene_start"]),str(eqtls[snp][gene]["gene_end"]),str(eqtls[snp][gene]["cis?"]),str(distance_from_snp)))
			eqtl_tissue.sort() #Sort by p-value
			for entry in eqtl_tissue:
				snp_summary.write(entry[1] + '\t' + entry[2] + '\t' + entry[3] + '\t' + str(entry[0]) + '\t' + entry[4] + '\t' + entry[5] + '\t' + entry[6] + '\t' + entry[7] + '\t' + entry[8] + '\t' + entry[9] + '\n')
		snp_summary.close()
	summary_table.close()
	if not suppress_graphs:
		produce_graphs(snps,genes,eqtls,snp_database_fp,output_dir)

def produce_graphs(snps,genes,eqtls,snp_database_fp,output_dir):
	print "\tProducing graphs"
	style.use("ggplot")
	if not os.path.isdir(output_dir + "/plots"):
		os.mkdir(output_dir + "/plots")
	int_colours = "rgb"
	eqtl_colours = "myc"
	int_colours = cycle(int_colours)
	eqtl_colours = cycle(eqtl_colours)
	snps_by_chr = {}
	for snp, info in snps.items():
		chrm = info[0]
		if not snps_by_chr.has_key(chrm):
			snps_by_chr[chrm] = []
		snps_by_chr[chrm].append((info[1],snp))
	
	int_colour_list = []
	eqtl_colour_list = []
	num_snpgenes = []
	num_eqtls = []
	rsIDs = []
	
	chrs = snps_by_chr.keys()
	chrs.sort(key=natural_keys) #So that the chromosomes are in a logical order on the graph
	chr_locs = []
	chr_ticks = []
	chrm_pos = 0
	last_count = 0
	count = 0
	
	print "\t\tChromosomes present:"
	for chrm in chrs:
		snp_list = snps_by_chr[chrm]
		snp_list.sort() #Sort by locus
		int_colour = int_colours.next()
		eqtl_colour = eqtl_colours.next()
		chr_locs.append(chrm_pos)
		chr_ticks.append(chrm)
		for snp in snp_list:
			num_snpgenes.append(len(genes[snp[1]]))
			num_eqtls.append(len(eqtls[snp[1]]) * -1)
			rsIDs.append(snp[1])
			int_colour_list.append(int_colour)
			eqtl_colour_list.append(eqtl_colour)
			count += 1
		chrm_pos = chrm_pos + count
		print "\t\t" + chrm + " (" + str(count) + " SNPs)"
		count = 0
	
	print "\t\tChromosome locations:"
	for i,_ in enumerate(chr_locs):
		print "\t\t" + str(chr_locs[i]) + ": " + chr_ticks[i]
	
	plt.clf()
	plt.bar(range(len(num_snpgenes)),num_snpgenes,color=int_colour_list)
	plt.bar(range(len(num_eqtls)),num_eqtls,color=eqtl_colour_list)
	axes = plt.gca()
	axes.yaxis.set_major_formatter(FuncFormatter(abs_value_ticks)) #So that eQTL values won't appear negative on plot
	plt.vlines(chr_locs,axes.get_ylim()[0],axes.get_ylim()[1],colors="k")
	axes.set_xticks(range(len(rsIDs)))
	axes.set_xticklabels(rsIDs,rotation='vertical')
	ax2 = axes.twiny()
	ax2.set_xlim(axes.get_xlim())
	ax2.set_xticks(chr_locs)
	ax2.set_xticklabels(chr_ticks)
	plt.tight_layout()
	plt.savefig(output_dir + "/plots/snpgene_and_eqtl_overview.png",dpi=300,format="png")
	plt.clf()

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def abs_value_ticks(x, pos):
	return abs(x)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--inputs",nargs='+',required=True,help="The the dbSNP IDs or loci of SNPs of interest in the format \"chr<x>:<locus>\"")
	parser.add_argument("-n","--include_cell_lines",nargs='+',help="Space-separated list of cell lines to include (others will be ignored). NOTE: Mutually exclusive with EXCLUDE_CELL_LINES.")
	parser.add_argument("-x","--exclude_cell_lines",nargs='+',help="Space-separated list of cell lines to exclude (others will be included). NOTE: Mutually exclusive with INCLUDE_CELL_LINES.")
	parser.add_argument("-d","--distance",type=int,default=500,help="The allowed distance from the locus of interest for a fragment to be considered (default: 500).")
	parser.add_argument("-s","--snp_database_fp",default="snpIndex.db",help="The database of SNPs to search for details on input dbSNP IDs (Optional - will try to build from files in SNP_DIR if non-existent.")
	parser.add_argument("-t","--snp_dir",default="snps",help="The directory containing (BED) files of SNPs to use as a reference (the dbSNP Build 144 (GRCh37.p13) BED files by default).")
	parser.add_argument("-c","--hic_data_dir",default="hic_data",help="The directory containing the directories representing cell lines of HiC experiment tables.")
	parser.add_argument("-f","--fragment_bed_fp",default="Homo_sapiens.ensembl.release-74.MboI.fragments.bed",help="Bed file detailing the start and end points of HiC experiment fragments (hg19 fragmented by MboI by default.)")
	parser.add_argument("-b","--fragment_database_fp",default="fragmentIndex.db",help="The database of fragments to search when assigning SNPs to fragments (GRCh37 digested by MboI by default.)")
	parser.add_argument("-g","--gene_bed_fp",default="hg19_genes.bed",help="Bed file detailing locations of genes in genome (the UCSC known gene list for hg19 by default).")
	parser.add_argument("-a","--gene_database_fp",default="geneIndex.db",help="Database constructed from gene BED file.")
	parser.add_argument("-e","--eqtl_data_dir",default="eQTLs",help="The directory containing databases of eQTL data from various tissues (the GTEx Analysis V6 by default).")
	parser.add_argument("-o","--output_dir",default="hiCquery_output",help="The directory in which to output results (\"hiCquery_output\" by default).")
	parser.add_argument("-l","--local_databases_only",action="store_true",default=False,help="Consider only local databases. Will only include cis-eQTLs if using downloadable GTEx dataset")
	parser.add_argument("-p","--num_processes",type=int,default=1,help="Desired number of processes for multithreading (default: 1).")
	parser.add_argument("-r","--suppress_graphs",action="store_true",default=False,help="Suppress graph output (default: False).")
	args = parser.parse_args()
	snps = process_inputs(args.inputs,args.snp_database_fp,args.snp_dir)
	interactions = find_interactions(snps,args.fragment_database_fp,args.hic_data_dir,args.distance,args.include_cell_lines,args.exclude_cell_lines)
	genes = find_genes(interactions,args.fragment_database_fp,args.gene_bed_fp)
	eqtls = find_eqtls(snps,genes,args.eqtl_data_dir,args.gene_database_fp,args.local_databases_only,args.num_processes)
	produce_output(snps,interactions,eqtls,args.snp_database_fp,args.output_dir,args.suppress_graphs)

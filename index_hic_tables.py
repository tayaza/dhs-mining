#!/usr/bin/python
import argparse
import sqlite3
import os

def build_table_index(table_fp,mapq_cutoff):
	chr_set = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"]
	if os.path.isdir(table_fp[:table_fp.rfind('.')] + "_db"):
		for f in os.listdir(table_fp[:table_fp.rfind('.')] + "_db"):
			os.remove(table_fp[:table_fp.rfind('.')] + "_db/" + f)
	else:
		os.mkdir(table_fp[:table_fp.rfind('.')] + "_db")
	
	for chr in chr_set:
		table_index_db = sqlite3.connect(table_fp[:table_fp.rfind('.')] + "_db/" +  table_fp[table_fp.rfind('/'):table_fp.rfind('.')] + '_' + chr + ".db")
		table_index = table_index_db.cursor()
		table_index.execute("CREATE TABLE interactions (chr1 text, locus1 integer, chr2 text, locus2 integer)")
		table_index.execute("CREATE INDEX loci ON interactions (chr1,locus1,chr2,locus2)") ##Broken, for some reason
	##Do line count for progress meter
	print "Determining table size..."
	rao_table = open(table_fp,'r')
	lines = 0
	for i in rao_table:
		lines += 1
	rao_table.close()
	lines = lines//100*100 #Get an approximation
	
	open_dbs = []
	table_index_db = sqlite3.connect(table_fp[:table_fp.rfind('.')] + "_db/" + "temp.db")
	table_index = table_index_db.cursor()
	rao_table = open(table_fp,'r')
	print "Indexing HiC interaction table..."
	for i,line in enumerate(rao_table):
		if i % (lines/100) == 0:
			print "\tProcessed %d%%..." % ((float(i)/float(lines))*100)
		interaction = line.strip().split(' ')
		chr1 = interaction[2]
		chr2 = interaction[6]
		if int(interaction[9]) >= mapq_cutoff and int(interaction[10]) >= mapq_cutoff:
			##Hacky caching of databases - should try to find a better way to do this
			if not chr1 in open_dbs:
				if len(open_dbs) >= 10:
					table_index.execute("DETACH db" + open_dbs.pop())
				table_index.execute("ATTACH '" + table_fp[:table_fp.rfind('.')] + "_db/" +  table_fp[table_fp.rfind('/'):table_fp.rfind('.')] + '_' + chr1 + ".db' AS db" + chr1)
			else:
				open_dbs.remove(chr1)
			open_dbs.insert(0,chr1)
			table_index.execute("INSERT INTO db" + chr1 + ".interactions VALUES (?,?,?,?)",[chr1,int(interaction[3]),chr2,int(interaction[7])])
			if not chr1 == chr2:
				if not chr2 in open_dbs:
					if len(open_dbs) >= 10:
						table_index.execute("DETACH db" + open_dbs.pop())
					table_index.execute("ATTACH '" + table_fp[:table_fp.rfind('.')] + "_db/" +  table_fp[table_fp.rfind('/'):table_fp.rfind('.')] + '_' + chr2 + ".db' AS db" + chr2)
				else:
					open_dbs.remove(chr2)
				open_dbs.insert(0,chr2)
				table_index.execute("INSERT INTO db" + chr2 + ".interactions VALUES (?,?,?,?)",[chr1,int(interaction[3]),chr2,int(interaction[7])])
	rao_table.close()
	table_index_db.commit()
	os.remove(table_fp[:table_fp.rfind('.')] + '_db/' + "temp.db")
	print "Done indexing HiC interaction tables."

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-t","--table")
	parser.add_argument("-m","--mapq_cutoff",type=int,default=150)
	args = parser.parse_args()
	build_table_index(args.table,args.mapq_cutoff)
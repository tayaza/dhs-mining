#!/usr/bin/python
import argparse
import sqlite3

def build_table_index(table_fp,mapq_cutoff):
	print "Indexing HiC interaction tables..."
	table_index_db = sqlite3.connect(table_fp[:table_fp.rfind('.')] + ".db")
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
	print "Done indexing HiC interaction tables."

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-t","--table")
	parser.add_argument("-m","--mapq_cutoff",type=int,default=150)
	args = parser.parse_args()
	build_table_index(args.table,args.mapq_cutoff)
The hiC_query script expects the following file structure:

[base hiC_query directory]/
	eQTLs/
		[databases created from tissue-specific eQTL tables, provided by GTEx in the file GTEx_Analysis_V6_eQTLs.tar.gz and created using index_tables.py with the -e option]
	hic_data/
		[one directory for each cell line]/
			[one directory for each replicate]/
				<replicate_name>_<chr>.db (Databases made from HiC interaction table, e.g. tables by Rao et al. or tables in similar format, using index_tables.py with the -i option. Expects 25 files, one for each chromosome including X, Y and MT.)
	snps/ (Optional if snpIndex.db is present in base hiC_query directory)
		[BED files for dbSNP SNPs for each chromosome]
	snpIndex.db (Will be created automatically from the contents of snps/ if it does not exist)
	hiC_query.py
	index_tables.py (Not necessary for querying, but required if wishing to build databases from source tables, or use custom tables)
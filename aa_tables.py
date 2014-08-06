#! /usr/bin/env python

import pandas as pd

aa_list = ['A','V','L','I','F','W','M','P','D','E','G','S','T','C','Y','N','Q','K','R','H','x']


# Works given dataframe like test_cals.  Now need to read in fasta format and make
# the sequence names the index of the rows.

def make_binary(df,aminos=aa_list):
	binary_table = pd.DataFrame()
	for i in aminos:
		func = lambda x: 1 if x == i else 0
		bin_df = df.applymap(func)
		bin_df.columns = [i+x for x in bin_df.columns]
		binary_table = pd.concat([binary_table,bin_df],axis=1)
		## assertion statement to assure num cols is correct
	return binary_table
	
## Functions for replacing alignment DataFrames with Biochemical values

def get_biochem_table(infile):
	return pd.read_csv(infile,index_col=0) # Index column is Amino acid letters
	
def make_biochem(infile,biochem_infile,feature):
	biochem_table = get_biochem_table(biochem_infile)
	df = pd.read_csv(infile)
	ef = infile.split('_')[0]
	ef_col = pd.DataFrame({"EF":[]+[ef]*len(df)}) #create column of EF label.
	func = lambda x: biochem_table.ix[x,feature]
	biochem_df = df.applymap(func) # could have EFs in infiles and
	# skip first column with this: ix[:,1:]
	return pd.concat([ef_col,biochem_df],axis=1)
	
def biochem_rows(file_list,biochem_infile,feature):
	full_table = pd.DataFrame()
	for f in file_list:
		df = make_biochem(f,biochem_infile,feature)
		full_table = pd.concat([full_table,df],ignore_index=True)
	return full_table.reindex(columns=['EF','1','2','3','4','5','6','7','8','9','10','11','12'])
	
def write_table(df,feature):
	df.to_csv("EF_"+feature+".csv",index=False)
	
# Wishlist: function to get all biochem features. Function to make a table of all
# biochem features and all EF hands in one table. Function to remove outliers,
# though maybe this is better suited to R
		


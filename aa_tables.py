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
	
def read_from_fasta(fas_infile):
	'''Read fasta file and return generator of sequence names and sequences.'''
	name = None
	seq = []
	with open(fas_infile) as f:
		for line in f:
			if line.startswith(">"):
				if name: # is next sequence, so yield values for last seq
					yield name, ''.join(seq)
				seq = []
				name = line[1:].strip() # skip ">"
			else:
				seq.append(line.strip())
		yield name, ''.join(seq)
	
def df_from_fasta(fasta_infile):
	'''
	Given fasta infile, return pandas DataFrame where rows are sequences.
	Also includes a column called "Id" which define the protein. These are
	designated by putting the id before the sequence name, separated by an
	underscore:
	EF1_human, EF1_mouse, EF2_human...etc.
	'''
	length = None
	seq_count = 0
	for seq_name,seq in read_from_fasta(fasta_infile):
		if length == None: # is first sequence
			length = len(seq)
			df_D = {i+1:[] for i in range(length)} 	# initialize AA columns
			df_D['Id'] = []						# initialize name columns
		assert len(seq) == length, "Sequences must all be the same length"
		seq_id = seq_name.split("_")[0]
		df_D["Id"] += [seq_id]
		for i,aa in enumerate(seq):
			df_D[i+1] += [aa]
	return pd.DataFrame(df_D,columns=['Id']+[i+1 for i in range(length)])
		
	
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
		


#! /usr/bin/env python

import pandas as pd
import sys

aa_list = ['A','V','L','I','F','W','M','P','D','E','G','S','T','C','Y','N','Q','K','R','H','x']

aa_properties = {'DeltaG': {'A': 1.8100000000000001,
            'C': 1.28,
            'D': -8.7200000000000006,                         
            'E': -6.8099999999999996,
            'F': 2.98,
            'G': 0.93999999999999995,
            'H': -4.6600000000000001,
            'I': 4.9199999999999999,
            'K': -5.5499999999999998,
            'L': 4.9199999999999999,
            'M': 2.3500000000000001,
            'N': -6.6399999999999997,
            'P': 3.5800000000000001,
            'Q': -5.54,
            'R': -14.92,
            'S': -3.3999999999999999,
            'T': -2.5699999999999998,
            'V': 4.04,
            'W': 2.3300000000000001,
            'Y': -0.14000000000000001},
 'Flexibility': {'A': 0.35999999999999999,
                 'C': 0.34999999999999998,
                 'D': 0.51000000000000001,
                 'E': 0.5,
                 'F': 0.31,
                 'G': 0.54000000000000004,
                 'H': 0.32000000000000001,
                 'I': 0.46000000000000002,
                 'K': 0.46999999999999997,
                 'L': 0.37,
                 'M': 0.29999999999999999,
                 'N': 0.46000000000000002,
                 'P': 0.51000000000000001,
                 'Q': 0.48999999999999999,
                 'R': 0.53000000000000003,
                 'S': 0.51000000000000001,
                 'T': 0.44,
                 'V': 0.39000000000000001,
                 'W': 0.31,
                 'Y': 0.41999999999999998},
 'Isoelectric': {'A': 6.1100000000000003,
                 'C': 5.0199999999999996,
                 'D': 2.98,
                 'E': 3.0800000000000001,
                 'F': 5.9100000000000001,
                 'G': 6.0599999999999996,
                 'H': 7.6399999999999997,
                 'I': 6.04,
                 'K': 9.4700000000000006,
                 'L': 6.04,
                 'M': 5.7400000000000002,
                 'N': 10.76,
                 'P': 6.2999999999999998,
                 'Q': 5.6500000000000004,
                 'R': 10.76,
                 'S': 5.6799999999999997,
                 'T': 5.5999999999999996,
                 'V': 6.0199999999999996,
                 'W': 5.8799999999999999,
                 'Y': 5.6299999999999999},
 'Volume': {'A': 90.099999999999994,
            'C': 103.5,
            'D': 117.09999999999999,
            'E': 140.80000000000001,
            'F': 193.5,
            'G': 63.799999999999997,
            'H': 159.30000000000001,
            'I': 164.90000000000001,
            'K': 170.0,
            'L': 164.59999999999999,
            'M': 167.69999999999999,
            'N': 127.5,
            'P': 123.09999999999999,
            'Q': 149.40000000000001,
            'R': 192.80000000000001,
            'S': 94.200000000000003,
            'T': 120.0,
            'V': 139.09999999999999,
            'W': 231.69999999999999,
            'Y': 197.09999999999999}}


## Functions for replacing alignment DataFrames with Biochemical values

class AATable:
	
	def __init__(self,alignment):
		self.aln_frame = self.df_from_fasta(alignment)
		self.biochem_table = pd.DataFrame(aa_properties)
		self.biochem_frame = None

	def read_from_fasta(self,fas_infile):
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
		
	def df_from_fasta(self,fasta_infile):
		'''
		Given fasta infile, return pandas DataFrame where rows are sequences.
		Also includes a column called "Id" which define the protein. These are
		designated by putting the id before the sequence name, separated by an
		underscore:
		EF1_human, EF1_mouse, EF2_human...etc.
		'''
		length = None
		seq_count = 0
		for seq_name,seq in self.read_from_fasta(fasta_infile):
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
		
	def get_biochem_table(self,infile):
		'''Read in user-defined biochem table (csv) and return pandas DataFrame'''
		return pd.read_csv(infile,index_col=0) # Index column is Amino acid letters
	
	def make_biochem(self,feature,user_table=None):
		'''
		Convert alignment DataFrame into amino-acid biophysical properties.
		The property used must be specified in feature and be in a user 
		specified table (as headers), or be one of the following:
		DeltaG - Flexibility - Volume - Isoelectric
		'''
		if user_table:
			self.biochem_table = self.get_biochem_table(biochem_infile)
		df = self.aln_frame # alignment as DataFrame
		assert df, "Alignment was not read in"
		func = lambda x: self.biochem_table.ix[x,feature]
		try:
			biochem_df = df.ix[:,1:].applymap(func) # skip "Id" column
		except KeyError, value:
			raise Exception("'%s' not found in properties table" % value[0])
		self.biochem_frame = pd.concat([self.aln_frame.ix[:,:'Id'],biochem_df],axis=1)
		return self.biochem_frame
		
	def write_tables(self,user_table=None,gene_name=''):
		'''Write all tables to csv files'''
		if user_table:
			self.biochem_table = self.get_biochem_table(biochem_infile)
		for feature in self.biochem_table:
			sys.stderr.write("Writing %s table\n" % feature)
			df = self.make_biochem(feature)
			df.to_csv('.'.join([gene_name,feature,"csv"]),index=False)

	
# Wishlist: function to get all biochem features. Function to make a table of all
# biochem features and all EF hands in one table. Function to remove outliers,
# though maybe this is better suited to R
		


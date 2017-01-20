import sys, os
import csv
import numpy as np

## Class - Read data from file to array, column by column ##
 # 
 # Using, eg:
 # cols = ['index, 'source', 'temperature']                  # Columns
 # fmt  = ['i', 's', 'f']                                    # Format of each column (eg: ['s', 'i', 'f'])
 # x    = restore('filename.txt', skip_lines=4, cols, fmt)
 # dat  = x.read()
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
class restore:
	# Formats of data columns #
	fmts = ['s', 'i', 'f']

	## Initiate function ##
	 #
	 # params str file Filename
	 # params int skip Lines to skip (comment lines)
	 # params list cols Names of column data
	 # params list format Format of each column (eg: ['s', 'i', 'f', 'd'])
	 # return void
	 #
	 # version 07/2016 
	 # author Nguyen Van Hiep ##
	def __init__(self, file='', skip=0, cols=[], formats=[]):
		self.ret  = {}
		self.file = file
		self.skip = skip
		self.cols = cols
		self.fm   = formats

		self.check()

		for col in self.cols:
			self.ret[col] = []

	## Check the formats of columns, filename ##
	 #
	 # params str file Filename
	 # return boolean
	 #
	 # version 07/2016 
	 # author Nguyen Van Hiep ##
	def check(self):
		fm = self.fm
		if ((self.file == '') or (type(self.file)!= str)):
			sys.exit('No file selected!')
			return False

		if (len(fm) != len(self.cols)):
			sys.exit('Columns and Formats do not have the same length!')
			return False

		if (any((True for x in fm if x not in self.fmts))):
			sys.exit('Wrong format!')
			return False

		return True

	## Read data to colums ##
	 #
	 # params Bool asarray Convert output into array or not
	 # return dict ret data for each column
	 #
	 # version 07/2016 
	 # author Nguyen Van Hiep ##
	def read(self, asarray = False):
		file = open (self.file,'r')

		for i in range(0, self.skip):
			file.readline() # comment line

		for line in file:
		    line    = line.strip()
		    columns = line.split()

		    for k in range(0, len(self.cols)) :
		    	if(self.fm[k] == 'i') :
		    		self.ret[self.cols[k]].append(int(columns[k]))
		    	elif(self.fm[k] == 'f'):
		    		self.ret[self.cols[k]].append(float(columns[k]))
		    	else:
		    		self.ret[self.cols[k]].append(columns[k])

		file.close()
		
		if(asarray):
			for col in self.cols:
				self.ret[col] = np.asarray(self.ret[col])

		return self.ret

	## Read csv data to colums ##
	 #
	 # params Bool asarray Convert output into array or not 
	 # return dict ret data for each column
	 #
	 # version 07/2016 
	 # author Nguyen Van Hiep ##
	def readcsv(self, asarray = False):
		with open(self.file, 'rb') as f:
			## Skip comment lines
			for i in range(self.skip):
				next(f, None)

			## Read the file
			reader = csv.reader(f)
			for row in reader:
			    for k in range(0, len(self.cols)) :
			    	if(self.fm[k] == 'i') :
			    		self.ret[self.cols[k]].append(int(row[k]))
			    	elif(self.fm[k] == 'f'):
			    		self.ret[self.cols[k]].append(float(row[k]))
			    	else:
			    		self.ret[self.cols[k]].append(row[k])
		if(asarray):
			for col in self.cols:
				self.ret[col] = np.asarray(self.ret[col])
		return self.ret

### ====== EXAMPLE ===== ###
# cols = ['idx', 'src', 'nhi_fk', 'err_fk', 'nh_pl','err_pl', 'nhi_hl', 'er_hl', 'err_hl', 'wnm', 'cnm', 'fk_hl', 'pl_hl', 'oh', 'nhi_thin']
# fmt  = ['i', 's', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f']
# x    = restore('result/nhi_and_uncertainties_full.txt', 4, cols, fmt)
# dat  = x.read()
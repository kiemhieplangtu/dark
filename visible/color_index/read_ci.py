import os, sys
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np

## Common class ##
from restore             import restore # Read txt file, csv file

## Read Color index infor #
 #
 # params string fname Filename
 # return dict info of color index
 # 
 # version 11/2016
 # Author Van Hiep ##
def read_info(fname = 'data/pleiadesdata.csv'):
	cols = ['idx', 'vmag', 'bmag', 'B-V', 'vmag1', 'absv', 'modulus']
	fmt  = ['i',    'f',   'f',    'f',    'f',    's',    's']
	data = restore(fname, 5, cols, fmt)
	dat  = data.readcsv(asarray=True)
	return dat

### MAIN ###
inf = read_info('data/pleiadesdata.csv')
ci  = inf['bmag'] - inf['vmag']
y   = inf['vmag']

plt.plot(ci, y, 'rd')
plt.ylim(14.0, 0.0)
plt.xlim(-0.2, 1.4)
plt.title('CM Diagram for Pleiades')
plt.xlabel('Color Index')
plt.grid()
plt.show()
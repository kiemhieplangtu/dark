import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import pylab             as pl
import operator

# Read info of 79 sources #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_info(fname = '79src_info.txt'):
	ret = {}

	ret['src'] = []
	ret['yn']  = []
	ret['l']  = []
	ret['b']  = []
	ret['ra_icrs']  = []
	ret['de_icrs']  = []
	ret['ra_j']  = []
	ret['de_j']  = []

	file = open (fname,'r')
	file.readline() # comment line
	file.readline() # comment line

	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    ret['src'].append(columns[1])
	    ret['yn'].append(int(columns[2]))
	    ret['l'].append(float(columns[3]))
	    ret['b'].append(float(columns[4]))
	    ret['ra_icrs'].append(float(columns[5]))
	    ret['de_icrs'].append(float(columns[6]))
	    ret['ra_j'].append(str(columns[7]))
	    ret['de_j'].append(str(columns[8]))

	file.close()

	return ret


# Read info of 26 sources #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_info_no_co(fname = '26src_no_co.dat'):
	ret = {}

	ret['src']      = []
	ret['l']        = []
	ret['b']        = []
	ret['ra_icrs']  = []
	ret['de_icrs']  = []
	ret['ra_j']     = []
	ret['de_j']     = []

	ret['nhi_heiles'] = []
	ret['nhi_warm']   = []
	ret['nhi_cold']   = []

	file = open (fname,'r')
	file.readline() # comment line
	file.readline() # comment line

	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    ret['src'].append(columns[1])
	    ret['l'].append(float(columns[2]))
	    ret['b'].append(float(columns[3]))
	    ret['ra_icrs'].append(float(columns[4]))
	    ret['de_icrs'].append(float(columns[5]))
	    ret['ra_j'].append(str(columns[6]))
	    ret['de_j'].append(str(columns[7]))
	    ret['nhi_heiles'].append(float(columns[8]))
	    ret['nhi_warm'].append(float(columns[9]))
	    ret['nhi_cold'].append(float(columns[10]))

	file.close()

	return ret

# Read N(HI)_fukui N(H)_Planck N(HI)_Heiles #
#
# params string fname Filename
#
# return dict info of N(H) and N(HI)
# 
# Author Van Hiep
##
def read_nhi_fukui_nh_planck(fname = 'result/nhi_and_uncertainties_full.txt'):
	ret = {}

	ret['idx']    = []
	ret['src']    = []
	ret['nhi_fk'] = []
	ret['err_fk'] = []

	ret['nh_pl']  = []
	ret['err_pl'] = []

	ret['nhi_hl'] = []
	ret['er_hl']  = []
	ret['err_hl'] = []

	ret['wnm']    = []
	ret['cnm']    = []

	ret['fk_hl']  = []
	ret['pl_hl']  = []

	ret['oh']     = []

	file = open (fname,'r')
	file.readline() # comment line
	file.readline() # comment line
	file.readline() # comment line
	file.readline() # comment line

	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    ret['idx'].append(int(columns[0]))
	    ret['src'].append(columns[1])
	    ret['nhi_fk'].append(float(columns[2]))
	    ret['err_fk'].append(float(columns[3]))

	    ret['nh_pl'].append(float(columns[4]))
	    ret['err_pl'].append(float(columns[5]))

	    ret['nhi_hl'].append(float(columns[6]))
	    ret['er_hl'].append(float(columns[7]))
	    ret['err_hl'].append(float(columns[8]))

	    ret['wnm'].append(float(columns[9]))
	    ret['cnm'].append(float(columns[10]))

	    ret['fk_hl'].append(float(columns[11]))
	    ret['pl_hl'].append(float(columns[12]))

	    ret['oh'].append(int(columns[13]))

	file.close()

	return ret

## Only to rearrange data #
 #
 # params void
 # return void
 # 
 # Author Van Hiep ##
def rearrange_data():
	inf = read_info_no_co('26src_no_co.dat')
	dat = read_nhi_fukui_nh_planck('result/nhi_and_uncertainties_full.txt')

	for k in range(0, 26):
		print('{0}\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}\t\t{7}\t\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}'.
			format(k, inf['src'][k], inf['l'][k], inf['b'][k], inf['ra_icrs'][k], inf['de_icrs'][k], inf['ra_j'][k],inf['de_j'][k], inf['nhi_heiles'][k], inf['nhi_warm'][k], inf['nhi_cold'][k], dat['er_hl'][k], dat['err_hl'][k], dat['oh'][k]))
#================= MAIN ========================#

# Define constants #
map_file = 'data/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
rearrange_data()
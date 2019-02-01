#
# steponeplus_v2.2-to-NAMI.py
#
#     Pedro Sousa Lacerda <pslacerda@gmail.com>
#     LaBiMM, UFBA
#
# This script converts XLS data from a StepOnePlus (v2.2.2) to CSV format
#
# You need pandas and xlrd installed in order to run this script.
#    pip install pandas xlrd
#
# To run type 
#    python steponeplus_v2.2-to-NAMI.py my_raw_data.xls
# output
#    is my_raw_data_NAMI.csv
#

import sys
import csv
from os.path import splitext
import pandas as pd

# IO file paths
xls_path = sys.argv[1]
output_path = splitext(xls_path)[0] + '_NAMI.csv'

# Read sheet
sheet_name = "Multicomponent Data"
data = pd.read_excel(xls_path, sheet_name, header=7)
data.columns = [col.strip() for col in data.columns]
data = data.set_index('Reading')
data.columns = [int(w) for w in range(1, len(data.columns)+1)]

with open(output_path, 'wb') as output_file:
	writer = csv.writer(output_file)
	writer.writerow(("Temp", "Well", "Intensities"))
	for well, intensities in data.iteritems():
		for temp, intensity in intensities.iteritems():
			writer.writerow((temp, well, intensity))

#This script converts raw data from a Rotorgene 6000 to single column format
#It was calibrated with data from a Corbett Research (now Qiagen) Rotor-Gene 
#6000 series 2Plex HRM software version is 1.7 (Build 87)
#To run type 
#python rotorgene6000_v1.7-to-NAMI.py my_raw_data.csv
#output is my_raw_data_NAMI.csv

#Errors are conserved by the script but not used by NAMI. 
#Also errors are offset by 0.5 degrees and has one less entry than intensity data.
#Basically they are terrible.

import csv
import sys

input_file = open(sys.argv[1], 'rt')

try:
	reader = csv.reader(input_file)
	rownum = 1
	temps = []
	you_go_glen_coco = True 
	for row in reader:
		if rownum == 10:
			wells = (row[1:])
			intensi = {key: [] for key in wells}
			errors = {key: [] for key in wells}
		if rownum > 11:
			if len(row) == 0: #reads from row 12 to the first empty row
				you_go_glen_coco = False
			elif you_go_glen_coco: 
				temps.append(row[0])
				for key in intensi:
					intensi[key].append(row[int(key)])
		if rownum > (19 + len(temps)) and rownum < (19 + 2 * len(temps)):
#reads further down and len(temps) increases as longs as the intensities are added to the intensi dictionary
			for key in errors:
					errors[key].append(row[int(key)])
		else:
			pass
		rownum += 1
finally:
	input_file.close()

for key in errors: #Adds one entry of 0 at the end since there are fewer errors than values
	errors[key].append(0)

if sys.version_info >= (3,0,0): #This is something about Windows issues
	output_file = open(sys.argv[1][:-4]+"_NAMI.csv", 'w', newline='')
else:
	output_file = open(sys.argv[1][:-4]+"_NAMI.csv", 'wb')

try:
	writer = csv.writer(output_file)
	writer.writerow( ("Temp", "Well", "Errors", "Intensities") )
	for wellnum in wells:
		for dummy_temp in range(len(temps)):
			writer.writerow((temps[dummy_temp], wellnum, errors[wellnum][dummy_temp], intensi[wellnum][dummy_temp]))


finally: 
	output_file.close()

import glob
import yaml

# Variables
DATE = "2019-05-01"
ECSCOLLECTION = "001"
DAY = 1
DAYOFYEAR = 1
ENDTIME = "2019-05-01-T00:00:00.00000z"
ORBITCOUNT = 15



mydir = "/Users/veergadodia/Documents/Summer-2020/Harvard-CFA-Internship/OMPS-NPP-O3-Month/files"

file_list = glob.glob(mydir + "/*") # Include slash or it will search in the wrong directory!!
for i in range(len(file_list)):
    file_list[i] = file_list[i][88:]


import sys
import ruamel.yaml

yaml = ruamel.yaml.YAML()
with open('test.yaml') as fp:
    data = yaml.load(fp)

LENGTH = 0
for elem in data:
    if elem == 'Input Files':
        LENGTH = len(data[elem])
        for fi in file_list:
            data[elem].append(fi)
        data[elem] = data[elem][LENGTH:]
    if elem == 'Output File':
        data[elem][0] = 'SAO_NPP_O3T_L3_{}.nc'.format(DATE)
    if elem == 'Metadata':
        data[elem]["Date"][0] = DATE
        data[elem]["StartDate"][0] = DATE
        data[elem]["Year"][0] = int(DATE[0:4])
        data[elem]["ECSCollection"][0] = ECSCOLLECTION
        data[elem]["Day"][0] = DAY
        data[elem]["DayOfYear"][0] = DAYOFYEAR
        data[elem]["EndTime"][0] = ENDTIME
        data[elem]["OrbitCount"][0] = ORBITCOUNT

with open("test.yaml", "w") as f:
    yaml.dump(data, f)
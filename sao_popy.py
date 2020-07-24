# This code has been tested with python 2.7
# Gonzalo Gonzalez Abad, March 25, 2019
# To run it using only one cpu:
# >taskset -c 1 python sao_popy.py
import yaml # conda install -c anaconda pyyaml
import logging
import sys
import numpy as np
from sao_popy_level3 import level3
from sao_popy_level2 import tropomi
from popy import popy
import os

# Get configuration filename
config_file = 'test.yaml'

# First read config.txt file
control=yaml.load(open(config_file,'r'),
                Loader=yaml.BaseLoader)

# create logger
logger = logging.getLogger('sao_popy')
# set up log format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(module)s - %(lineno)s - %(message)s')

# set logger level and output destination
# high logging detail output to file
if control['Runtime Parameters']['DEBUG'] == '1':
    debug=True
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler('RunApp.log',mode='w')
    logger.addHandler(fh)
    fh.setFormatter(formatter)
# high logging detail output to standar output
elif control['Runtime Parameters']['DEBUG'] == '2':
    debug=True
    ch = logging.StreamHandler()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(ch)
    ch.setFormatter(formatter)
# Only output to file errors (for production)
else:
    debug=False
    logger.setLevel(logging.WARNING)
    fh = logging.FileHandler('RunApp.log',mode='w')
    logger.addHandler(fh)
    fh.setFormatter(formatter)

# Log starting message
logger.warning('Start SAO popy run')

# Log runtime parameters
logger.info('runtime parameters:\n%s' % '\n'.join('   {}={}'.format(*el) 
                                                  for el in control['Runtime Parameters'].items()))

# Log L2 file list
logger.info('generate list of input L2 files')
indir = control['Runtime Parameters']['Lv2Dir']
flist = []
for f in control['Input Files']:
    flist.append(indir+'/'+f)

# Initialize l2 object
logger.info('initialize l2 object')
try:
    l2 = tropomi(flist,
                 lon=control['Runtime Parameters']['l2_lon'],
                 lat=control['Runtime Parameters']['l2_lat'],
                 lon_bounds=control['Runtime Parameters']['l2_lon_bounds'],
                 lat_bounds=control['Runtime Parameters']['l2_lat_bounds'],
                 weights=control['Runtime Parameters']['l2_weights'])
except Exception as e:
    logger.exception('in tropomi')
    sys.exit()

out_flist=[]

# Add fields to regrid
logger.info('fields to grid:')
for key in control['Regrid Fields']:
    try:
      l2.add_grid_field(key)
      out_flist.append(control['Regrid Fields'][key])
    except Exception as e:
      logger.exception('failed adding field {0}'.format(fld))
      sys.exit()

# Add Filters
logger.info('adding filters:')
for key, conditions in control['Filter Fields'].items():
    for pair in conditions:
        try:
            l2.add_filter_field(key,pair[0],float(pair[1]))
        except Exception as e:
            logger.exception('failed adding {0} {1} {2:6.2f} filter'.
            format(key,pair[0],pair[1]))

# Init the new popy
logger.info('initialize popy')
try:
    popy_obj = popy(instrum=control['Runtime Parameters']['Instrument']
                   ,grid_size=np.float32(control['Runtime Parameters']['res'])
                   ,west=np.float32(control['Runtime Parameters']['minLon'])
                   ,east=np.float32(control['Runtime Parameters']['maxLon'])
                   ,south=np.float32(control['Runtime Parameters']['minLat'])
                   ,north=np.float32(control['Runtime Parameters']['maxLat']))
except Exception as e:
    logger.exception('initializing popy')
    sys.exit()

# Do the gridding
logger.info('performing gridding...')
try:
    popy_obj.regrid(l2)
except Exception as e:
    logger.exception('performing gridding')
    sys.exit()

# Create level 3 file
logger.info('create L3 file')
try:
    f=level3()
    f.create(control['Output File'][0],
             popy_obj.xgrid,popy_obj.ygrid)
except Exception as e:
    logger.exception('in create L3 file')
    if debug:
        sys.exit()

# Create level 3 file variables
logger.info('create L3 variables')
try:
    for key, values in control['Regrid Fields'].items():
        f.create_variable(values['var_def'],values['var_met'])
except Exception as e:
    logger.exception('in create L3 variables')
    if debug:
        sys.exit()

# Save latitude and longitude grids to level 3 file
logger.info('save lat/lon grid to L3 file')
try:
    f.set_grid(f.ncid.variables['longitude'],popy_obj.xgrid.astype(np.float32))
    f.set_grid(f.ncid.variables['latitude'],popy_obj.ygrid.astype(np.float32))
except Exception as e:
    logger.exception('in save lat/lon grid to L3 file')
    if debug:
        sys.exit()

# Save data to output variables
# Start by the variables that need to be present:
# data_quality_flag
logger.info('save data quality flag')
try:
    f.set_essential_variable(f.ncqa.variables['data_quality_flag'],
                             (popy_obj.quality_flag.astype(np.int8)),
                             popy_obj.quality_flag)
except Exception as e:
    logger.exception('in save data quality flag')
    sys.exit()

# sample_weight
logger.info('save sample weight')
try:
    f.set_essential_variable(f.ncsup.variables['sample_weight'],
                             popy_obj.total_sample_weight.astype(np.float32),
                             popy_obj.quality_flag)
except Exception as e:
    logger.exception('in save sample_weight')
    sys.exit()

# num_samples
logger.info('save number of samples')
try:
    f.set_essential_variable(f.ncqa.variables['num_samples'],
                             popy_obj.num_samples.astype(np.float32),
                             popy_obj.quality_flag)
except Exception as e:
    logger.exception('in save number of samples')
    sys.exit()

# Save the rest of the fields by looping over cotrol['Regrid Fields']
for key,values in control['Regrid Fields'].items():
    logger.info('save {0}'.format(values['var_def'][1]))
    try:
        f.set_variable(values['var_def'],
                       popy_obj.C[key],
                       popy_obj.quality_flag)
    except Exception as e:
        logger.exception('saving {0}'.format(values['var_def'][1]))
        sys.exit()

# Generate compulsory metadata: history, production_time, GranuleID
logger.info('generate default metadata dictionary')
try:
    f.set_metadata(control,popy_obj)
    logger.info('metadata:\n%s' % '\n'.join('   {}={}'.format(*el) for el in f.metadata.items()))
except Exception as e:
    logger.exception(' in generate metadata dictionary')
    if debug:
        sys.exit()

# Save metadata to L3 file
logger.info('save metadata to L3 file')
try:
    f.ncid.setncatts(f.metadata)
except Exception as e:
    logger.info(' in save metadata to L3 file')
    if debug:
        sys.exit()

# Output control file to level 3 file if need by.
if control['Runtime Parameters']['output_control?'] == 'True':
    grp=f.ncid.createGroup('control_file')
    grp.setncattr('control_file',yaml.safe_dump(control,indent=5,
                   default_flow_style=False,line_break=True,explicit_start=True))

# Close level 3 file
logger.info('close L3 file')
try:
    f.close()
except Exception as e:
    logger.exception(' in close L3 file')
    if debug:
        sys.exit()

# Rename output file to match granuleID where processing time
# has been properly set

# src=control['Output File'][0]
# dst=f.metadata['GranuleID']
# logger.info('rename file from {0} to {1}'.format(
#         src,dst))
# try:
#     os.rename(src,dst)
# except Exception as e:
#     logger.exception(' renaming L3 file')
#     sys.exit()

# Log end message
logger.warning('End SAO popy run')

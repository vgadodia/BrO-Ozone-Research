import logging

##########################################
def string_to_type(string):
    ''' Convert string to numpy type
        ARGS:
             string (str): String indicating type.
                           Accepted options are:
                           char, int8, int16, 
                           int32, float32, float64
        RETURNS:
             np.type
    '''
    import numpy as np
    import sys
    log=logging.getLogger('sao_popy.string_to_type')
    switcher = {
        'char': type('char'),
        'int8': np.int8,
        'int16': np.int16,
        'int32': np.int32,
        'float32': np.float32,
        'float64': np.float64,
    }
    # Get the type from swithcher dictionary
    output = switcher.get(string,'error')
    if output == 'error':
        log.info('These are types currently implemented:')
        log.info('char, int8, int16, float32, float64')
        log.error('can not convert {0}'.format(string))
        sys.exit()
    else:
        return output

class level3(object):

    '''
    The level3 class is used to create and read
    SAO level 3 files

    [long descript goes here]

    Note:
   
    Attributes:
    '''

########################################
    def __init__(self):
    
        ''' Initialize the level 3 SAO popoy file
            ARGS:
                NONE
        '''
        self.logger=logging.getLogger('sao_popy.level3')
        self.logger.info('creating an instance of level3')

################################################################
    def create(self,filename,lon,lat):
        ''' Create level 3 file filling up global attributes values
            ARGS:
                filename (string): Filename of the output file
                lon( 1D numpy array float32): Number of longitudes grid cells
                lat( 1D numpy array float32): Number of latitudes grid cells
        '''

        from netCDF4 import Dataset
        import numpy as np
        import sys

        # Logger
        self.logger.info('create {} level 3 file'.format(filename))
        # Save copy of initialization argument and get dimension sizes
        self.filename=filename
        self.nlon=lon.size
        self.nlat=lat.size
    
        # ######### #
        # Open file #
        # ######### #
        self.ncid = Dataset(self.filename,'w',format='NETCDF4')

        # ################# #
        # Create dimensions #
        # ################# #
        self.ncid.createDimension('longitude',self.nlon)
        self.ncid.createDimension('latitude',self.nlat)

        # ########################## #
        # Create dimension variables #
        # ########################## #
        self.lons = self.ncid.createVariable('longitude',np.float32,dimensions=('longitude'),
                                             fill_value=-1.0e+30,zlib=True,complevel=6)
        self.lons.comment='longitude at grid box center'
        self.lons.long_name='longitude'
        self.lons.units='degrees_east'
        self.lons.valid_min=np.float32(-180.)
        self.lons.valid_max=np.float32(180.)
        self.lons._Storage='contiguous'

        self.lats = self.ncid.createVariable('latitude',np.float32,dimensions=('latitude'),
                                             fill_value=-1.0e+30,zlib=True,complevel=6)
        self.lats.comment='latitude at grid box center'
        self.lats.long_name='latitude'
        self.lats.units='degrees_north'
        self.lats.valid_min=np.float32(-90.)
        self.lats.valid_max=np.float32(90.)
        self.lats._Storage='contiguous'
 
        # ####################################### #
        # Create quality assurance data variables #
        # ####################################### #
        self.ncqa=self.ncid.createGroup('qa_statistics')
        self.num_samples=self.ncqa.createVariable('num_samples',np.float32,
                                                  dimensions=('latitude','longitude'),
                                                  fill_value=-1.0,zlib=True,complevel=6)
        self.num_samples.comment=('number of samples in the calculation considering the '
                                  'summed spatial sensitivity of all satellite pixels in '
                                  'each level 3 grid box')
        self.num_samples.long_name='number of samples'
        self.num_samples.units='1'
        self.num_samples.valid_min=np.float32(0)
        self.num_samples.valid_max=np.float32(1.e+3)
        self.num_samples.coordinates='longitude latitude'
        self.num_samples._Storage='contiguous'

        self.data_quality_flag=self.ncqa.createVariable('data_quality_flag',np.int8,
                                                        dimensions=('latitude','longitude'),
                                                        fill_value=2,zlib=True,complevel=6)
        self.data_quality_flag.comment=('main data quality flag. 0 (good, number of samples > 0.1)'
                                        '1 (good, number of samples < 0.1) 2 (bad / not computed)')
        self.data_quality_flag.long_name='main data quality flag'
        self.data_quality_flag.falg_values=np.int8([0,1,2])
        self.data_quality_flag.flag_meanings=('good_number_of_samples_greater_than_0.1 '
                                              'good_number_of_samples_less_than_0.1 bad_or_not_computed')
        self.data_quality_flag.valid_min=np.int8(0)
        self.data_quality_flag.valid_max=np.int8(2)
        self.data_quality_flag.coordinates='longitude latitude'
        self.data_quality_flag._Storage='contiguous'
 
        # ################################################### #
        # Create sample_weight variable in support_data group #
        # ################################################### #
        self.ncsup=self.ncid.createGroup('support_data') 
        self.sample_weight=self.ncsup.createVariable('sample_weight',np.float32,
                                                     dimensions=('latitude','longitude'),
                                                     fill_value=-1.0e+30,zlib=True,complevel=6)
        self.sample_weight.comment='sample weight'
        self.sample_weight.long_name='sample weight'
        self.sample_weight.units='1'
        self.sample_weight.valid_min=np.float32(0.0)
        self.sample_weight.coordinates='longitude latitude'
        self.sample_weight._Storage='contiguous'

        # ############### #
        # Create metadata #
        # ############### #
        # self.ncmet=self.ncid.createGroup('metadata')
        

##########################################
    def create_variable(self,var_def,var_met):
        ''' Create variable in netCDF file
            ARGS:
                 var_def: Is a list of four members
                          --->var_def[0] group where variable is to reside (string)
                          --->var_def[1] name of the variable (string)
                          --->var_def[2] variable type (converted to numpy type by string_to_type)
                          --->var_def[3] variable fill value
                 var_met: Is a list of lists which contains:
                          --->var_met[i][0]: metadata name
                          --->var_met[i][1]: metadata value
                          --->var_met[i][2]: metadata value type
                          ...
        '''   
        from netCDF4 import Dataset
        import sys
        import numpy as np

        # Logging
        self.logger.info('create level 3 variable {0} in group {1}'.
                          format(var_def[1],var_def[0]))
        # Create group even if it exists:
        try:
            grp=self.ncid.createGroup(var_def[0])
        except:
            self.logger.error('failed creating group {0}'.format(var_def[0]))
        # Create variable
        try:
            var=grp.createVariable(var_def[1],string_to_type(var_def[2]),
                                   dimensions=('latitude','longitude'),
                                   fill_value=np.array(var_def[3]).astype(string_to_type(var_def[2])),
                                   zlib=True,complevel=6)
        except Exception as e: 
            self.logger.exception('failed creating variable {0}'.format(var_def[1]))
        # Add metadata to variable:
        for attr in var_met:
            try:
                if attr[2] == 'char':
                    var.setncattr(attr[0],attr[1])
                else:
                    var.setncattr(attr[0],np.array(attr[1]).astype(string_to_type(attr[2])))
            except:
                self.logger.warning('failed creating attribute {0} in variable {1}'.
                                     format(attr[0],var_def[1]))
                        

##########################################
    def set_grid(self,ncfvar,npvar):
        ''' Set values to netCDF variables
            ARGS:
                ncfvar: netCDF variable
                npvar: numpy array
        '''  
        self.logger.info('Set variable {}'.format(ncfvar.name))
        # Check consistency of data types
        self.check_dtype(ncfvar,npvar)
        # Check consistency of data shape and size
        self.check_dim(ncfvar,npvar)
        # Set data values
        ncfvar[:]=npvar
        
##########################################
    def set_essential_variable(self,ncfvar,npvar,mask):
        ''' Set values to netCDF variables
            ARGS:
                ncfvar: netCDF variable
                npvar: numpy array
                mask: omi_popy.quality_flag (numpy array)
        '''
        self.logger.info('Set variable {}'.format(ncfvar.name))
        # Check consistency of data types
        self.check_dtype(ncfvar,npvar)
        # Check consistency of data shape and size
        self.check_dim(ncfvar,npvar)
        # Get variable fill value
        fill=ncfvar._FillValue
        # Set variable fill value for data_quality_flag not equal to 0 or 1
        f1 = mask == 2
        npvar[f1] = fill
        # Set data values
        ncfvar[:]=npvar


##########################################
    def set_variable(self,var_def,npvar,mask):
        ''' Set values to netCDF variables
            ARGS:
                 var_def: Is a list of four members
                    --->var_def[0] group where variable is to reside (string)
                    --->var_def[1] name of the variable (string)
                    --->var_def[2] variable type (converted to numpy type by string_to_type)
                    --->var_def[3] variable fill value
                npvar: numpy array
                mask: popy_obj.quality_flag (numpy array)
        '''  
        self.logger.info('Set variable {0}'.format(var_def[1]))
        # Get level 3 variable id using var_def[0] (group) and var_def[1] (name)
        if(var_def[0] == '/'):
            ncvar=self.ncid[var_def[1]]
        else:
            dummy=self.ncid
            for grp in var_def[0].split('/'):
                dummy=dummy[grp]
            ncvar=dummy[var_def[1]]
        # Check consistency of data types
        self.check_dtype(ncvar,npvar)
        # Check consistency of data shape and size
        self.check_dim(ncvar,npvar)
        # Get variable fill value
        fill=ncvar._FillValue
        # Set variable fill value for data_quality_flag not equal to 0 or 1
        f1 = mask == 2
        npvar[f1] = fill
        # Set data values
        ncvar[:]=npvar

#######################################
    def check_dtype(self,ncfvar,npvar):
        ''' Check compatibility between data type of numpy array
            and netCDF variable
            ARGS:
                ncfvar: netCDF file variable
                npvar: numpy array (holding data to be saved in file )
        '''

        import sys
        self.logger.info('check {} variable type'.format(ncfvar.name))
        if (ncfvar.dtype.type == npvar.dtype.type):
            return True
        else:
            self.logger.error(("numpy array type {} is not compatible with "
                               "netCDF file '{}' variable {} type. Abort write ouput!!!")
                              .format(npvar.dtype.type,ncfvar.name,ncfvar.dtype.type))
            sys.exit()

#####################################
    def check_dim(self,ncfvar,npvar):
        ''' Check compatibility between data type of numpy array
            and netCDF variable
            ARGS:
                ncfvar (netCDF file variable)
                npvar (numpy variable (holding data to be saved in file )
        '''
        import sys
        self.logger.info('check {} variable dimensions'.format(ncfvar.name))
        if (ncfvar.size == npvar.size and ncfvar.shape == npvar.shape):
            return True
        else:
            self.logger.error(("numpy array shape/size {}/{} is not compatible with"
                               " netCDF file '{}' shape/size {}/{}. Abort write ouput!!!")
                              .format(npvar.shape,npvar.size,ncfvar.name,ncfvar.shape,ncfvar.size))
            sys.exit()

#####################################
    def set_metadata(self,control,omi_popy):
        ''' Create metadata dictionary
            ARGS:
                config_params (dictionary containing control.txt info)
                omi_popoy (omi_popy object)
        '''
        import sys
        import numpy as np
        import time
        self.logger.info('setting metadata values')
        granule_id='{0}_{1}m{2}{3}_v{4}-{5}.nc'.format(
            control['Metadata']['L3prodID'][0],
            control['Metadata']['StartDate'][0][0:4],
            control['Metadata']['StartDate'][0][5:7],
            control['Metadata']['StartDate'][0][8:10],
            control['Metadata']['ECSCollection'][0],
            time.strftime('%Ym%m%dt%H%M%S',time.gmtime()))
        production_time=time.strftime('%Ym%m%dt%H%M%SZ',time.gmtime())
        self.metadata={'GranuleID':granule_id,
                       'ProductionDateTime':production_time}
        for key,value in control['Metadata'].items():
            try:
                if value[1] == 'char':
                    self.metadata[key]=value[0]
                else:
                    self.metadata[key]=np.array(value[0].astype(string_to_type(value[1])))
            except:
                self.logger.info('could not set {0} metadata'.format(key))

####################
    def close(self):
        # ########## #
        # Close file #
        # ########## #
        self.logger.info('close level 3 file')
        self.ncid.close()


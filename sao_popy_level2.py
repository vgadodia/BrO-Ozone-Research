# -*- coding: iso-8859-1 -*-
# ------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------- #
import logging
def compare_fld(data,op,val):

    if(op == '!='):
        return data != val
    elif(op == '=='):
        return data == val
    elif(op == '<'):
        return data < val
    elif(op == '<='):
        return data <= val
    elif(op == '>'):
        return data > val
    elif(op == '>='):
        return data >= val
    else:
        raise Exception(op+' is not a valid operation')

def order_pts(x_pts,y_pts):

    import numpy as np

    # We need to check for +/- 180 overlap
    dx_max = np.max(x_pts,axis=1) - np.min(x_pts,axis=1) 
    for n in np.where( dx_max > 180.0 )[0]:
        x_sub = x_pts[n,:] 
        x_sub[x_sub < 0.0] = x_sub[x_sub < 0.0] + 360.0
        x_pts[n,:] = x_sub
            
    # Sort by x
    npt = x_pts.shape[0]
    idx = np.argsort(x_pts,axis=1)
    x_xsrt = np.array( [x_pts[n,:].squeeze()[idx[n,:].squeeze()] for n in range(npt)] )
    y_xsrt = np.array( [y_pts[n,:].squeeze()[idx[n,:].squeeze()] for n in range(npt)] )
    

    # Left and right most
    x_left = x_xsrt[:,:2] ; y_left = y_xsrt[:,:2]
    x_rght = x_xsrt[:,2:] ; y_rght = y_xsrt[:,2:]

    # Sort by y
    idy = np.argsort(y_left,axis=1)
    x_left = np.array( [x_left[n,:].squeeze()[idy[n,:].squeeze()] for n in range(npt)] )
    y_left = np.array( [y_left[n,:].squeeze()[idy[n,:].squeeze()] for n in range(npt)] )
    
    x_bl = x_left[:,0].squeeze() ; y_bl = y_left[:,0].squeeze()
    x_tl = x_left[:,1].squeeze() ; y_tl = y_left[:,1].squeeze()

    # do the same for the right
    idy = np.argsort(y_rght,axis=1)
    x_rght = np.array( [x_rght[n,:].squeeze()[idy[n,:].squeeze()] for n in range(npt)] )
    y_rght = np.array( [y_rght[n,:].squeeze()[idy[n,:].squeeze()] for n in range(npt)] )
    x_br = x_rght[:,0].squeeze() ; y_br = y_rght[:,0].squeeze()
    x_tr = x_rght[:,1].squeeze() ; y_tr = y_rght[:,1].squeeze()

    # Convert back to [-180,180]
    x_sort  = np.array([x_bl,x_tl,x_tr,x_br]).T 
    x_sort[x_sort > 180.0] = x_sort[x_sort > 180.0] - 360.0

    return x_sort,np.array([y_bl,y_tl,y_tr,y_br]).T

def list_available_flds(ncid):
    ''' The function has one argument (ncid) to pass a netCDF4 id as returned
    by netCDF4 library.
    '''
    # Echo all the variables fields
    try:
        for fld in ncid.variables.keys():
            print(ncid.path+':'+fld,ncid.variables[fld].shape)
    except:
        pass        
    for grp in ncid.groups.keys():
        list_available_flds(ncid[grp])

def find_group(fldname,ncid):
    ''' The function has two arguments:
       fldname (string): Name of the variable to find
       ncid (object): id object returned by netCDF4 Dataset or one of its groups
    '''
    group=''; found=False
    if (fldname in ncid.variables.keys()):
        group=ncid.path; found=True
        return group, found
    else:
        for grp in ncid.groups:
            group, found = find_group(fldname,ncid[grp])
            if found:
                return group, found
        return group, found

# ------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------- #

class measures_l2(object):

    ''' The level 2 class is used to allow measures l2 
    data to interact with the gridding code

    [long descript goes here]

    Note:
   
    Attributes:
    '''
    
    def __init__(self,file_list,time='time',lon='longitude',lat='latitude',          \
                 clon='longitude_bounds',clat='latitude_bounds',vcd='column_amount', \
                 dvcd='column_uncertainty',cf='cloud_fraction', pcld='cloud_pressure'):

        ''' Initialize the level 3 SAO OMHCHOd file
            ARGS:
                file_list[f] (str): List of L2 files
            
            OPTIONAL:
                time (str): Name of time field


                where the dimensions are:
                  f - # files
            

        '''
        from netCDF4 import Dataset

        # Initialize
        self.file_list = file_list
        self.nfile     = len(file_list)
        self.fld_list  = [] # Fields to Load
        self.grid_list = [] # Fields to grid
        self.filt_list = [] # Fields used as filters
        
        # Attach first file to get metadata
        self.fid = 0
        self.ncid = Dataset(self.file_list[0])
        
        # Set field names for required groups
        self.fld_list.append({'fldname':time,'group':self.find_fld_group(time)})
        self.fld_list.append({'fldname': lon,'group':self.find_fld_group( lon)})
        self.fld_list.append({'fldname': lat,'group':self.find_fld_group( lat)})
        self.fld_list.append({'fldname':clon,'group':self.find_fld_group(clon)})
        self.fld_list.append({'fldname':clat,'group':self.find_fld_group(clat)})
        self.fld_list.append({'fldname': vcd,'group':self.find_fld_group( vcd)})
        self.fld_list.append({'fldname':dvcd,'group':self.find_fld_group(dvcd)})
        self.fld_list.append({'fldname':  cf,'group':self.find_fld_group(  cf)})
        
        # Save popy definitions for special treatment - output transformations
        self.popydef = {}
        self.popydef[time] = 'UTC_matlab_datenum'
        self.popydef[lon ] = 'lonc'
        self.popydef[lat ] = 'latc'
        self.popydef[clon] = 'lonr'
        self.popydef[clat] = 'latr'
        self.popydef[vcd ] = 'column_amount'
        self.popydef[dvcd] = 'column_uncertainty'
        self.popydef[cf  ] = 'cloud_fraction'
        self.popydef[pcld] = 'cloud_pressure'

        # And the inverse
        self.l2def = {} 
        for name in self.popydef.keys():
            self.l2def[self.popydef[name]] = name

        # Valid operations
        self.valid_op = ['!=','==','<','<=','>','>=']

    def is_popyfld(self, l2fldname, popyfldname):

        this_is_popyfld = False

        # First check if its a popy defined field
        if( l2fldname in self.popydef ):
            if( self.popydef[l2fldname] == popyfldname):
                  this_is_popyfld = True

        return this_is_popyfld

    def print_fields(self):
        list_available_flds(self.ncid)

    def find_fld_group(self, fldname):

        # Set logical
        found = False; group = ''

        # Look for fldname
        group, found = find_group(fldname,self.ncid)
 
        # Raise exception if the group could not be found
        if(not found):
            raise Exception('Cound not find group in l2 file for '+fldname)
        else:
            print(fldname,group)
        return group

    def add_grid_field(self,fldname,group=None):
        
        ''' Add a field to the list for gridding

            ARGS:
                fldname(str): Name of field in L2 file

            OPTIONAL:
                group(str): The NetCDF group the field belongs to
        '''
        
        # Determine group if not present
        if(group is None):
            group = self.find_fld_group(fldname)

        # Dict for field to add
        ncfld = {'fldname':fldname,'group':group}

        # Append to overall load list
        if(not ncfld in self.fld_list):
            self.fld_list.append(ncfld)

        # Append to grid list
        if(not fldname in self.grid_list):
            self.grid_list.append(fldname)
        
    def add_filter_field(self,fldname,op,val,group=None):

        ''' Add a filtering criteria
            
            # Data that meets "fldname <op> val" is kept for gridding

            ARGS:
                fldname(str): Name of field in L2 file
                op(str): Operation to compare to. Follows python syntax - 
                         '!=' Not equal
                         '==' Equal
                         '<'  Less than
                         '<=' Less than or equal to
                         '>'  Greater than
                         '>=' Greater than or equal to
                val(*): threshold value to filter. Data type should be same as field

            OPTIONAL:
                group(str): The NetCDF group the field belongs to
            
            If group is not set the field will be searched for within the groups in the file
        '''

        # Determine group if not present
        if(group is None):
            group = self.find_fld_group(fldname)

        # Dict for field to add
        ncfld = {'fldname':fldname,'group':group}

        # Append to overall load list
        if(not ncfld in self.fld_list):
            self.fld_list.append(ncfld)
        
        # filter dict
        filt = {'fldname':fldname,'op':op,'val':val}

        # Append to filter list
        if(not filt in self.filt_list):
            self.filt_list.append(filt)
    
    def read_full_fields(self,file_index):
        
        from netCDF4 import Dataset
        
        # Attach file if necessary
        if( file_index != self.fid ):
            self.ncid.close()
            self.ncid = Dataset(self.file_list[file_index],'r')
            self.fid = file_index

        # Load fields
        data = {}
        for fld in self.fld_list:

            if(fld['group'] == '/'):
                data[fld['fldname']] = self.ncid.variables[fld['fldname']][:]
            else:
                data[fld['fldname']] = self.ncid.groups[fld['group']].variables[fld['fldname']][:]
        
        return data
        
    def get_filters(self,data):
        
        import numpy as np
        
        # Determine shape from longitude field
        shp = data[self.l2def['lonc']].shape
        
        # Apply the filters
        idv = np.ones(shp,dtype=bool)
        for f in self.filt_list:
            idv = idv & compare_fld(data[f['fldname']],f['op'],f['val'])
        
        return idv
    
    def apply_filters(self,data,idv):
        
        import numpy as np
        
        # Subset the data
        for fld in self.fld_list:
            
            # Special treatment for 
            if(fld['fldname'] in self.popydef.keys()):
                
                # The popy defined field name
                popy_name = self.popydef[fld['fldname']]
                
                # In addition to filtering may need to perform transformation
                if(popy_name == 'lonr' or popy_name == 'latr'):
                    ll = data[fld['fldname']][:,:,0].squeeze()[idv]
                    ul = data[fld['fldname']][:,:,1].squeeze()[idv]
                    lr = data[fld['fldname']][:,:,2].squeeze()[idv]
                    ur = data[fld['fldname']][:,:,3].squeeze()[idv]
                    data[fld['fldname']] = np.column_stack((ll,ul,ur,lr))
                elif(popy_name == 'UTC_matlab_datenum'):
                    
                    # Just make empty field for now
                    tmp = np.ones(idv.shape)[idv]
                    data[self.l2def['UTC_matlab_datenum']] = np.zeros(tmp.shape)

                    
                # Default behaviour
                else:
                    data[fld['fldname']] = data[fld['fldname']][idv]
            else:
                data[fld['fldname']] = data[fld['fldname']][idv]

        # Check order of the corners - vectorize this later
        popy_x = self.l2def['lonr'] ; popy_y = self.l2def['latr'] 
        data[popy_x][:,:],data[popy_y][:,:] = \
          order_pts(data[popy_x][:,:],data[popy_y][:,:])
          
        return data
        
    def load_subset(self,file_index):

        ''' Load subset of l2 data that meets filter thresholds

            ARGS:
                file_index: index of file to load
            
        '''
        
        # Load full data fields
        data = self.read_full_fields(file_index)
        
        # Get filters
        idv = self.get_filters(data)
        
        # Apply filters
        data = self.apply_filters(data,idv)
        
        return data
        
class omi_he5_l2(object):

    ''' The level 2 class is used to allow omi hdf-eos5 l2 
    data to interact with the gridding code

    [long descript goes here]

    Note:
   
    Attributes:
    '''

    def __init__(self,file_list,time='TimeUTC',lon='Longitude',lat='Latitude',          \
                 clon='PixelCornerLongitudes',clat='PixelCornerLatitudes',           \
                 vcd='ReferenceSectorCorrectedVerticalColumn', \
                 dvcd='ColumnUncertainty',cf='AMFCloudFraction', pcld='AMFCloudPressure'):

        ''' Initialize the level 3 SAO OMHCHOd file
            ARGS:
                file_list[f] (str): List of L2 files
            
            OPTIONAL:
                time (str): Name of time field


                where the dimensions are:
                  f - # files
            

        '''
        import h5py
        self.logger=logging.getLogger('omhchod.omi_he5_l2')
        self.logger.info('initializing l2 object')
        # Initialize
        self.file_list = file_list
        self.nfile     = len(file_list)
        self.fld_list  = [] # Fields to Load
        self.grid_list = [] # Fields to grid
        self.filt_list = [] # Fields used as filters
        
        # Attach first file to get metadata
        self.fid = 0
        self.ncid = h5py.File(self.file_list[0],'r')      
        
        # Save the swaths
        self.swth = list( self.ncid['HDFEOS']['SWATHS'].keys() )

        # Set field names for required groups
        self.fld_list.append({'fldname':time,'group':self.find_fld_group(time)})
        self.fld_list.append({'fldname': lon,'group':self.find_fld_group( lon)})
        self.fld_list.append({'fldname': lat,'group':self.find_fld_group( lat)})
        self.fld_list.append({'fldname':clon,'group':self.find_fld_group(clon)})
        self.fld_list.append({'fldname':clat,'group':self.find_fld_group(clat)})
        self.fld_list.append({'fldname': vcd,'group':self.find_fld_group( vcd)})
        self.fld_list.append({'fldname':dvcd,'group':self.find_fld_group(dvcd)})
        self.fld_list.append({'fldname':  cf,'group':self.find_fld_group(  cf)})
        
        # Save popy definitions for special treatment - output transformations
        self.popydef = {}
        self.popydef[time] = 'UTC_matlab_datenum'
        self.popydef[lon ] = 'lonc'
        self.popydef[lat ] = 'latc'
        self.popydef[clon] = 'lonr'
        self.popydef[clat] = 'latr'
        self.popydef[vcd ] = 'column_amount'
        self.popydef[dvcd] = 'column_uncertainty'
        self.popydef[cf  ] = 'cloud_fraction'
        self.popydef[pcld] = 'cloud_pressure'

        # And the inverse
        self.l2def = {} 
        for name in self.popydef.keys():
            self.l2def[self.popydef[name]] = name

        # Valid operations
        self.valid_op = ['!=','==','<','<=','>','>=']

    def is_popyfld(self, l2fldname, popyfldname):


        this_is_popyfld = False

        # First check if its a popy defined field
        if( l2fldname in self.popydef ):
            if( self.popydef[l2fldname] == popyfldname):
                  this_is_popyfld = True

        return this_is_popyfld
    
    def find_fld_group(self, fldname):


        # Set logical
        found = False
        
        # Loop over swaths
        for swth in self.swth:

            dfld = list(self.ncid['HDFEOS']['SWATHS'][swth]['Data Fields'].keys())
            gfld = list(self.ncid['HDFEOS']['SWATHS'][swth]['Geolocation Fields'].keys())

            if(fldname in (dfld)):
                group = [swth,'Data Fields']
                found = True
            if(fldname in (gfld)):
                group = [swth,'Geolocation Fields']
                found = True
        
        # Raise exception if the group could not be found
        if(not found):
            raise Exception('Could not find group in l2 file for '+fldname)
        return group
    
    def add_grid_field(self,fldname,group=None):
        
        ''' Add a field to the list for gridding

            ARGS:
                fldname(str): Name of field in L2 file

            OPTIONAL:
                group(str): The NetCDF group the field belongs to
        '''

        self.logger.info('in add_grid_field: %s',fldname)

        # Determine group if not present
        if(group is None):
            group = self.find_fld_group(fldname)

        # Dict for field to add
        ncfld = {'fldname':fldname,'group':group}

        # Append to overall load list
        if(not ncfld in self.fld_list):
            self.fld_list.append(ncfld)

        # Append to grid list
        if(not fldname in self.grid_list):
            self.grid_list.append(fldname)
    
    def add_filter_field(self,fldname,op,val,group=None):

        ''' Add a filtering criteria
            
            # Data that meets "fldname <op> val" is kept for gridding

            ARGS:
                fldname(str): Name of field in L2 file
                op(str): Operation to compare to. Follows python syntax - 
                         '!=' Not equal
                         '==' Equal
                         '<'  Less than
                         '<=' Less than or equal to
                         '>'  Greater than
                         '>=' Greater than or equal to
                val(*): threshold value to filter. Data type should be same as field

            OPTIONAL:
                group(str): The NetCDF group the field belongs to
            
            If group is not set the field will be searched for within the groups in the file
        '''

        self.logger.info('in add_filter_field: %s',fldname)

        # Determine group if not present
        if(group is None):
            group = self.find_fld_group(fldname)

        # Dict for field to add
        ncfld = {'fldname':fldname,'group':group}

        # Append to overall load list
        if(not ncfld in self.fld_list):
            self.fld_list.append(ncfld)
        
        # filter dict
        filt = {'fldname':fldname,'op':op,'val':val}

        # Append to filter list
        if(not filt in self.filt_list):
            self.filt_list.append(filt)
    
    
    def read_full_fields(self, file_index):
        
        import numpy as np
        import h5py


        # Attach file if necessary
        if( file_index != self.fid ):
            self.ncid.close()
            self.ncid = h5py.File(self.file_list[file_index],'r')
            self.fid = file_index
        
        # Load fields
        data = {}
        for fld in self.fld_list:
            swth = fld['group'][0] ; fldn = fld['fldname'] ; dgf = fld['group'][1]
            data[fldn] = self.ncid['HDFEOS']['SWATHS'][swth][dgf][fldn][:]
        
        return data
        
    def get_filters(self,data):
        
        import numpy as np
        

        # Determine shape from longitude field
        shp = data[self.l2def['lonc']].shape
        
        # Apply the filters
        idv = np.ones(shp,dtype=bool)
        for f in self.filt_list:
            idv = idv & compare_fld(data[f['fldname']],f['op'],f['val'])
        
        return idv
    
    def apply_filters(self,data,idv):
        
        import numpy as np
        
        # Subset the data
        for fld in self.fld_list:
            
            # Special treatment for 
            if(fld['fldname'] in self.popydef.keys()):
                
                # The popy defined field name
                popy_name = self.popydef[fld['fldname']]
                
                # In addition to filtering may need to perform transformation
                if(popy_name == 'lonr' or popy_name == 'latr'): 
                    ll = data[fld['fldname']][0:-1,0:-1][idv]
                    ul = data[fld['fldname']][1:,0:-1][idv]
                    lr = data[fld['fldname']][0:-1,1:][idv]
                    ur = data[fld['fldname']][1:,1:][idv]
                    data[fld['fldname']] = np.column_stack((ll,ul,ur,lr))
                elif(popy_name == 'UTC_matlab_datenum'):
                    
                    # Just make empty field for now
                    tmp = np.ones(idv.shape)[idv]
                    data[self.l2def['UTC_matlab_datenum']] = np.zeros(tmp.shape)

                    
                # Default behaviour
                else:
                    data[fld['fldname']] = data[fld['fldname']][idv]
            else:
                data[fld['fldname']] = data[fld['fldname']][idv]
            
        # Check order of the corners - vectorize this later
        popy_x = self.l2def['lonr'] ; popy_y = self.l2def['latr'] 

        if idv.any():
            data[popy_x][:,:],data[popy_y][:,:] = \
                order_pts(data[popy_x][:,:],data[popy_y][:,:])          
        return data
        
    def load_subset(self,file_index):

        ''' Load subset of l2 data that meets filter thresholds

            ARGS:
                file_index: index of file to load
            
        '''
        
        import numpy as np
        import h5py
          
        # Load the native datasets
        data = self.read_full_fields(file_index)
        
        # Get the filter indices for the data
        idv = self.get_filters(data)
        
        # Apply the filters to the dataset
        data = self.apply_filters(data,idv)

        return data

class tropomi(object):

    ''' The level 2 class is used to allow S5P??? L2 
    data to interact with the gridding code

    [long descript goes here]

    Note:
   
    Attributes:
    '''
    
    def __init__(self,file_list,lon='longitude',lat='latitude', \
                lon_bounds='longitude_bounds', lat_bounds='latitude_bounds',
                weights='column_uncertainty'):

        ''' Initialize the level 3 SAO popy file
            ARGS:
                file_list[f] (str): List of L2 files
            
            OPTIONAL: (if any of the optional parameters is not provided
                and the default is not found in the level 2 file execution
                will stop with an error message)
                lon (str): longitude centers variable in level 2 file 
                lat (str): latitude centers variable in level 2 file
                lon_bounds (str): longitude bounds variable in level 2 file
                lat_bounds (str): latitude bounds variable in level 2 file
                untertainty (str): uncertainty variable in level 2 file

            where the dimensions are:
                f - number of files
        '''
        from netCDF4 import Dataset
        
        self.logger=logging.getLogger('sao_popy.tropomi')
        self.logger.info('initializing l2 object')

        # Initialize
        self.file_list = file_list
        self.nfile     = len(file_list)
        self.fld_list  = [] # Fields to Load
        self.grid_list = [] # Fields to grid
        self.filt_list = [] # Fields used as filters
        
        # Attach first file to get metadata
        self.fid = 0
        self.ncid = Dataset(self.file_list[0])
        
        # Set field names for required groups
        self.fld_list.append({'fldname':       lon,'group':self.find_fld_group(       lon)})
        self.fld_list.append({'fldname':       lat,'group':self.find_fld_group(       lat)})
        self.fld_list.append({'fldname':lon_bounds,'group':self.find_fld_group(lon_bounds)})
        self.fld_list.append({'fldname':lat_bounds,'group':self.find_fld_group(lat_bounds)})
        self.fld_list.append({'fldname':   weights,'group':self.find_fld_group(   weights)})
        
        # Save popy definitions for special treatment - output transformations
        self.popydef = {}
        self.popydef[lon] = 'lon'
        self.popydef[lat] = 'lat'
        self.popydef[lon_bounds] = 'lon_bounds'
        self.popydef[lat_bounds] = 'lat_bounds'
        self.popydef[weights] = 'weights'

        # And the inverse
        self.l2def = {} 
        for name in self.popydef.keys():
            self.l2def[self.popydef[name]] = name

        # Valid operations
        self.valid_op = ['!=','==','<','<=','>','>=']

        # Add weights to list of grid fields
        self.add_grid_field(weights)

    def is_popyfld(self, l2fldname, popyfldname):

        this_is_popyfld = False

        # First check if its a popy defined field
        if( l2fldname in self.popydef ):
            if( self.popydef[l2fldname] == popyfldname):
                  this_is_popyfld = True

        return this_is_popyfld

    def print_fields(self):

        self.logger.info('printing fields in file {0}'.format(self.ncid.filename))
        list_available_flds(self.ncid)

    def find_fld_group(self, fldname):

        # Logging
        self.logger.info('looking for {0} variable in level 2 file'.format(fldname))

        # Set logical
        found = False; group = ''

        # Look for fldname
        group, found = find_group(fldname,self.ncid)
 
        # Raise exception if the group could not be found
        if(not found):
            self.logger.error('{0} not found in {1}'.format(fldname,self.ncid.filename))
            raise Exception('Cound not find group in l2 file for '+fldname)
        else:
            self.logger.info('{0} found in {1}'.format(fldname,group))
        return group      

    def add_grid_field(self,fldname,group=None):
        
        ''' Add a field to the list for gridding

            ARGS:
                fldname(str): Name of field in L2 file

            OPTIONAL:
                group(str): The NetCDF group the field belongs to
        '''     
        # Logging
        self.logger.info('adding grid field --> {0}'.format(fldname))

        # Determine group if not present
        if(group is None):
            group = self.find_fld_group(fldname)

        # Dict for field to add
        ncfld = {'fldname':fldname,'group':group}

        # Append to overall load list
        if(not ncfld in self.fld_list):
            self.fld_list.append(ncfld)

        # Append to grid list
        if(not fldname in self.grid_list):
            self.grid_list.append(fldname)
        
    def add_filter_field(self,fldname,op,val,group=None):

        ''' Add a filtering criteria    
            # Data that meets "fldname <op> val" is kept for gridding
            ARGS:
                fldname(str): Name of field in L2 file
                op(str): Operation to compare to. Follows python syntax - 
                         '!=' Not equal
                         '==' Equal
                         '<'  Less than
                         '<=' Less than or equal to
                         '>'  Greater than
                         '>=' Greater than or equal to
                val(*): threshold value to filter. Data type should be same as field
            OPTIONAL:
                group(str): The NetCDF group the field belongs to
            If group is not set the field will be searched for within the groups in the file
        '''
        # Logging
        self.logger.info('adding filter --> {0} {1} {2}'.
            format(fldname,op,val))

        # Determine group if not present
        if(group is None):
            group = self.find_fld_group(fldname)

        # Dict for field to add
        ncfld = {'fldname':fldname,'group':group}

        # Append to overall load list
        if(not ncfld in self.fld_list):
            self.fld_list.append(ncfld)
        
        # filter dict
        if (op == '!=' or op == '==' or op == '<' or
        op == '<=' or op == '>' or op == '>='):
            filt = {'fldname':fldname,'op':op,'val':val}
        else:
            self.logger.error('operation {0} not allowed for {1} field'.
                            format(fldname,op))
            sys.exit()

        # Append to filter list
        if(not filt in self.filt_list):
            self.filt_list.append(filt)
    
    def read_full_fields(self,file_index):
        
        from netCDF4 import Dataset
        
        # Attach file if necessary
        if( file_index != self.fid ):
            self.ncid.close()
            self.ncid = Dataset(self.file_list[file_index],'r')
            self.fid = file_index

        # Load fields
        data = {}
        for fld in self.fld_list:
            if(fld['group'] == '/'):
                data[fld['fldname']] = self.ncid.variables[fld['fldname']][:]
            else:
                dummy=self.ncid
                for grp in fld['group'].split('/')[1:]:
                    dummy=dummy[grp]
                data[fld['fldname']] = dummy[fld['fldname']][:].squeeze()
                del dummy
        
        return data
        
    def get_filters(self,data):
        
        import numpy as np
        
        # Determine shape from longitude field
        shp = data[self.l2def['lon']].shape
        
        # Apply the filters
        idv = np.ones(shp,dtype=bool)
        for f in self.filt_list:
            idv = idv & compare_fld(data[f['fldname']],f['op'],f['val'])
        
        return idv
    
    def apply_filters(self,data,idv):
        
        import numpy as np
        
        # Subset the data
        for fld in self.fld_list:
            self.logger.info('appliying filter for {0}'.format(fld['fldname']))
            # Special treatment for pixel bounds
            if(fld['fldname'] in self.popydef.keys()):
                # The popy defined field name
                popy_name = self.popydef[fld['fldname']]    
                # In addition to filtering may need to perform transformation
                if(popy_name == 'lon_bounds' or popy_name == 'lat_bounds'):
                    ll = data[fld['fldname']][:,:,0].squeeze()[idv]
                    ul = data[fld['fldname']][:,:,1].squeeze()[idv]
                    lr = data[fld['fldname']][:,:,2].squeeze()[idv]
                    ur = data[fld['fldname']][:,:,3].squeeze()[idv]
                    data[fld['fldname']] = np.column_stack((ll,ul,ur,lr))
                # Default behaviour
                else:
                    data[fld['fldname']] = data[fld['fldname']][idv]
            else:
                data[fld['fldname']] = data[fld['fldname']][idv]

        # Check order of the corners - vectorize this later
        self.logger.info('checking order of corners')
        popy_x = self.l2def['lon_bounds'] ; popy_y = self.l2def['lat_bounds'] 
        data[popy_x][:,:],data[popy_y][:,:] = \
          order_pts(data[popy_x][:,:],data[popy_y][:,:])
          
        return data
        
    def load_subset(self,file_index):

        ''' Load subset of l2 data that meets filter thresholds

            ARGS:
                file_index: index of file to load
            
        '''
        self.logger.info('loading data fields')
        # Load full data fields
        data = self.read_full_fields(file_index)
        self.logger.info('loading filters')
        # Get filters
        idv = self.get_filters(data)
        self.logger.info('appliying filters')
        # Apply filters
        data = self.apply_filters(data,idv)
        
        return data



# flist = ['.//S5P_OFFL_L2__O3_____20200423T173931_20200423T192102_13100_01_010108_20200425T103720.nc']

# l2 = tropomi(flist)
# print(l2)

# from netCDF4 import Dataset
# root = Dataset("S5P_OFFL_L2__O3_____20200423T173931_20200423T192102_13100_01_010108_20200425T103720.nc")
# print(root.groups["PRODUCT"].groups["SUPPORT_DATA"].groups["GEOLOCATIONS"].variables.keys())


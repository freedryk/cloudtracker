import numpy, ConfigParser, glob, os

config = ConfigParser.RawConfigParser()
config.read('config.cfg')
model_config = {}

for option in ('ug', 'vg', 'dt', 'dz', 'dy', 'dx'):
	model_config[option] = config.getfloat('modelconfig', option)
for option in ('nz', 'ny', 'nx'):
	model_config[option] = config.getint('modelconfig', option)
for option in ('case_name', 'input_directory', 'data_directory', 'sam_directory'):
	model_config[option] = config.get('modelconfig', option)

model_config['do_entrainment'] = config.getboolean('modelconfig', 'do_entrainment')

nz, ny, nx = model_config['nz'], model_config['ny'], model_config['nx']
dt, dx, dy, dz = model_config['dt'], model_config['dz'], model_config['dy'], model_config['dz']

ug, vg = model_config['ug'], model_config['vg']

case_name = model_config['case_name']
do_entrainment = model_config['do_entrainment']

input_directory = model_config['input_directory']
data_directory = model_config['data_directory']
sam_directory = model_config['sam_directory']

nt = len(glob.glob('%s/%s_[!A-Z]*.bin3D' % (input_directory, case_name)))

def get_stat():
	filename = next(glob.iglob(data_directory + '/*_stat.nc'))
	return filename

def time_picker(file_name):
	f = file_name.split('/')[-1].split('_')
	
	if(f[1] == 'CORE'):
		filelist = glob.glob('%s/core_entrain/*.nc' % data_directory)
	elif (f[1] == 'CLOUD'):
		filelist = glob.glob('%s/condensed_entrain/*.nc' % data_directory)
	else:
		filelist = glob.glob('%s/variables/*.nc' % data_directory)
	
	filelist.sort()
	index = filelist.index(file_name)
	return index

def index_to_zyx(index):
    z = index / (ny*nx)
    index = index % (ny*nx)
    y = index / nx
    x = index % nx
    return numpy.array((z, y, x))
                
def zyx_to_index(z, y, x):
    return ny*nx*z + nx*y + x

def index_to_array_3d(index):
    temp_array = numpy.zeros((nz*ny*nx,), numpy.bool) 
    temp_array[index] = 1
    return temp_array.reshape((nz, ny, nx))

def calc_com(mask):
    pts = index_to_zyx( mask )

    z = pts[0,:].astype(float).mean()
    # Correct Center of Mass for reentrant domain
    y1 = pts[1,:].astype(float)
    x1 = pts[2,:].astype(float)
    y2 = (y1 < ny/2.)*y1 + (y1>= ny/2.)*(y1 - ny)
    x2 = (x1 < nx/2.)*x1 + (x1>= nx/2.)*(x1 - nx)
    y1m = y1.mean()
    y2m = y2.mean()
    x1m = x1.mean()
    x2m = x2.mean()
    
    if numpy.var(y2 - y2m) > numpy.var(y1 - y1m):
        y = y1m
    else:
        y = (y2m + .5)%ny - .5
        
    if numpy.var(x2 - x2m) > numpy.var(x1 - x1m):
        x = x1m
    else:
        x = (x2m + .5)%nx - .5
        
    return numpy.array((z, y, x))

def calc_distance(point1, point2):
    # Calculate distances corrected for reentrant domain
    delta_x = numpy.abs(point2[2] - point1[2])
    if delta_x >= (nx/2): delta_x = nx - delta_x
    delta_y = numpy.abs(point2[1] - point1[1])
    if delta_y >= (ny/2): delta_y = ny-delta_y
    delta_z = point2[0] - point1[0]
    return numpy.sqrt(delta_x**2 + delta_y**2 + delta_z**2)

def expand_indexes(indexes):
    # Expand a given set of indexes to include the nearest
    # neighbour points in all directions.
    # indexes is an array of grid indexes
                    
    K_J_I = index_to_zyx( indexes )

    stack_list = [K_J_I, ]
    for item in ((-1, 0, 0), (1, 0, 0),
                 (0, -1, 0), (0, 1, 0), 
                 (0, 0, -1), (0, 0, 1)):
        stack_list.append( K_J_I + numpy.array(item)[:, numpy.newaxis] )
    
    expanded_index = numpy.hstack(stack_list)

    # re-entrant domain
    expanded_index[0, expanded_index[0, :] == nz] = nz-1
    expanded_index[0, expanded_index[0, :] < 0] = 0
    expanded_index[1, :] = expanded_index[1, :]%ny
    expanded_index[2, :] = expanded_index[2, :]%nx

    # convert back to indexes
    expanded_index = zyx_to_index(expanded_index[0, :],
                                  expanded_index[1, :],
                                  expanded_index[2, :])
                                  
    expanded_index = numpy.unique(expanded_index)
    
    return expanded_index



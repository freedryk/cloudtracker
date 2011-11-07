NOTE: The most up-to-date version of cloudtracker may be found at 

https://github.com/freedryk/cloudtracker

-----------------------------------

To run cloudtracker, the following Python modules must be installed:

numpy
networkx

And *either*:

netcdf4-python or pupynere

---------------------------

The model_config.cfg file contains model-specific parameters that need to 
be set for your model, and the location of a directory containing the input 
NetCDF files.

--------------------------

cloudtracker takes a set of NetCDF files as input.  There should be one
file per model time step, and the files should be named so that an ASCII 
sort will put them in order of increasing time.  The easiest way to do 
this is simply by appending an index to the files, like so:

file_0000001.nc
file_0000002.nc
file_0000003.nc
...

cloudtracker will treat every file in the input directory as a valid 
input file, so ensure only the NetCDF files for processing are present
in the directory.

Each NetCDF file should contain the following integer variables:

core: 1 for buoyant, upward moving model points with condensed liquid water,
      0 otherwise
condensed: 1 for model points with condensed liquid water 
      0 otherwise
plume: 1 for model points conforming to the tracer criteria defined in 
       Couvreax et al. (2010)
      0 otherwise

And the following floating point variables:
u: east-west velocity
v: north-south velocity
w: up-down velocity

All variables should be of size (nz, ny, nx), as defined in model_config.cfg.

-------------------------

To run cloudtracker, type the following at a command line:

./track_clouds.py model_config.cfg

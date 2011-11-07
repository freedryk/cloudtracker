#!/usr/bin/env python

import numpy
import glob
import cPickle
import sys
import os
import os.path

from generate_cloudlets import generate_cloudlets
from cluster_cloudlets import cluster_cloudlets
from make_graph import make_graph
from output_cloud_data import output_cloud_data

try:
   from netCDF4 import Dataset
except:
   try:
       from netCDF3 import Dataset
   except:
       from pupynere import netcdf_file as Dataset

#-------------------

def load_data(filename):
    input_file = Dataset(filename)
        
    core = input_file.variables['core'][:].astype(bool)
    condensed = input_file.variables['condensed'][:].astype(bool)
    plume = input_file.variables['plume'][:].astype(bool)
    u = input_file.variables['u'][:].astype(double)
    v = input_file.variables['v'][:].astype(double)
    w = input_file.variables['w'][:].astype(double)
        
    input_file.close()

    return core, condensed, plume, u, v, w

#---------------

def main(MC, save_all=True):
    sys.setrecursionlimit(100000)
    input_dir = MC['input_directory']
    nx = MC['nx']
    ny = MC['ny']
    nz = MC['nz']
    nt = MC['nt']
    
    filelist = glob.glob('%s/*' % input_dir)
    filelist.sort()
    
    if (len(filelist) != nt):
        raise "Only %d files found, nt=%d files expected" % (len(filelist), nt)

    if not os.path.exists('pkl'):
        os.mkdir('pkl')
    if not os.path.exists('output'):
        os.mkdir('output')

    for n, filename in enumerate(filelist):
        print "generate cloudlets; time step: %d" % n
        core, condensed, plume, u, v, w = load_data(filename)

        cloudlets = generate_cloudlets(core, condensed, plume, 
                                       u, v, w, MC)
    
        cPickle.dump(cloudlets, open('pkl/cloudlets_%08g.pkl' % n,'wb'))

#----cluster----

    cluster_cloudlets(MC)
        
#----graph----

    print "make graph"

    cloud_graphs, cloud_noise = make_graph(MC)
    
    print "\tFound %d clouds" % len(cloud_graphs)

    if save_all:
        cPickle.dump((cloud_graphs, cloud_noise), open('pkl/graph_data.pkl', 'wb'))
            
#----output----

    for n in range(nt):
        print "output cloud data, time step: %d" % n
        output_cloud_data(cloud_graphs, cloud_noise, n, MC)
            

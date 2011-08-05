#!/usr/bin/env python
# Runtime (690, 130, 128, 128): 1 hour 20 minutes

import cPickle
import networkx
import numpy
from utility_functions import zyx_to_index, index_to_zyx, calc_radii, \
    expand_indexes
import sys
#import scipy.io

def calc_shell(index, MC):
    # Expand the cloud points outward
    maskindex = expand_indexes(index, MC)
    
    # From the expanded mask, select the points outside the cloud
    shellindex = numpy.setdiff1d(maskindex, index, assume_unique=True)

    return shellindex

def calc_edge(index, shellindex, MC):
    # Find all the points just inside the clouds
    maskindex = expand_indexes(shellindex, MC)
    
    # From the expanded mask, select the points inside the cloud
    edgeindex = numpy.intersect1d(maskindex, index, assume_unique=True)

    return edgeindex


def calc_env(index, shellindex, edgeindex, MC):
    if len(shellindex) > 0:
        K_J_I = index_to_zyx(shellindex, MC)

        n = 6
        for i in range(n):
            # Find all the points just outside the clouds
            stacklist = [K_J_I, ]
            for item in ((0, -1, 0), (0, 1, 0),
                         (0, 0, -1), (0, 0, 1)):
                stacklist.append( K_J_I + numpy.array(item)[:, numpy.newaxis] )

            maskindex = numpy.hstack(stacklist)
            maskindex[1, :] = maskindex[1, :] % MC['ny']
            maskindex[2, :] = maskindex[2, :] % MC['nx']            
            maskindex = numpy.unique( zyx_to_index(maskindex[0, :],
                                                   maskindex[1, :],
                                                   maskindex[2, :],
                                                   MC) )

            # From the expanded mask, select the points outside the cloud
            envindex = numpy.setdiff1d(maskindex, index, assume_unique = True)

            K_J_I = index_to_zyx(envindex, MC)
            
        # Select the points within 4 grid cells of cloud
        r = calc_radii(envindex, edgeindex, MC)
        mask = r < 4.5
        envindex = envindex[mask]
    else:
        envindex = []

    return envindex


def calculate_data(cluster, MC):
    result = {}

    result['plume'] = cluster['plume']

    condensed = cluster['condensed']
    result['condensed'] = condensed
    condensed_shell = calc_shell(condensed, MC)
    result['condensed_shell'] = condensed_shell
    condensed_edge = calc_edge(condensed, condensed_shell, MC)
    result['condensed_edge'] = condensed_edge
    result['condensed_env'] = calc_env(condensed, condensed_shell, condensed_edge, MC)

    core = cluster['core']
    result['core'] = core
    core_shell = calc_shell(core, MC)
    result['core_shell'] = core_shell
    core_edge = calc_edge(core, core_shell, MC)
    result['core_edge'] = core_edge
    result['core_env'] = calc_env(core, core_shell, core_edge, MC)

    return result

def save_text_file(clouds, t, MC):
    count = 0

    for id in clouds:
        for point_type in clouds[id]:
            count = count + len(clouds[id][point_type])

    recarray = numpy.zeros(count, dtype=[('id', 'i4'),('type', 'a14'),('x','i4'),('y', 'i4'), ('z', 'i4')])

    count = 0
    for id in clouds:
        for point_type in clouds[id]:
            data = clouds[id][point_type]
            n = len(data)
            if n == 0: continue
            recarray['id'][count:n + count] = id
            recarray['type'][count:n + count] = point_type
            z, y, x = index_to_zyx(data, MC)
            recarray['x'][count:n + count] = x
            recarray['y'][count:n + count] = y
            recarray['z'][count:n + count] = z
            count = count + n
            
    recarray.tofile(open('output/clouds_at_time_%08g.txt' % t, 'w'), '\r\n')


def output_cloud_data(cloud_graphs, cloud_noise, t, MC):

    print 'Timestep:', t

    # Load the cluster zyx data for the current time
    clusters = {}
    cluster_dict = cPickle.load(open('pkl/clusters_%08g.pkl' % t, 'rb'))
    for id in cluster_dict:
        key = "%08g|%08g" % (t, id)
        clusters[key] = cluster_dict[id]

    # If t
    clouds = {}
    id = 0
    for subgraph in cloud_graphs:
        # Grab the nodes at the current time 
        # that all belong to subgraph 'id'
        nodes = [item for item in subgraph.nodes() 
                      if item[:8] == ('%08g' % t)]
                      
        if nodes:
            # Pack them into a single cloud object
            core = []
            condensed = []
            plume = []
            for node in nodes:
                core.append(clusters[node]['core'])
                condensed.append(clusters[node]['condensed'])
                plume.append(clusters[node]['plume'])
                
            cloud = {'core': numpy.hstack(core),
                     'condensed': numpy.hstack(condensed),
                     'plume': numpy.hstack(plume)}

            # Calculate core/cloud, env, shell and edge
            clouds[id] = calculate_data(cloud, MC)
        id = id + 1

    # Add all the noise to a noise cluster
    noise_clust = {'core': [], 'condensed': [], 'plume': []}
    for subgraph in cloud_noise:
        nodes = [item for item in subgraph.nodes() 
                       if item[:8] == ('%08g' % t)]
        if nodes:
            for node in nodes:
                noise_clust['core'].append(clusters[node]['core'])
                noise_clust['condensed'].append(clusters[node]['condensed'])
                noise_clust['plume'].append(clusters[node]['plume'])
                    
    if noise_clust['core']:                    
        noise_clust['core'] = numpy.hstack(noise_clust['core'])         
    if noise_clust['condensed']: 
        noise_clust['condensed'] = numpy.hstack(noise_clust['condensed'])
    if noise_clust['plume']:
        noise_clust['plume'] = numpy.hstack(noise_clust['plume'])

    # Only save the noise if it contains cloud core
    clouds[-1] = calculate_data(noise_clust, MC)
            
    print "Number of Clouds at Current Timestep: ", len(clouds.keys())

    cPickle.dump(clouds, open('pkl/cloud_data_%08g.pkl' % t,'wb'))

    save_text_file(clouds, t, MC)

#   save .mat file for matlab
#    new_dict = {}
#    for key in clouds:
#        new_dict['cloud_%08g' % key ] = clouds[key]
    
#    savedict = {'clouds': new_dict} 
#    scipy.io.savemat('mat/cloud_data_%08g.mat' % t, savedict)
    


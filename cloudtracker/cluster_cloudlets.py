#!/usr/bin/env python
# Runtime (690, 130, 128, 128): 1.5 hours

import numpy
import cPickle

from cloud_objects import Cloudlet, Cluster
from utility_functions import index_to_zyx, zyx_to_index

saveit = True

#--------------------------

def make_spatial_cloudlet_connections(cloudlets, MC):
    # Find all the cloudlets which have adjacent cores or plumes
    # Store the information in the cloudlet.adjacent dict.
    
    # This function does this by constructing a 3d array of the 
    # cloudlet numbers of each core and plume point
    # Then it pulls the edge points from these 3d arrays
    # to see if the cloudlets are bordering other cloudlets
    
    condensed_array = -1*numpy.ones((MC['nz']*MC['ny']*MC['nx'],), numpy.int)
    plume_array = -1*numpy.ones((MC['nz']*MC['ny']*MC['nx'],), numpy.int)

    # label the cloud core and plume points using the list index of the
    # cloudlet
    for cloudlet in cloudlets:
        condensed_array[cloudlet.condensed_mask()] = cloudlet.id
        plume_array[cloudlet.plume_mask()] = cloudlet.id

    for cloudlet in cloudlets:
        
        # Find all cloudlets that have adjacent clouds
        adjacent_condensed = condensed_array[cloudlet.condensed_halo()]
        adjacent_condensed = adjacent_condensed[adjacent_condensed > -1]
        if len(adjacent_condensed) > 0:
            volumes = numpy.bincount(adjacent_condensed)
            adjacent_condensed = numpy.unique(adjacent_condensed)
            for id in adjacent_condensed:
                cloudlet.adjacent['condensed'].append((volumes[id], cloudlets[id]))
            cloudlet.adjacent['condensed'].sort()
            cloudlet.adjacent['condensed'].reverse()

        # Find all cloudlets that have adjacent plumes
        adjacent_plumes = plume_array[cloudlet.plume_halo()]
        adjacent_plumes = adjacent_plumes[adjacent_plumes > -1]
        if len(adjacent_plumes) > 0:
            volumes = numpy.bincount(adjacent_plumes)
            adjacent_plumes = numpy.unique(adjacent_plumes)
            for id in adjacent_plumes:
                cloudlet.adjacent['plume'].append((volumes[id], cloudlets[id]))
            cloudlet.adjacent['plume'].sort()
            cloudlet.adjacent['plume'].reverse()

    return cloudlets

#-----------------

def advect_indexes(indexes, u, v, w, MC):
    K_J_I = index_to_zyx(indexes, MC)
    K_J_I[0, :] = K_J_I[0, :] - w
    K_J_I[1, :] = K_J_I[1, :] - v
    K_J_I[2, :] = K_J_I[2, :] - u
                               
    K_J_I[0, K_J_I[0, :] >= MC['nz']] = MC['nz']-1
    K_J_I[0, K_J_I[0, :] < 0] = 0
    K_J_I[1, :] = K_J_I[1, :] % MC['ny']
    K_J_I[2, :] = K_J_I[2, :] % MC['nx']

    advected_indexes = zyx_to_index(K_J_I[0,:], K_J_I[1,:], K_J_I[2,:], MC)
    
    return advected_indexes

def count_overlaps(key, overlaps, cloudlet):
    bin_count = numpy.bincount(overlaps)
    indexes = numpy.arange(len(bin_count))
    indexes = indexes[bin_count > 0]
    bin_count = bin_count[bin_count > 0]
    for n, index in enumerate(indexes):
        cloudlet.overlap[key].append( (bin_count[n],  index) )

def make_temporal_connections(cloudlets, old_clusters, MC):
    # For each cloudlet, find the previous time's
    # cluster that overlaps the cloudlet the most
    condensed_array = -1*numpy.ones((MC['nz']*MC['ny']*MC['nx'],), numpy.int)
    plume_array = -1*numpy.ones((MC['nz']*MC['ny']*MC['nx'],), numpy.int)

    # label the cloud core and plume points using the list index of the
    # cloud cluster
    for id, cluster in old_clusters.iteritems():
        condensed_array[cluster.condensed_mask()] = id
        plume_array[cluster.plume_mask()] = id

    for cloudlet in cloudlets:
        # Find cloud-cloud overlaps
        if cloudlet.has_condensed():
            # Correct for cloud advection
            advected_condensed_mask = advect_indexes(cloudlet.condensed_mask(),
                                                cloudlet.u['condensed'],
                                                cloudlet.v['condensed'],
                                                cloudlet.w['condensed'],
                                                MC)

            # Get indexes of previous cores from the advected array
            overlapping_condenseds = condensed_array[advected_condensed_mask]
            overlapping_condenseds = overlapping_condenseds[overlapping_condenseds > -1]
            
            if len(overlapping_condenseds) > 0:
                count_overlaps('condensed->condensed', overlapping_condenseds, cloudlet)
                            
            # Find core-plume overlaps
            overlapping_plumes = plume_array[advected_condensed_mask]
            overlapping_plumes = overlapping_plumes[overlapping_plumes > -1]

            if len(overlapping_plumes) > 0:
                count_overlaps('plume->condensed', overlapping_plumes, cloudlet)
            
        # Find plume-core overlaps
        advected_plume_mask = advect_indexes(cloudlet.plume_mask(),
                                             cloudlet.u['plume'],
                                             cloudlet.v['plume'],
                                             cloudlet.w['plume'],
                                             MC)
                                             
        overlapping_condenseds = condensed_array[advected_plume_mask]
        overlapping_condenseds = overlapping_condenseds[overlapping_condenseds > -1]
        if len(overlapping_condenseds) > 0:
            count_overlaps('condensed->plume', overlapping_condenseds, cloudlet)

        # Find plume-plume overlaps
        overlapping_plumes = plume_array[advected_plume_mask]
        overlapping_plumes = overlapping_plumes[overlapping_plumes > -1]
        if len(overlapping_plumes) > 0:
            count_overlaps('plume->plume', overlapping_plumes, cloudlet)
                
        for item in cloudlet.overlap:
            cloudlet.overlap[item].sort()
            cloudlet.overlap[item].reverse()
            
#---------------------

def create_new_clusters(cloudlets, clusters, max_id, MC):

    core_list = []
    condensed_list = []
    plume_list = []
    for cloudlet in cloudlets:
        if cloudlet.has_core():
            core_list.append(cloudlet)
        elif cloudlet.has_condensed():
            condensed_list.append(cloudlet)
        else:
            plume_list.append(cloudlet)
    
    n = 0

    # Make clusters out of the cloudlets with core points
    while core_list:
        cloudlet = core_list.pop()
        cluster = Cluster(max_id, [cloudlet], MC)
        cluster.events.append('NCOR')
        # Add cloudlets with adjactent clouds to the cluster
        # Adding cloudlets may bring more cloudlets into cloud
        # contact with the cluster, so we loop until acloud is empty
        acondenseds = cluster.adjacent_cloudlets('condensed')
        while acondenseds:
            n = n + len(acondenseds)
            for cloudlet in acondenseds:
                try:
                    core_list.remove( cloudlet )
                except:
                    cloud_list.remove( cloudlet )
            cluster.add_cloudlets( acondenseds )
            acondenseds = cluster.adjacent_cloudlets('condensed')

        clusters[max_id] = cluster
        max_id = max_id + 1

    # Make clusters out of the cloudlets without core points
    while condensed_list:
        cloudlet = condensed_list.pop()
        cluster = Cluster(max_id, [cloudlet], MC)
        cluster.events.append('NCLD')
        if (len(cluster.adjacent_cloudlets('condensed')) > 0): print "        condensed connection ERROR"

    # Make clusters out of the cloudlets without core points
    while plume_list:
        cloudlet = plume_list.pop()
        cluster = Cluster(max_id, [cloudlet], MC)
        cluster.events.append('NP')
        clusters[max_id] = cluster
        max_id = max_id + 1

    return clusters

#---------------------

def associate_cloudlets_with_previous_clusters(cloudlets, old_clusters, MC):
    clusters = {}
    new_cloudlets = []

    for cloudlet in cloudlets:
        back_conns = set()
        max_conn = -1

        if cloudlet.overlap['condensed->condensed']:
            conns = cloudlet.overlap['condensed->condensed']
            max_conn = conns[0][1]
            conns = conns[1:]
            for conn in conns:
                back_conns.add(conn[1])
        elif cloudlet.overlap['plume->condensed']:  
            conns = cloudlet.overlap['plume->condensed']
            for conn in conns:
                if not old_clusters[conn[1]].has_condensed():
                    if max_conn > -1:
                        back_conns.add(max_conn)
                    else:
                        max_conn = conn[1]
        elif cloudlet.overlap['plume->plume']:
            if not cloudlet.has_condensed():
                conns = cloudlet.overlap['plume->plume']
                for conn in conns:
                    if not old_clusters[conn[1]].has_condensed():
                        if max_conn > -1:
                            back_conns.add(max_conn)
                        else:
                            max_conn = conn[1]
 

        # If there are back connections, add the cloudlet to
        # a cluster
        if max_conn > -1:
            if max_conn in clusters:
                clusters[max_conn].add_cloudlet(cloudlet)
            else:
                clusters[max_conn] = Cluster(max_conn, [cloudlet], MC)
                clusters[max_conn].events.append('O%d' % max_conn)
                clusters[max_conn].past_connections.add(max_conn)
            for conn in back_conns:
                clusters[max_conn].merge_connections.add(conn)
                clusters[max_conn].events.append('M%d' % conn)
        else:
            new_cloudlets.append( cloudlet )
 
    return new_cloudlets, clusters

#---

def check_for_adjacent_cloudlets(new_cloudlets, clusters):
    n = 0
    # Checks the clusters list to see if any of the cloudlets which did not
    # overlap previous clusters are connected to the current clusters
    for cluster in clusters.values():
        condensed_connections = cluster.adjacent_cloudlets('condensed')
        while condensed_connections:
            n = n + 1
            connected_cloudlet = condensed_connections.pop()
            if connected_cloudlet in new_cloudlets:
                cluster.add_cloudlet( connected_cloudlet )
                new_cloudlets.remove( connected_cloudlet )
                condensed_connections = cluster.adjacent_cloudlets('condensed')

#---

def split_clusters(clusters, max_id, MC):

    for cluster in clusters.values():
        groups = cluster.connected_cloudlet_groups()
        if len(groups) > 1:
            sizes = []
            for group in groups:
                size = 0
                for cloudlet in group:
                    size = size + cloudlet.volume
                sizes.append( (size, group) )

            sizes.sort()

            # Turn the smaller groups into new clusters
            for size, group in sizes[:-1]:
                cluster.remove_cloudlets(group)
                new_cluster = Cluster(max_id, group, MC)
                new_cluster.events.append('S%d' % cluster.id)
                new_cluster.split_connections.add(cluster.id)
                clusters[max_id] = new_cluster
                max_id = max_id + 1
                
    return max_id
                
#----------

def make_clusters(cloudlets, old_clusters, MC):
    # make_clusters generates a dictionary of clusters

    max_id = max(old_clusters.keys()) + 1

    # Find the horizontal connections between cloudlets
    make_spatial_cloudlet_connections(cloudlets, MC)

    # associate cloudlets with previous timestep clusters
    # cloudlets that can't be associated are assumed to be newly created
    new_cloudlets, current_clusters = associate_cloudlets_with_previous_clusters(cloudlets, old_clusters, MC)

    # See if any of the new cloudlets are touching a cluster
    check_for_adjacent_cloudlets(new_cloudlets, current_clusters)

    # See if the cloudlets in a cluster are no longer touching
    max_id = split_clusters(current_clusters, max_id, MC)

    # Create new clusters from any leftover new cloudlets
    final_clusters = create_new_clusters(new_cloudlets, current_clusters, max_id, MC)

    return final_clusters

#---------------------

def load_cloudlets(t, MC):
    cloudlets = cPickle.load(open('pkl/cloudlets_%08g.pkl' % t,'rb'))

    result = []
    n = 0
    for cloudlet in cloudlets:
        if ((len(cloudlet['plume']) > 7) 
            or (len(cloudlet['condensed']) > 1)
            or (len(cloudlet['core']) > 0)):
            result.append( Cloudlet( n, t, cloudlet, MC ) )
            n = n + 1

    return result

def save_clusters(clusters, t):
    new_clusters = {}
    for id, clust in clusters.iteritems():
        new_dict = {}
        new_dict['past_connections'] = clust.past_connections
        new_dict['merge_connections'] = clust.merge_connections
        new_dict['split_connections'] = clust.split_connections
        new_dict['events'] = clust.events
        new_dict['core'] = clust.core_mask()
        new_dict['condensed'] = clust.condensed_mask()
        new_dict['plume'] = clust.plume_mask()
        new_clusters[id] = new_dict
    cPickle.dump(new_clusters, open('pkl/clusters_%08g.pkl' % t, 'wb'))
    cPickle.dump(clusters, open('pkl/cluster_objects_%08g.pkl' % t, 'wb'))

def cluster_cloudlets(MC):

    print "cluster cloudlets; time step: 0"
    cloudlets = load_cloudlets(0, MC)    
    make_spatial_cloudlet_connections( cloudlets, MC )
    new_clusters = create_new_clusters(cloudlets, {}, 0, MC)
    print "\t%d clusters" % len(new_clusters)
    save_clusters(new_clusters, 0)
        
    for t in range(1, MC['nt']):
        print "cluster cloudlets; time step: %d" % t
        old_clusters = new_clusters
        cloudlets = load_cloudlets(t, MC)

        # Finds the ids of all the previous timestep's cloudlets that overlap
        # the current timestep's cloudlets.
        make_temporal_connections(cloudlets, old_clusters, MC)

        # Uses the previous timestep overlap info to group
        # current cloudlets into clusters.
        new_clusters = make_clusters(cloudlets, old_clusters, MC)
        print "\t%d clusters" % len(new_clusters)

        save_clusters(new_clusters, t)

    
if __name__ == "__main__":
    main()


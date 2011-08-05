#!/usr/bin/env python

import numpy
from utility_functions import expand_indexes, \
        index_to_zyx, find_halo, calc_distance

def calc_com(mask, MC):
    ny, nx = MC['ny'], MC['nx']
    pts = index_to_zyx( mask, MC )

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

    if numpy.var(y2 - y2m) > numpy.var(y1 - y1m):
        y = y1m
    else:
	y = (y2m + .5)%ny - .5

    if numpy.var(x2 - x2m) > numpy.var(x1 - x1m):
        x = x1m
    else:
	x = (x2m + .5)%nx - .5

    return numpy.array(([z], [y], [x]))

#----------------------------

#def calc_distance(point1, point2, MC):
#    # Calculate distances corrected for reentrant domain
#    point1 = numpy.atleast_2d(point1)
#    point2 = numpy.atleast_2d(point2)
#    ny, nx = MC['ny'], MC['nx']
#    
#    delta_x = numpy.abs(point2[2, :] - point1[2, :])
#    mask = delta_x >= (nx/2)
#    delta_x[mask] = nx - delta_x[mask]
#    
#    delta_y = numpy.abs(point2[1, :] - point1[1, :])
#    mask = delta_y >= (ny/2)
#    delta_y[mask] = ny - delta_y[mask]
#    
#    delta_z = point2[0, :] - point1[0, :]
#    
#    return numpy.sqrt(delta_x**2 + delta_y**2 + delta_z**2)

#----------------------------

class Cloudlet:
    def __init__(self, id, time, cloudlet_dict, MC):
        self.id = id
        self.time = time
        self.MC = MC

        self.masks = {}
        for item in ('core', 'condensed', 'plume'):
            self.masks[item] = cloudlet_dict[item]

        self.adjacent = {'core': [],
                         'condensed': [],
                         'plume': []}

        self.overlap = {'condensed->condensed': [],
                        'condensed->plume': [],
                        'plume->condensed': [],
                        'plume->plume': []}
       
        self.u = {'condensed': cloudlet_dict['u_condensed'],
                  'plume': cloudlet_dict['u_plume']}
        self.v = {'condensed': cloudlet_dict['v_condensed'],
                  'plume': cloudlet_dict['v_plume']}
        self.w = {'condensed': cloudlet_dict['w_condensed'],
                  'plume': cloudlet_dict['w_plume']}
                  
        self.cluster = None

        if self.has_core():
            self.volume = len(self.core_mask())
        elif self.has_condensed():
            self.volume = len(self.condensed_mask())
        else:
            self.volume = len(self.plume_mask())

    def has_core(self):
        return len(self.masks['core']) > 0

    def core_mask(self):
        return self.masks['core']

    def has_condensed(self):
        return len(self.masks['condensed']) > 0

    def condensed_mask(self):
        return self.masks['condensed']

    def condensed_halo(self):
        return find_halo( self.condensed_mask(), self.MC )

    def plume_mask(self):
        return self.masks['plume']

    def plume_halo(self):
        return find_halo( self.plume_mask(), self.MC )

#----------------------------

class Cluster:
    def __init__(self, cluster_id, initial_cloudlets, MC):
        self.id = cluster_id
        self.cloudlets = set()
        self.past_connections = set()
        self.split_connections = set()
        self.merge_connections = set()
        self.add_cloudlets(initial_cloudlets)
        self.events = []
        self.MC = MC

    def add_cloudlet(self, cloudlet):
        if not cloudlet.cluster:
            cloudlet.cluster = self
            self.cloudlets.add(cloudlet)
        else:
            raise "Cloudlet already belongs to a cluster!"

    def add_cloudlets(self, cloudlets):
        for cloudlet in cloudlets:
            self.add_cloudlet( cloudlet )

    def remove_cloudlets(self, cloudlets):
        for cloudlet in cloudlets:
            if cloudlet.cluster:
                cloudlet.cluster = None
                self.cloudlets.remove( cloudlet )
            else:
                raise "Cloudlet does not belong to cluster!"

    def has_core(self):
        for cloudlet in self.cloudlets:
            if cloudlet.has_core(): return True
        return False

    def core_mask(self):
        clist = []
        for cloudlet in self.cloudlets:
            clist.append(cloudlet.masks['core'])
        return numpy.hstack(clist)

    def has_condensed(self):
        for cloudlet in self.cloudlets:
            if cloudlet.has_condensed(): return True
        return False

    def condensed_mask(self):
        clist = []
        for cloudlet in self.cloudlets:
            clist.append(cloudlet.masks['condensed'])
        return numpy.hstack(clist)


    def plume_mask(self):
        clist = []
        for cloudlet in self.cloudlets:
            clist.append(cloudlet.masks['plume'])
        return numpy.hstack(clist)

    def condensed_halo(self):
        return find_halo( self.condensed_mask(), self.MC )

    def plume_halo(self):
        return find_halo( self.plume_mask(), self.MC )

    def adjacent_cloudlets(self, key):
        result = {}
        for cloudlet in self.cloudlets:
            for volume, adjacent_cloudlet in cloudlet.adjacent[key]:
                if adjacent_cloudlet not in self.cloudlets:
                    if adjacent_cloudlet in result:
                        result[adjacent_cloudlet] += volume
                    else:
                        result[adjacent_cloudlet] = volume
                        
        final = [(result[cloudlet], cloudlet) for cloudlet in result]
        final.sort()
        final.reverse()
        result = [item[1] for item in final]
        return result

    def connected_cloudlet_groups(self):
        # only split if both components are connected to ground
        plume_cloudlets = []
        condensed_cloudlets = []
        for cloudlet in self.cloudlets:
            if cloudlet.has_condensed():
                condensed_cloudlets.append(cloudlet)
            else:
                plume_cloudlets.append(cloudlet)

        groups = []
        while condensed_cloudlets:
            cloudlet = condensed_cloudlets.pop()
            group = [cloudlet, ]
            conns = cloudlet.adjacent['condensed'][:]
            for vol, cloudlet in conns:
                if cloudlet in condensed_cloudlets:
                    group.append(cloudlet)
                    condensed_cloudlets.remove(cloudlet)
                    conns.extend(cloudlet.adjacent['condensed'][:])
            groups.append(group)
        
        if len(groups) > 1:
            detached_groups = []
            attached_groups = []
            for group in groups:
                mask = []
                for cloudlet in group:
                    mask.append(cloudlet.plume_mask())
                mask = numpy.hstack(mask)
                    
                if (mask < self.MC['nx']*self.MC['ny']).any():
                    attached_groups.append(group)
                else:
                    detached_groups.append(group)

            if attached_groups:
                volumes = []
                volume = 0
                for group in attached_groups:
                    for cloudlet in group:
                        volume += cloudlet.volume
                    volumes.append((volume, group))                    
                volumes.sort()
                
                if detached_groups:
                    attached_group_masks = []
                    for group in attached_groups:
                        mask_list = [cloudlet.condensed_mask() for cloudlet in group]
                        mask = numpy.hstack(mask_list)
                        attached_group_masks.append((calc_com(mask, self.MC), group))
    
                    for item in detached_groups:
                        mask_list = [cloudlet.condensed_mask() for cloudlet in item]
                        item_mask = numpy.hstack(mask_list)
                        item_com = calc_com(item_mask, self.MC)

                        com_list = [(calc_distance(item_com, current_group[0], self.MC), current_group[1]) 
                                    for current_group in attached_group_masks]
                        com_list.sort()
                        com_list[0][1].extend(item)
                    
                groups = attached_groups
            else:
                groups = detached_groups


        while plume_cloudlets:
            cloudlet = plume_cloudlets.pop()
            group = [cloudlet, ]
            conns = cloudlet.adjacent['plume'][:]
            for volume, cloudlet in conns:
                if cloudlet in plume_cloudlets:
                        group.append(cloudlet)
                        plume_cloudlets.remove(cloudlet)
                        conns.extend(cloudlet.adjacent['plume'][:])
            groups.append(group)

        return groups



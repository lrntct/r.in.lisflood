#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
MODULE:    r.in.lisflood

AUTHOR(S): Laurent Courty

PURPOSE:   Reads LISFLOOD-FP input files and create relevant GRASS raster maps.
           Output could be used as input of t.sim.flood

COPYRIGHT: (C) 2015 by Laurent Courty
           This program is free software under the GNU General Public
           License (v3). Read the LICENCE file for details.
"""

#%module
#% description: Import relevant maps from LISFLOOD-FP entry datas
#% keywords: raster
#% keywords: import
#%end

#%option G_OPT_F_INPUT
#% key: par_file
#% description: LISFLOOD-FP input *.par file
#% required: yes
#%end

#%option G_OPT_R_OUTPUT
#% key: dem
#% description: Name of output DEM raster map
#% required: yes
#%end

#%option G_OPT_R_OUTPUT
#% key: friction
#% description: Name of output Manning's n raster map
#% required: yes
#%end

#%option G_OPT_R_OUTPUT
#% key: start_file
#% description: Name of output starting water depth raster map
#% required: yes
#%end

#%option G_OPT_R_INPUT
#% key: bc
#% description: Name of output boundary conditions type raster map
#% required: no
#%end

#%option G_OPT_R_INPUT
#% key: bcval
#% description: Name of output boundary conditions value raster map
#% required: no
#%end

#%option G_OPT_R_INPUT
#% key: user_inflow
#% description: Name of output user inflow raster map
#% required: no
#%end

import grass.script as grass
import grass.temporal as tgis
from grass.pygrass import raster
from grass.pygrass import utils
from grass.pygrass.gis.region import Region
from grass.pygrass.messages import Messenger
from grass.pygrass.modules import Module

import numpy as np
import chardet
import codecs
import sys
import os

def main():
    # start messenger
    msgr = Messenger()

    # Store current region
    region = Region()

    # reads CLI options
    rast_n_file_name = options['friction']
    rast_dem_name = options['dem']
    rast_start_file_name = options['start_file']
    rast_bc_name = options['bc']
    rast_bcval_name = options['bcval']
    rast_user_name = options['user_inflow']

    # Load *.par file
    par = Par(msgr, options['par_file'])
    par.read()

    # Write DEM file
    par.write_dem(rast_dem_name, grass.overwrite())
    # set computational region to match DEM
    par.set_region_from_map(rast_dem_name)
    # Write friction
    par.write_n_map(rast_n_file_name, grass.overwrite())
    # Write start file
    par.write_start_h(rast_start_file_name, grass.overwrite())

    # *.bci file
    if par.bci_file and options['bcval'] and options['bc']:
        bci_full_path = os.path.join(par.directory, par.bci_file)
        bci = Bci(msgr, bci_full_path, region=par.region)
        bci.read()
        bci.write_fixed_bc(rast_bc_name, rast_bcval_name,
                        rast_user_name, grass.overwrite())

    # *.bdy file
    if par.bdy_file:
        bdy_full_path = os.path.join(par.directory, par.bdy_file)
        bdy = Bdy(msgr, bdy_full_path, bci.bdy_kwd,
            bcvar=bci.bcvar, region=par.region)
        bdy.read()



    # Make sure the original region is restored
    region.write()
    return 0


class Par(object):
    '''
    '''

    # valid *.par keywords
    kwd = {'dem_file':'DEMfile',
            'bci_file':'bcifile',
            'bdy_file':'bdyfile',
            'start_file':'startfile',
            'friction':'fpfric',
            'n_file':'manningfile',
            'sim_time':'sim_time'}


    def __init__(self, msgr, par_file):
        self.msgr = msgr
        self.par_file = par_file  # full path
        self.directory = os.path.dirname(self.par_file)
        self.dem_file = ''
        self.n_file = ''
        self.start_file = ''
        self.bci_file = ''
        self.bdy_file = ''
        self.sim_time = ''
        self.n = ''


    def read(self):
        '''
        '''
        file_encoding = chardet.detect(open(self.par_file, "r").read())
        file_open = codecs.open(
                    self.par_file, encoding=file_encoding['encoding'])
        with file_open as input_file:
            for line in input_file:
                # remove starting or trailing white spaces
                line = line.strip()
                if line.startswith(self.kwd['dem_file']):
                    self.dem_file = line.split()[1]
                if line.startswith(self.kwd['bci_file']):
                    self.bci_file = line.split()[1]
                if line.startswith(self.kwd['bdy_file']):
                    self.bdy_file = line.split()[1]
                if line.startswith(self.kwd['start_file']):
                    self.start_file = line.split()[1]
                if line.startswith(self.kwd['sim_time']):
                    self.sim_time = line.split()[1]
                if line.startswith(self.kwd['n_file']):
                    self.n_file = line.split()[1]
                if line.startswith(self.kwd['friction']):
                    self.n = line.split()[1]
                if line.startswith('latlong'):
                    self.msgr.fatal('lat-long import is not supported')
        return self


    def write_dem(self, raster_name, overwrite):
        '''
        '''
        if self.dem_file:
            r_in_gdal = Module("r.in.gdal",
                input=os.path.join(self.directory, self.dem_file),
                output=raster_name, overwrite=overwrite)
        else:
            self.msgr.fatal('No {} found'.format(self.kwd['dem_file']))
        return self


    def write_n_map(self, raster_name, overwrite):
        '''Write a GRASS map from informations of *.par file
        in case n_file is provided, any fpfric value is discarded
        raster_name = name of grass raster to be written (string)
        overwrite = boolean value
        '''
        if not self.n_file and not self.n:
            self.msgr.fatal('Either {} or {} should be provided'.format(
                            par_kwd['n_file'], par_kwd['friction']))
        if self.n and not self.n_file:
            mapcalc_expression = '{map} = {value}'.format(
                            map=raster_name, value=self.n)
            Module("r.mapcalc",
                expression=mapcalc_expression,
                overwrite=overwrite)
        if self.n_file:
            r_in_gdal = Module("r.in.gdal",
                input=os.path.join(self.directory, self.n_file),
                output=raster_name, overwrite=overwrite)
        return self


    def write_start_h(self,  raster_name, overwrite):
        '''
        '''
        if self.start_file:
            r_in_gdal = Module("r.in.gdal",
                input=os.path.join(self.directory, self.start_file),
                output=raster_name, overwrite=overwrite)
        return self


    def set_region_from_map(self, map_name):
        '''
        '''
        # get region from the DEM
        self.region = Region()
        self.region.from_rast(map_name)
        if self.region.ewres != self.region.nsres:
            self.mesgr.fatal('Non-square cells not supported by lisflood!')
        self.region.write()
        return self


class Bci(object):
    '''
    '''

    # valid key letters for boundary location
    klt = ['N', 'S', 'E', 'W', 'P', 'F']
    # valid boundary types
    bc_type = ['CLOSED', 'FREE',
                'QFIX', 'QVAR',
                'HFIX', 'HVAR']
    # Correspondance between LISFLOOD-FP boundary conditions and t.sim.flood
    # QFIX and QVAR are taken into account in another way
    bc_conv = {'CLOSED':1,
                'FREE':2,
                'HFIX':3,
                'HVAR':3}


    def __init__(self, msgr, bci_file, region=None):
        self.bci_file = bci_file  # full path to the file
        self.msgr = msgr
        self.content = []
        self.region = region
        self.bdy_kwd = []
        self.bcvar = {}


    def read(self):
        '''
        '''
        file_encoding = chardet.detect(open(self.bci_file, "r").read())
        file_open = codecs.open(self.bci_file,
                        encoding=file_encoding['encoding'])
        with file_open as input_file:
            for line_num, line in enumerate(input_file, 1):
                # transform in a list without leading and trailing spaces
                line = line.strip().split()
                if not line or '#' in line[0]:
                    continue
                if line[0] in self.klt:
                    self.content.append(line)
                if line[3] in ('HVAR', 'QVAR'):
                    self.bdy_kwd.append(line[4])
                if line[3] not in self.bc_type:
                    self.msgr.fatal(
                        'Unknown boundary type {} at line {}'.format(
                                            line[3], line_num))
        return self


    def write_fixed_bc(self, rast_type_name, rast_value_name,
                        rast_user_name, overwrite):
        '''Generate and write fixed input flow and
        boundary conditions type and value maps
        '''
        # Boundary conditions type arrays
        arr_bctype = grass.array.array(dtype=np.uint8)
        arr_bc_N = arr_bctype[0,:]
        arr_bc_S = arr_bctype[-1,:]
        arr_bc_E = arr_bctype[:,0]
        arr_bc_W = arr_bctype[:,-1]

        # Boundary conditions value arrays
        arr_bcvalue = grass.array.array(dtype=np.float32)
        arr_bcval_N = arr_bcvalue[0,:]
        arr_bcval_S = arr_bcvalue[-1,:]
        arr_bcval_E = arr_bcvalue[:,0]
        arr_bcval_W = arr_bcvalue[:,-1]

        # User inflow arrays
        arr_user = grass.array.array(dtype=np.float32)
        arr_user_N = arr_user[0,:]
        arr_user_S = arr_user[-1,:]
        arr_user_E = arr_user[:,0]
        arr_user_W = arr_user[:,-1]

        for line in self.content:
            if line[0] == 'N':
                self.set_fixed_bc(bc_pos=self.region.north, line=line,
                    reg_min=self.region.west, reg_max=self.region.east,
                    arr_type=arr_bc_N, arr_value=arr_bcval_N,
                    arr_user=arr_user_N)
            if line[0] == 'S':
                self.set_fixed_bc(bc_pos=self.region.south, line=line,
                    reg_min=self.region.west, reg_max=self.region.east,
                    arr_type=arr_bc_S, arr_value=arr_bcval_S,
                    arr_user=arr_user_S)
            if line[0] == 'E':
                self.set_fixed_bc(bc_pos=self.region.east, line=line,
                    reg_min=self.region.south, reg_max=self.region.north,
                    arr_type=arr_bc_E, arr_value=arr_bcval_E,
                    arr_user=arr_user_E)
            if line[0] == 'W':
                self.set_fixed_bc(bc_pos=self.region.west, line=line,
                    reg_min=self.region.south, reg_max=self.region.north,
                    arr_type=arr_bc_W, arr_value=arr_bcval_W,
                    arr_user=arr_user_W)
            if line[0] == 'P' and line[3] == 'QFIX':
                coord = utils.coor2pixel(
                        (float(line[1]), float(line[2])), self.region)
                arr_user[coord] = float(line[4])
            if line[0] == 'P' and line[3] == 'QVAR':
                coord = utils.coor2pixel(
                        (float(line[1]), float(line[2])), self.region)
                self.bcvar[line[4]] = {'loc':line[0],
                                    'type':line[3],
                                    'bc_len':0,
                                    'coord':coord}

        # write maps in GRASS
        arr_bctype.write(mapname=rast_type_name, overwrite=overwrite)
        if not np.count_nonzero(arr_bcvalue) == 0:
            arr_bcvalue.write(mapname=rast_value_name, overwrite=overwrite)
        if not np.count_nonzero(arr_user) == 0:
            arr_user.write(mapname=rast_user_name, overwrite=overwrite)
        return self


    def set_fixed_bc(self, bc_pos, reg_min, reg_max, line,
                arr_type, arr_value, arr_user):
        '''bc_pos: coordinate of the considered boundary in coordinate.
                    north, east, west or south
        reg_min: region min boundary
        reg_max: region min boundary
        line: line content
        arr_res_type: numpy array
        arr_res_value: numpy array
        arr_user: numpy array
        '''
        # crop coordinates to fit in region
        coord_geo1 = max(float(line[1]), reg_min)
        coord_geo2 = min(float(line[2]), reg_max)
        bc_len = coord_geo2 - coord_geo1
        if bc_len < 0:
             self.msgr.fatal(
                'Incoherent coordinates \n {}'.format(line))

        # transform into array coordinates
        coord_arr_1 = utils.coor2pixel(
                            (coord_geo1, bc_pos), self.region)[1]
        coord_arr_2 = utils.coor2pixel(
                            (coord_geo2, bc_pos), self.region)[1]
        # assign type
        if line[3] not in ('QFIX', 'QVAR'):
            arr_type[coord_arr_1:coord_arr_2] = self.bc_conv[line[3]]

        # flow in m/s for QFIX
        if line[3] == 'QFIX':
            arr_user[coord_arr_1:coord_arr_2] = (
                        float(line[4]) * bc_len /
                        self.region.ewres * self.region.nsres)
        # assign value
        if line[3] == 'HFIX':
            arr_value[coord_arr_1:coord_arr_2] = float(line[4])
        # store location of variable conditions
        if line[3] in ('HVAR', 'QVAR'):
            self.bcvar[line[4]] = {'loc':line[0],
                                    'type':line[3],
                                    'bc_len':bc_len,
                                    'coord':(coord_arr_1, coord_arr_2)}
        return self


class Bdy(object):
    """
    """

    kwd_units = ['days', 'hours', 'seconds']

    def __init__(self, msgr, bdy_file, bdy_kwd, bcvar, region=None):
        self.bdy_file = bdy_file
        self.bdy_kwd = bdy_kwd
        self.region = region
        self.bcvar = bcvar

        self.content = {}
        self.unit = {}


    def read(self):
        '''
        '''
        file_content = []
        file_encoding = chardet.detect(open(self.bdy_file, "r").read())
        file_open = codecs.open(self.bdy_file,
                        encoding=file_encoding['encoding'])
        # Read the file and transform it to a list of lists
        with file_open as input_file:
            for line_num, line in enumerate(input_file, 1):
                if line_num == 1:
                    continue
                line = line.strip().split()
                file_content.append(line)
        for line in file_content:
            if line[0] in self.bdy_kwd and line[0] not in self.content:
                current_section_name = line[0]
                current_section_line = line_num
                self.content[current_section_name] = []
                self.unit[current_section_name] = ''
            if line_num == current_section_line + 1:
                current_section_unit = Line[1]
                self.unit[current_section_name] = current_section_unit
                if current_section_unit not in self.kwd_units:
                    self.msgr.fatal(
                'Unknown unit: {} for boundary {}'.format(
                    current_section_unit, current_section_name))
            if is_number(line[0]) and is_number(line[1]):
                # likely a value line
                self.content[current_section_name].append(
                    (float(line[0]), int(line[1])))
        return self


def is_number(s):
    try:
        float(s)
    except ValueError:
        return False
    else:
        return True


if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())

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
from grass.pygrass.gis import Mapset
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

    # Use temporary GRASS region
    grass.use_temp_region()

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

    # boundary conditions
    bc = BoundaryConditions(msgr, sim_time=par.sim_time,
                        region=par.region)
    if par.bci_file:
        bci_full_path = os.path.join(par.directory, par.bci_file)
        bc.read_bci(bci_full_path)
    if par.bdy_file:
        bdy_full_path = os.path.join(par.directory, par.bdy_file)
        bc.read_bdy(bdy_full_path)
    print bc.content

    # Restore original region
    grass.del_temp_region()
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
            'rain_file':'rainfall',
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
        self.rain_file = ''
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
                if line.startswith(self.kwd['rain_file']):
                    self.rain_file = line.split()[1]
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


class BoundaryConditions(object):
    """
    """
    # definition of valid LISFLOOD-FP keywords
    valid_unit = ['days', 'hours', 'seconds']
    type_need_value = ['QFIX', 'QVAR', 'HFIX', 'HVAR']
    type_var_value = ['QVAR', 'HVAR']
    valid_type = ['CLOSED', 'FREE'].extend(type_need_value)
    valid_location = ['N', 'S', 'E', 'W', 'P', 'F']

    # Correspondance between LISFLOOD-FP boundary type and t.sim.flood
    # QFIX and QVAR are not registered as boundary conditions in t.sim.flood
    # but rather as user defined inflow in m/s
    bc_conv = {'CLOSED':1,
                'FREE':2,
                'HFIX':3,
                'HVAR':3}

    def __init__(self, msgr, sim_time, region):
        '''
        content: a dict as follow:
        {'bc_name':{'type':'HVAR', 'start_coor':(x,y), 'end_coor':(x,y)},
        'bc_name2':{'type':'CLOSED', 'start_coor':(x,y), 'end_coor':(x,y),
        'time_unit':'', 'values':[(,),(,)]}  # with data read from bdy
        ...}
        '''
        self.msgr = msgr
        self.region = region
        self.content = {}
        # lisflood sim_time in seconds, from .par file
        self.sim_time = sim_time
        self.mapset = set_mapset()

    def set_mapset(self):
        mapset = Mapset()
        return mapset.name

    def read_bci(self, bci_file):
        '''
        '''
        # Read the file and transform it to a list of lists
        file_content = []
        file_encoding = chardet.detect(open(bci_file, "r").read())
        file_open = codecs.open(bci_file,
                        encoding=file_encoding['encoding'])
        with file_open as input_file:
            for line in input_file:
                line = line.strip().split()
                if not line:
                    continue
                file_content.append(line)

        # read the file content and add it to the content dict
        for line_num, line in enumerate(file_content):
            # Check validity
            if not is_bci_line_valid(line):
                msgr.fatal(
                    'File {}, line {}: Invalid format'.format(
                    os.path.basename(bci_file), line_num))
            elif line[3] in self.type_need_value and not line[4]:
                msgr.fatal(
                    'File {}, line {}: Value needed'.format(
                    os.path.basename(bci_file), line_num))
            elif line[3] in self.type_var_value and (
                not idx_exist(line, 4) or is_number(line[4])):
                msgr.fatal(
                    'File {}, line {}: Time-series name needed'.format(
                    os.path.basename(bci_file), line_num))
            # Assign boundary conditions data
            elif line[3] in self.type_var_value:
                self.content[line[4]] = {'type': line[3],
                    'start_coor': self.get_array_coordinates(line)[0],
                    'end_coor': self.get_array_coordinates(line)[1]}
            elif line[3] not in self.type_var_value:
                self.content[str(line_num)] = {'type': line[3],
                    'start_coor': self.get_array_coordinates(line)[0],
                    'end_coor': self.get_array_coordinates(line)[1]
            # Add value for fixed boundary condition
            if (line[3] in self.type_need_value and
                        not in self.type_var_value):
                self.content[str(line_num)] = {
                   'time_unit':'seconds',
                    'values':[(line[4], 0), (line[4], self.sim_time)]}

        return self

    def get_array_coordinates(self, line):
        if line[0] == 'P':
            coord_geo1_row = float(line[2])
            coord_geo1_col = float(line[1])
            coord_geo2_row = float(line[2])
            coord_geo2_col = float(line[1])
        if line[0] == 'E':
            coord_geo1_row = float(line[2])
            coord_geo1_col = self.region.east
            coord_geo2_row = float(line[1])
            coord_geo2_col = self.region.east
        if line[0] == 'W':
            coord_geo1_row = float(line[2])
            coord_geo1_col = self.region.west
            coord_geo2_row = float(line[1])
            coord_geo2_col = self.region.west
        if line[0] == 'N':
            coord_geo1_row = self.region.north
            coord_geo1_col = float(line[2])
            coord_geo2_row = self.region.north
            coord_geo2_col = float(line[1])
        if line[0] == 'S':
            coord_geo1_row = self.region.south
            coord_geo1_col = float(line[2])
            coord_geo2_row = self.region.south
            coord_geo2_col = float(line[1])

        assert is_number(coord_geo1_row)
        assert is_number(coord_geo1_col)
        assert is_number(coord_geo2_row)
        assert is_number(coord_geo2_col)

        # crop coordinates to fit in region
        coord_geo1_row = max(coord_geo1_row, self.region.south)
        coord_geo1_col = max(coord_geo1_col, self.region.west)
        coord_geo2_row = min(coord_geo2_row, self.region.north)
        coord_geo2_col = min(coord_geo2_col, self.region.east)

        bc_len_row = coord_geo1_row - coord_geo2_row
        bc_len_col = coord_geo1_col - coord_geo2_col
        if bc_len_row < 0 or bc_len_col < 0:
             self.msgr.fatal(
                'Incoherent coordinates \n {}'.format(line))

        # transform into array coordinates
        coord_arr_1 = utils.coor2pixel(
                        (coord_geo1_col, coord_geo1_row), self.region)
        coord_arr_2 = utils.coor2pixel(
                        (coord_geo2_col, coord_geo2_row), self.region)
        return (coord_arr_1, coord_arr_2)

    def read_bdy(self, bdy_file):
        '''
        '''
        # Read the file and transform it to a list of lists
        file_content = []
        file_encoding = chardet.detect(open(bdy_file, "r").read())
        file_open = codecs.open(bdy_file,
                        encoding=file_encoding['encoding'])
        with file_open as input_file:
            for line_num, line in enumerate(input_file, 1):
                line = line.strip().split()
                if line_num == 1 or not line:
                    continue
                file_content.append(line)

        # read the file content and add it to the content dict
        for line_num, line in enumerate(file_content):
            # detect the start of a section
            if line[0] in self.content:
                # if key 'unit' already in content dict, it's likely that
                # the section has already been encountered
                if 'time_unit' in self.content[line[0]]:
                    msgr.fatal(
                        'Line {}: Boundary condition {} already read'.format(
                            line_num, line[0]))
                section_name = line[0]
                section_line_num = line_num
                current_section = self.content[section_name]

            # only if at least 2 elements in line
            if idx_exist(line, 1):
                # line after the section is supposed to hold the unit
                # the first element should be a number
                if line_num == section_line_num + 1 and is_number(line[0]):
                    # check validity of the unit
                    if line[1] not in self.valid_unit:
                        msgr.fatal(
                            'Line {}: unit {} unknown'.format(
                                line_num, line[1]))
                    current_section['time_unit'] = line[1]

                # if both elements are numbers, it's likely the values
                if is_number(line[0]) and is_number(line[1]):
                    if 'values' not in current_section:
                        current_section['values'] = [(
                            float(line[0]), float(line[1]))]
                    else:
                        current_section['values'].append((
                            float(line[0]), float(line[1])))

        return self

    def create_stds(self, stds_name, overwrite):
        stds_id = tgis.AbstractMapDataset.build_id(stds_name, self.mapset)
        stds_type="strds"
        temporal_type="relative"
        tgis.init()

        self.dbif = tgis.SQLDatabaseInterfaceConnection()
        self.dbif.connect()

        self.stds_h = tgis.open_new_stds(name=stds_id, type=stds_type,
                    temporaltype=temporal_type, title='', descr='',
                    semantic="mean", dbif=self.dbif,
                    overwrite=overwrite)
        return self

    def populate_user_flow_stds(self, rast_quser_name, overwrite):
        '''rast_quser_name: name of user flow raster map
        '''
        arr_qfix = grass.array.array(dtype=np.float32)
        arr_user_var = grass.array.array(dtype=np.float32)
        
        # fixed boundary condition
        

        map_list = []
        # iterate in variable boundary conditions
        
        for bc_key, bc_value in self.bcvar.iteritems():
            # iterate in content of a specific boundary condition
            for content in self.content[bc_key]:
                if bc_value['loc'] == 'P':
                    start_coord = bc_value['coord']
                    end_coord = bc_value['coord']
                    # set values in GRASS maps in m/s
                    arr_qvar = self.populate_array(
                            start_coord, end_coord,
                            value=(content[0] / self.region.ewres))
                if bc_value['loc'] == 'W':
                    start_coord = bc_value['coord']
                    end_coord = bc_value['coord']
                    arr_qvar = self.populate_array(
                            start_coord, end_coord,
                            value=(content[0] / self.region.ewres))

        # Create GRASS map including all QFIX and QVAR boundary conditions
        
        # write GRASS map
        rast_name_var = '{}_{}'.format(rast_user_name, str(content[1]))
        rast_id_var = tgis.AbstractMapDataset.build_id(
                        rast_name_var, self.mapset)
        # add temporal informations
        rast_var = tgis.RasterDataset(rast_id_var)
        rast_var.set_relative_time(start_time=content[1],
                    end_time=None, unit=self.unit[bc_key])
        map_list.append(rast_var)

        tgis.register.register_map_object_list('raster',
                        map_list, output_stds=self.stds_h,
                             delete_empty=True, unit=self.unit[bc_key],
                             dbif=self.dbif)
        return self


def is_number(s):
    try:
        float(s)
    except ValueError:
        return False
    else:
        return True

def is_bci_line_valid(line):
    if not idx_exist(line, 3):
        return False
    else:
        return (not is_number(line[0])
                and is_number(line[1])
                and is_number(line[2])
                and not is_number(line[3]))

def idx_exist(lst, idx):
    try:
        lst[idx]
    except IndexError:
        return False
    else:
        return True

def populate_array(start_coord, end_coord, value):
    start_row = start_coord[0]
    start_col = start_coord[1]
    end_row = end_coord[0]
    end_col = end_coord[1]
    assert start_row >= end_row
    assert start_col >= end_col

    # Make sure slices have at least one cell
    if start_row == end_row:
        row_slice = start_row
    else:
        row_slice = slice(start_row, end_row)
    if start_col == end_col:
        col_slice = start_col
    else:
        col_slice = slice(start_col, end_col)
    # value affectation
    arr = grass.array.array(dtype=np.float32)
    arr[row_slice, col_slice] = value
    return arr

if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())

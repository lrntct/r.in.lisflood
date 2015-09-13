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

import grass.script as grass
import grass.temporal as tgis
from grass.pygrass import raster
from grass.pygrass.gis.region import Region
from grass.pygrass.messages import Messenger
from grass.pygrass.modules import Module

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

    # Load *.par file
    par = Par(msgr, options['par_file'])
    par.read()

    # Write DEM file
    par.write_dem(rast_dem_name, grass.overwrite())
    # set rcomputational region to match DEM
    par.set_region_from_map(par.dem_file)
    # Write friction
    par.write_n_map(rast_n_file_name, grass.overwrite())
    # Write start file
    par.write_start_h(rast_start_file_name, grass.overwrite())

    # *.bci file
    if par.bci_file:
        bci_full_path = os.path.join(par.directory, par.bci_file)
        bci = Bci(msgr, bci_full_path)
        bci.read()
        print bci.content

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
        self.region.write()
        return self


class Bci(object):
    '''
    '''

    # valid key letters for boundary location
    klt = ['N', 'S', 'E', 'W', 'P', 'F']
    # valid boundary types
    bc_type = ['CLOSED', 'FREE',
                'HFIX', 'QFIX',
                'HFIX', 'HVAR']


    def __init__(self, msgr, bci_file):
        self.bci_file = bci_file  # full path to the file
        self.msgr = msgr
        self.content = []


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
                if not line or line[0] == '#':
                    continue
                if line[0] in self.klt:
                    self.content.append(line)
                if line[3] not in self.bc_type:
                    self.msgr.fatal(
                        'Unknown boundary type {} at line {}'.format(
                                            line[3], line_num))
        return self


if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())

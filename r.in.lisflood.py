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
    par_file = options['par_file']
    par_directory = os.path.dirname(par_file)
    par_encoding = chardet.detect(open(par_file, "r").read())
    par_file_open = codecs.open(par_file, encoding=par_encoding['encoding'])
    n_file = None
    friction = None
    with par_file_open as input_file:
        for line in input_file:
            # remove starting or trailing white spaces
            line = line.strip()
            if line.startswith(par_kwd['dem_file']):
                dem_file = line.split()[1]
            if line.startswith(par_kwd['bci_file']):
                bci_file = line.split()[1]
            if line.startswith(par_kwd['bdy_file']):
                bdy_file = line.split()[1]
            if line.startswith(par_kwd['start_file']):
                start_file = line.split()[1]
            if line.startswith(par_kwd['sim_time']):
                sim_time = line.split()[1]
            if line.startswith(par_kwd['n_file']):
                n_file = line.split()[1]
            if line.startswith(par_kwd['friction']):
                friction = line.split()[1]

    # load DEM file
    if dem_file:
        r_in_gdal = Module("r.in.gdal",
            input=os.path.join(par_directory, dem_file),
            output=rast_dem_name, overwrite=grass.overwrite())
        # get region from the DEM
        dem_region = Region()
        dem_region.from_rast(rast_dem_name)
    else:
        msgr.fatal('No {} found'.format(par_kwd['dem_file']))

    # friction
    if not n_file and not friction:
        msgr.fatal('Either {} or {} should be provided'.format(
                            par_kwd['n_file'], par_kwd['friction']))
    if friction and not n_file:
        write_n_map(friction, rast_n_file_name, dem_region)
    if n_file:
        r_in_gdal = Module("r.in.gdal",
            input=os.path.join(par_directory, n_file),
            output=rast_n_file_name, overwrite=grass.overwrite())

    # Start file
    if start_file:
        r_in_gdal = Module("r.in.gdal",
            input=os.path.join(par_directory, start_file),
            output=rast_start_file_name, overwrite=grass.overwrite())

    # *.bci file
    if bci_file:
        bci_full_path = os.path.join(par_directory, bci_file)
        bci_content = read_bci(bci_full_path, msgr)

    # Make sure the original region is restored
    region.write()
    return 0


# keywords of the *.par file
par_kwd = {'dem_file':'DEMfile',
            'bci_file':'bcifile',
            'bdy_file':'bdyfile',
            'start_file':'startfile',
            'friction':'fpfric',
            'n_file':'manningfile',
            'sim_time':'sim_time'}

# valid key letters for boundary location
bc_klt = ['N', 'S', 'E', 'W', 'P', 'F']
# valid boundary types
bc_type = ['CLOSED', 'FREE',
            'HFIX', 'QFIX',
            'HFIX', 'HVAR']


def write_n_map(friction, file_name, dem_region):
    '''write an uniform friction map from a given value
    '''
    # set the region
    dem_region.write()
    mapcalc_expression = '{map} = {value}'.format(map=file_name,
                                                value=friction)
    Module("r.mapcalc",
        expression=mapcalc_expression,
        overwrite=grass.overwrite())
    return 0


def read_bci(bci, msgr):
    '''Read a boundary condition *.bci file
    '''
    bci_file = bci
    bci_directory = os.path.dirname(bci_file)
    bci_encoding = chardet.detect(open(bci_file, "r").read())
    bci_file_open = codecs.open(bci_file, encoding=bci_encoding['encoding'])
    bci_content = []
    with bci_file_open as input_file:
        for line_num, line in enumerate(input_file, 1):
            # transform in a list without leading and trailing spaces
            line = line.strip().split()
            if not line or line[0] == '#':
                continue
            if line[0] in bc_klt:
                bci_content.append(line)
            if line[3] not in bc_type:
                msgr.fatal('Unknown boundary type {} at line {}'.format(
                                                line[3], line_num))
    return bci_content


if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())

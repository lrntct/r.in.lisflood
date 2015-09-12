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

    # Load *.par file
    par_file = options['par_file']
    par_directory = os.path.dirname(par_file)
    par_encoding = chardet.detect(open(par_file, "r").read())
    rast_n_file_name = options['friction']
    rast_dem_name = options['dem']
    n_file = None
    friction = None
    with codecs.open(par_file, encoding=par_encoding['encoding']) as input_file:
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


if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())

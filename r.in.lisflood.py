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

import grass.script as grass
import grass.temporal as tgis
from grass.pygrass import raster
from grass.pygrass.gis.region import Region
from grass.pygrass.messages import Messenger

import chardet
import codecs
import sys
import os

def main():
    par_file = options['par_file']
    par_directory = os.path.dirname(par_file)
    par_encoding = chardet.detect(open(par_file, "r").read())
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
    print dem_file, bci_file, bdy_file, start_file
    return 0
 
# keywords of the *.par file
par_kwd = {'dem_file':'DEMfile',
            'bci_file':'bcifile',
            'bdy_file':'bdyfile',
            'start_file':'startfile',
            'sim_time':'sim_time'}

if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())

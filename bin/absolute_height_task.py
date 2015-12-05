#!/usr/bin/env python
import argparse
from absoluteHeightTask import absoluteHeightTask

parser = argparse.ArgumentParser(description='Perform absolute height analysis on OGP data.')
parser.add_argument('ogp_file', help='OGP data file of xyz tuples containing gauge block and sensor metrology points')
parser.add_argument('sensor_id', help='LSST ID number of sensor')
parser.add_argument('--zoffset', type=float, default=-1.2,
                    help='z-offset')

args = parser.parse_args()

absoluteHeightTask(args.sensor_id, args.ogp_file, zoffset=args.zoffset)

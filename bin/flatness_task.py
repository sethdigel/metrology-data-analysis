#!/usr/bin/env python
import argparse
from flatnessTask import flatnessTask

parser = argparse.ArgumentParser(description='Sensor flatness analysis.')
parser.add_argument('infile', help='Data file of xyz tuples containing sensor metrology points')
parser.add_argument('sensor_id', help='LSST ID number of sensor')
parser.add_argument('--datatype', type=str, default='OGP',
                    help='Formatting of input file: OGP, e2v, ITL')

args = parser.parse_args()

flatnessTask(args.sensor_id, args.infile, dtype=args.datatype)

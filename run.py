#!/usr/bin/env python3

import subprocess
import argparse

# parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("--ps", type=int)
ap.add_argument("--t", type=float)
ap.add_argument("--o", type=int)
ap.add_argument("--run", type=int)
args = ap.parse_args()

#run original Efros-Leung synthesis
if not args.run:
    p = subprocess.run(['efros_leung', 'input_0.png', '-p', str(args.ps), '-t',  str(args.t), '-s',  str(args.o), 'output.png','map.png',
                            'copy_map.png'])
#run accelerated Efros-Leung synthesis
if args.run:
    p = subprocess.run(['efrosLeung_accel', 'input_0.png', '-p', str(args.ps), '-t', str(args.t), '-s', str(args.o), \
            'output.png', 'map.png', 'copy_map.png'])
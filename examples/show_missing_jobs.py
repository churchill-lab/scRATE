#!/usr/bin/env python
import argparse
import os
import sys
import glob

def main():
    # Organize user input
    parser = argparse.ArgumentParser()
    parser.add_argument('inbase', action='store', type=str, help='basename of input files')
    parser.add_argument('outbase', action='store', type=str, help='basename of output files')
    parser.add_argument('--csv', dest='csv_wanted', action='store_true', help='print the Job IDs comma separated')
    args = parser.parse_args()

    inlist  = glob.glob(args.inbase + '*')
    outlist = glob.glob(args.outbase + '*')
    inset  = set([ os.path.basename(i).split('.')[1] for i in inlist ])
    outset = set([ os.path.basename(i).split('.')[1] for i in outlist ])
    difflist = list(inset - outset)
    difflist.sort()
    if args.csv_wanted:
        print(','.join([ i.lstrip('0') for i in difflist ]))
    else:
        for i in difflist:
            print('{}.{}'.format(args.inbase, i))
    print

    return 0

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover

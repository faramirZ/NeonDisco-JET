#!/usr/bin/env python
"""
RNA Editing Site Filter Script

Usage:
    python script.py [OPTIONS]

Options:
    --prefix PREFIX          Prefix for directory names (default: '')
    --known FILE             Path to knownEditing file (default: 'knownEditing')
    --pos FILE               Path to pos.txt file (default: 'pos.txt')
    --posalu FILE            Path to posalu.txt file (default: 'posalu.txt')
    --output FILE            Output file name (default: 'editing.txt')
    -h, --help               Show this help message
"""

import sys
import os
import glob
import argparse

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Filter RNA editing sites based on Alu element regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('--prefix', type=str, default='',
                        help='Prefix for directory names (default: "")')
    parser.add_argument('--known', type=str, default='knownEditing',
                        help='Path to knownEditing file (default: "knownEditing")')
    parser.add_argument('--pos', type=str, default='pos.txt',
                        help='Path to pos.txt file (default: "pos.txt")')
    parser.add_argument('--posalu', type=str, default='posalu.txt',
                        help='Path to posalu.txt file (default: "posalu.txt")')
    parser.add_argument('--output', type=str, default='editing.txt',
                        help='Output file name (default: "editing.txt")')
    
    args = parser.parse_args()
    
    # Build directory paths with prefix
    prefix = args.prefix
    if prefix and not prefix.endswith('-'):
        prefix += '-'
    
    firstalu_pattern = '%sfirstalu/DnaRna_*/outTable_*' % prefix
    second_pattern = '%ssecond/DnaRna_*/outTable_*' % prefix
    
    # Find table files using glob
    atab_matches = glob.glob(firstalu_pattern)
    ftab_matches = glob.glob(second_pattern)
    
    if not atab_matches:
        sys.exit('No files found matching pattern: %s' % firstalu_pattern)
    
    atab = atab_matches[0]  # alu refined
    ftab = ftab_matches[0] if ftab_matches else None
    
    # Check required input files exist
    if not os.path.exists(args.known):
        sys.exit('%s file not found.' % args.known)
    if not os.path.exists(args.pos):
        sys.exit('%s file not found.' % args.pos)
    if not os.path.exists(args.posalu):
        sys.exit('%s file not found.' % args.posalu)
    
    # Main algorithm (unchanged)
    o = open(args.output, 'w')
    f = open(args.known)
    for i in f:
        o.write(i)
    f.close()
    
    if ftab and os.path.exists(ftab):
        f = open(ftab)
        d = {}
        for i in f:
            if i.startswith('Region'):
                continue
            l = (i.strip()).split('\t')
            d[(l[0], l[1])] = 0
        f.close()
        f = open(args.pos)
        for i in f:
            if i.startswith('Region'):
                continue
            l = (i.strip()).split('\t')
            if d.has_key((l[0], l[1])):
                o.write(i)
        f.close()
    
    f = open(atab)
    d = {}
    for i in f:
        if i.startswith('Region'):
            continue
        l = (i.strip()).split('\t')
        d[(l[0], l[1])] = 0
    f.close()
    f = open(args.posalu)
    for i in f:
        if i.startswith('Region'):
            continue
        l = (i.strip()).split('\t')
        if d.has_key((l[0], l[1])):
            o.write(i)
    f.close()
    o.close()
    
    print 'Successfully created %s' % args.output
    print 'Used firstalu table: %s' % atab
    if ftab:
        print 'Used second table: %s' % ftab

if __name__ == '__main__':
    main()

    
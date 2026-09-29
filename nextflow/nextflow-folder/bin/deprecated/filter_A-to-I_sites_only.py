#!/usr/bin/env python
"""
Filter A-to-I (AG) Editing Sites

Extracts only ADAR-mediated A-to-I editing sites (appearing as AG substitutions)
from REDItools editing output and adds a descriptive header.

Usage:
    python filterAG.py [OPTIONS]

Options:
    --input FILE             Input editing file (default: 'editing.txt')
    --output FILE            Output AG-only file (default: 'editing_AG_only.txt')
    -h, --help               Show this help message
"""

import sys
import os
import argparse


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Filter for A-to-I (AG) editing sites only',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('--input', type=str, default='editing.txt',
                        help='Input editing file (default: "editing.txt")')
    parser.add_argument('--output', type=str, default='editing_AG_only.txt',
                        help='Output AG-filtered file (default: "editing_AG_only.txt")')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not os.path.exists(args.input):
        sys.exit('%s file not found.' % args.input)
    
    # Define REDItools column headers
    header = [
        'Region',
        'Position',
        'Reference',
        'Strand',
        'Coverage-q',
        'MeanQ',
        'BaseCount[A,C,G,T]',
        'AllSubs',
        'Frequency',
        'gCoverage-q',
        'gMeanQ',
        'gBaseCount[A,C,G,T]',
        'gAllSubs',
        'gFrequency',
        'RepeatType',
        'RepeatName',
        'SNPFlag',
        'dbSNP_ID',
        'REDIPortalKnownEditingSites'
    ]
    
    # Filter for AG sites
    print 'Reading %s...' % args.input
    total_sites = 0
    ag_sites = 0
    
    f_out = open(args.output, 'w')
    
    # Write header
    f_out.write('\t'.join(header) + '\n')
    
    # Process input file
    f_in = open(args.input)
    for line in f_in:
        line = line.strip()
        if not line:
            continue
        
        total_sites += 1
        cols = line.split('\t')
        
        # Column 7 (index 7) contains AllSubs
        if len(cols) > 7 and cols[7] == 'AG':
            ag_sites += 1
            f_out.write(line + '\n')
    
    f_in.close()
    f_out.close()
    
    print 'Total editing sites: %d' % total_sites
    print 'A-to-I (AG) sites: %d (%.1f%%)' % (ag_sites, (ag_sites / float(total_sites)) * 100)
    print 'Output saved to: %s' % args.output


if __name__ == '__main__':
    main()

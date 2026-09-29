#!/usr/bin/env python
"""
RNA Editing Statistics Calculator

Calculates substitution type distributions across different repeat categories
(ALU, non-ALU repeats, non-repetitive regions) from RNA editing data.

Usage:
    python getStatistics.py [OPTIONS]

Options:
    --input FILE             Input editing file (default: 'editing.txt')
    --output FILE            Output statistics file (default: 'editingStats.txt')
    -h, --help               Show this help message
"""

import sys
import os
import ast
import argparse


def getDistro(lines):
    """Calculate substitution type distribution from editing lines.
    
    Args:
        lines: List of parsed editing data lines (split by tabs)
        
    Returns:
        Dictionary mapping substitution types (e.g. 'AG') to tuple of
        (count, total, percentage)
    """
    # Initialize all possible substitution types
    s = {}
    for i in 'ACGT':
        for j in 'ACGT':
            if i != j:
                s[i + j] = 0
    
    # Nucleotide to index mapping
    n = {}
    x = 0
    for i in 'ACGT':
        n[i] = x
        x += 1
    
    # Count substitutions
    all = 0
    for i in lines:
        sub = i[7].split()[0]  # Substitution type (e.g. 'AG')
        nuc = ast.literal_eval(i[6])  # Nucleotide counts array (safe evaluation)
        nv = nuc[n[sub[1]]]    # Count for the variant nucleotide
        s[sub] += nv
        all += nv
    
    # Calculate percentages
    d = {}
    for i in s:
        try:
            v = (s[i] / float(all)) * 100
        except:
            v = 0.0
        d[i] = (s[i], all, v)
    
    return d


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Calculate RNA editing statistics by repeat category',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('--input', type=str, default='editing.txt',
                        help='Input editing file (default: "editing.txt")')
    parser.add_argument('--output', type=str, default='editingStats.txt',
                        help='Output statistics file (default: "editingStats.txt")')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not os.path.exists(args.input):
        sys.exit('%s file not found.' % args.input)
    
    # Parse editing data into categories
    alu, nonalu, nonrep, kn = [], [], [], 0
    
    print 'Reading %s...' % args.input
    f = open(args.input)
    for i in f:
        if i.startswith('Reg'):
            continue
        l = (i.strip()).split('\t')
        
        # Count known editing sites
        if l[18] == 'ed':
            kn += 1
        
        # Categorize by repeat type
        if l[14] == 'SINE' and l[15][:3] == 'Alu':
            alu.append(l)
        elif l[14] != '-' and l[15][:3] != 'Alu':
            nonalu.append(l)
        elif l[14] == '-' and l[15] == '-':
            nonrep.append(l)
    f.close()
    
    print 'Found %d ALU sites' % len(alu)
    print 'Found %d non-ALU repeat sites' % len(nonalu)
    print 'Found %d non-repetitive sites' % len(nonrep)
    print 'Found %d known editing sites' % kn
    
    # Calculate distributions
    print 'Calculating substitution distributions...'
    alust = getDistro(alu)
    nonalust = getDistro(nonalu)
    nonrepst = getDistro(nonrep)
    all = getDistro(alu + nonalu + nonrep)
    
    # Write output
    print 'Writing statistics to %s...' % args.output
    f = open(args.output, 'w')
    h = ['SubType', 'ALU', 'REPnonALU', 'NONREP', 'ALL']
    f.write('\t'.join(h) + '\n')
    
    for i in alust:
        r = [i, alust[i][2], nonalust[i][2], nonrepst[i][2], all[i][2]]
        r = [str(x) for x in r]
        f.write('\t'.join(r) + '\n')
    f.close()
    
    print 'Done! Statistics saved to %s' % args.output


if __name__ == '__main__':
    main()

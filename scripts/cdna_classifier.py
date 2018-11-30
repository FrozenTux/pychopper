#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from Bio import SeqIO
import os
import tqdm
from pychopper import seq_utils as seu
from pychopper import chopper
from pychopper import report

"""
Parse command line arguments.
"""
parser = argparse.ArgumentParser(
    description='Tool to identify full length cDNA reads. Primers have to specified as they are on the forward strand.!!!!!DEVELOPMENT VERSION!!!!')
parser.add_argument(
    '-b', metavar='barcodes', type=str, default=None, help="Primers fasta.", required=True)
parser.add_argument(
    '-i', metavar='input_format', type=str, default='fastq', help="Input/output format (fastq).")
parser.add_argument('-g', metavar='aln_params', type=str,
                    help="Alignment parameters (match, mismatch,gap_open,gap_extend).", default="1,-1,2,1")
parser.add_argument(
    '-t', metavar='target_length', type=int, default=200, help="Number of bases to scan at each end (200).")
parser.add_argument(
    '-s', metavar='score_percentile', type=float, default=98, help="Score cutoff percentile (98).")
parser.add_argument(
    '-n', metavar='sample_size', type=int, default=100000, help="Number of samples when calculating score cutoff (100000).")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, default=None, help="Report PDF.")
parser.add_argument(
    '-u', metavar='unclass_output', type=str, default=None, help="Write unclassified reads to this file.")
parser.add_argument('-a', metavar='aln_report', type=str, default=None, help="Alignment report file")
parser.add_argument('input_fastx', metavar='input_fastx', type=str, help="Input file.")
parser.add_argument('output_fastx', metavar='output_fastx', type=str, help="Output file.")


def _revcomp_seq(read):
    """ Reverse complement sequence and fix Id. """
    rev_read = read.reverse_complement()
    rev_read.id, rev_read.description = read.id, read.description
    return rev_read


def _filter_and_annotate(read, match):
    """ Filter sequences by match and annotate with direction. """
    if match is not None:
        direction = '+' if match == 'fwd_match' else '-'
        read.id = read.id + "|" + direction
        if match == 'rev_match':
            read = _revcomp_seq(read)
        return read, True
    else:
        return read, False


def _parse_aln_params(pstr):
    """ Parse alignment parameters. """
    res = {}
    tmp = [int(x) for x in pstr.split(',')]
    res['match'] = tmp[0]
    res['mismatch'] = tmp[1]
    res['gap_open'] = tmp[2]
    res['gap_extend'] = tmp[3]
    return res


def _record_size(read, in_format):
    """ Calculate record size. """
    dl = len(read.description)
    sl = len(read.seq)
    if in_format == 'fastq':
        bl = dl + 2 * sl + 6
    elif in_format == 'fasta':
        bl = dl + sl + 3
    else:
        raise Exception("Unkonwn format!")
    return bl


if __name__ == '__main__':
    args = parser.parse_args()

    ALIGN_PARAMS = _parse_aln_params(args.g)

    barcodes = chopper.load_barcodes(args.b)
    barcodes = chopper.calculate_score_cutoffs(
        barcodes, aln_params=ALIGN_PARAMS, target_length=args.t, percentile=args.s, nr_samples=args.n)

    output_handle = open(args.output_fastx, "w")
    if args.u is not None:
        unclass_handle = open(args.u, "w")
    if args.a:
        aln_handle = open(args.a, 'w')

    unclass_nr_hits = []
    fwd_matches = 0
    rev_matches = 0
    unclassified = 0

    # Get the size of input file:
    input_size = os.stat(args.input_fastx).st_size

    pbar = tqdm.tqdm(total=input_size)

    for read in seu.read_seq_records(args.input_fastx, args.i):
        pbar.update(_record_size(read, args.i))
        match, nr_hits, alns = list(chopper.score_barcode_groups(read, barcodes, args.t, ALIGN_PARAMS).values())[0]
        if match is not None:
            if match == 'fwd_match':
                fwd_matches += 1
            if match == 'rev_match':
                rev_matches += 1

            if args.a:
                # Compute values for alignment report

                # Write lines
                aln_line = read.id + "\t" # 1 : read id
                aln_line += match.split("_")[0] # 2: fwd or rev
                # TODO

        read, match = _filter_and_annotate(read, match)

        if match is True:
            SeqIO.write(read, output_handle, args.i)
        else:
            unclassified += 1
            if args.r is not None:
                unclass_nr_hits.append(nr_hits)
            if args.u is not None:
                SeqIO.write(read, unclass_handle, args.i)

    output_handle.flush()
    output_handle.close()

    if args.u is not None:
        unclass_handle.flush()
        unclass_handle.close()
    if args.a:
        aln_handle.flush()
        aln_handle.close()

    if args.r is not None:
        plotter = report.Report(args.r)
        plotter.plot_bars_simple({'Classified': fwd_matches + rev_matches, 'Unclassified': unclassified}, title="Basic statistics", ylab="Count")
        plotter.plot_histograms({'nr_hits': unclass_nr_hits}, title="Number of hits in unclassified reads", xlab="Number of hits", ylab="Count")
        plotter.plot_bars_simple({'+': fwd_matches, '-': rev_matches}, title="Strandedness of classified reads", xlab="Strand", ylab="Count")
        plotter.close()

    pbar.close()

#!/usr/bin/env python

import csv
import os
import sys
import gzip
from collections import defaultdict
import argparse

# See http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

parser = argparse.ArgumentParser(description ="""This script parses input
                gff3 file and generates a table of introns""")
parser.add_argument('-i', '--infile', help="input file; default STDIN")
parser.add_argument('-o', '--outfile', help="output file; default STDOUT")
parser.add_argument('-z', '--gzip', help="input is gzip compressed", \
                    action = 'store_true')
args = parser.parse_args()

if args.infile:
    if args.gzip:
        gff3_file = gzip.open(args.infile, 'rt')
    else:
        gff3_file = open(args.infile, 'r')
else:
    gff3_file = sys.stdin

if args.outfile:
    introns_file = open(args.outfile, 'w')
else:
    introns_file = sys.stdout

def get_exons(gff3_file):
    tbl = csv.reader(gff3_file, delimiter = '\t')
    exons_dict = defaultdict(list)
    for line in tbl:
        if not line[0].startswith('#'):
            [
                chrom,
                feat_source,
                feat_type,
                start,
                stop,
                score,
                strand,
                phase,
                attribs
            ] = line
            if (feat_type == "exon"
                and 'transcript_id=' in attribs):
                new_attribs = process_attribs(attribs)
                #print(new_attribs)
                if 'gene_id' in new_attribs:
                    gene_id = new_attribs['gene_name']
                elif 'gene' in new_attribs:
                    gene_id = new_attribs['gene']
                else:
                    gene_id = 'UNKNOWN'
                tx = new_attribs['transcript_id']
                start, stop = int(start), int(stop)
                exons_dict[(chrom, strand, gene_id, tx)].append((start, stop))
            elif(feat_type == "exon"
                 and 'Parent=transcript' in attribs):
                attribs = attribs.split(';')
                new_attribs = {}
                for attrib in attribs:
                    attrib = attrib.split('=')
                    new_attribs[attrib[0]] = attrib[1]
                gene_id = new_attribs['Parent'].strip('transcript:')
                tx = new_attribs['Parent'].strip('transcript:')
                start, stop = int(start), int(stop)
                exons_dict[(chrom, strand, gene_id, tx)].append((start, stop))
    gff3_file.close()
    return exons_dict

def process_attribs(attribs):
    new_attribs = {}
    attribs = list(filter(None, attribs.split(';'))) ## removes empty strings, needed because some gff3 lines have ";;"
    for attrib in attribs:
        k, v = attrib.split('=')
        if k == 'Dbxref':
            xrefs = v.split(',')
            for xref in xrefs:
                terms = xref.split(':')
                new_attribs[terms[-2]] = terms[-1]
        else:
            new_attribs[k] = v
    return new_attribs

def get_introns(exons_dict):
    introns_dict = defaultdict(set)
    for tx_info, exons in exons_dict.items():
        if len(exons) > 1:
            [
                chrom,
                strand,
                gene_id,
                tx
            ] = tx_info
            exons = sorted(exons)
            i = 0
            while (i + 1) < len(exons):
                introns = (exons[i][1] + 1, exons[i+1][0] - 1)
                i = i + 1
                introns_dict[(chrom, strand, gene_id, tx)].add(introns)
    return introns_dict

def tabulate_introns(introns_dict, introns_file):
    tbl = csv.writer(introns_file,
                     delimiter = '\t',
                     lineterminator = os.linesep)
    tbl.writerow(['#chrom',
                  'intron_start',
                  'intron_end',
                  'strand',
                  'gene_id',
                  'tx_acc',
                  'intron_num',
                  'intron_ct'])
    for tx_info, introns in introns_dict.items():
        [
            chrom,
            strand,
            gene_id,
            tx_acc
        ] = tx_info
        if strand == '+':
            introns = sorted(introns)
        elif strand == '-':
            introns = sorted(introns, reverse = True)
        num_introns = len(introns)
        for intron in introns:
            tbl.writerow([chrom,
                         intron[0],
                         intron[1],
                         strand,
                         gene_id,
                         tx_acc,
                         introns.index(intron) + 1,
                         num_introns])
    introns_file.close()

exons_dict = get_exons(gff3_file)
introns_dict = get_introns(exons_dict)
tabulate_introns(introns_dict, introns_file)


'''
{'ID': 'exon:ENSMUST00000192146.2:9', 'Parent': 'ENSMUST00000192146.2', 'gene_id': 'ENSMUSG00000091476.9', 'transcript_id': 'ENSMUST00000192146.2', 'gene_type': 'protein_coding', 'gene_name': 'Catspere2', 'transcript_type': 'nonsense_mediated_decay', 'transcript_name': 'Catspere2-203', 'exon_number': '9', 'exon_id': 'ENSMUSE00001336495.2', 'level': '2', 'protein_id': 'ENSMUSP00000142187.2', 'transcript_support_level': '1', 'mgi_id': 'MGI:5589632', 'havana_gene': 'OTTMUSG00000050375.4', 'havana_transcript': 'OTTMUST00000127854.3'}
{'ID': 'exon:ENSMUST00000192146.2:10', 'Parent': 'ENSMUST00000192146.2', 'gene_id': 'ENSMUSG00000091476.9', 'transcript_id': 'ENSMUST00000192146.2', 'gene_type': 'protein_coding', 'gene_name': 'Catspere2', 'transcript_type': 'nonsense_mediated_decay', 'transcript_name': 'Catspere2-203', 'exon_number': '10', 'exon_id': 'ENSMUSE00001345739.2', 'level': '2', 'protein_id': 'ENSMUSP00000142187.2', 'transcript_support_level': '1', 'mgi_id': 'MGI:5589632', 'havana_gene': 'OTTMUSG00000050375.4', 'havana_transcript': 'OTTMUST00000127854.3'}
{'ID': 'exon:ENSMUST00000192146.2:11', 'Parent': 'ENSMUST00000192146.2', 'gene_id': 'ENSMUSG00000091476.9', 'transcript_id': 'ENSMUST00000192146.2', 'gene_type': 'protein_coding', 'gene_name': 'Catspere2', 'transcript_type': 'nonsense_mediated_decay', 'transcript_name': 'Catspere2-203', 'exon_number': '11', 'exon_id': 'ENSMUSE00001338542.2', 'level': '2', 'protein_id': 'ENSMUSP00000142187.2', 'transcript_support_level': '1', 'mgi_id': 'MGI:5589632', 'havana_gene': 'OTTMUSG00000050375.4', 'havana_transcript': 'OTTMUST00000127854.3'}
{'ID': 'exon:ENSMUST00000192146.2:12', 'Parent': 'ENSMUST00000192146.2', 'gene_id': 'ENSMUSG00000091476.9', 'transcript_id': 'ENSMUST00000192146.2', 'gene_type': 'protein_coding', 'gene_name': 'Catspere2', 'transcript_type': 'nonsense_mediated_decay', 'transcript_name': 'Catspere2-203', 'exon_number': '12', 'exon_id': 'ENSMUSE00001340383.2', 'level': '2', 'protein_id': 'ENSMUSP00000142187.2', 'transcript_support_level': '1', 'mgi_id': 'MGI:5589632', 'havana_gene': 'OTTMUSG00000050375.4', 'havana_transcript': 'OTTMUST00000127854.3'}
{'ID': 'exon:ENSMUST00000192146.2:13', 'Parent': 'ENSMUST00000192146.2', 'gene_id': 'ENSMUSG00000091476.9', 'transcript_id': 'ENSMUST00000192146.2', 'gene_type': 'protein_coding', 'gene_name': 'Catspere2', 'transcript_type': 'nonsense_mediated_decay', 'transcript_name': 'Catspere2-203', 'exon_number': '13', 'exon_id': 'ENSMUSE00001342144.2', 'level': '2', 'protein_id': 'ENSMUSP00000142187.2', 'transcript_support_level': '1', 'mgi_id': 'MGI:5589632', 'havana_gene': 'OTTMUSG00000050375.4', 'havana_transcript': 'OTTMUST00000127854.3'}
{'ID': 'exon:ENSMUST00000192146.2:14', 'Parent': 'ENSMUST00000192146.2', 'gene_id': 'ENSMUSG00000091476.9', 'transcript_id': 'ENSMUST00000192146.2', 'gene_type': 'protein_coding', 'gene_name': 'Catspere2', 'transcript_type': 'nonsense_mediated_decay', 'transcript_name': 'Catspere2-203', 'exon_number': '14', 'exon_id': 'ENSMUSE00001339204.2', 'level': '2', 'protein_id': 'ENSMUSP00000142187.2', 'transcript_support_level': '1', 'mgi_id': 'MGI:5589632', 'havana_gene': 'OTTMUSG00000050375.4', 'havana_transcript': 'OTTMUST00000127854.3'}
'''

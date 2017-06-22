import sys
import os
import argparse
import json
import re

program = "tax2format"
version = "0.1.0"
author = "Darcy Jones, Evan Krell"
date = "21 June 2017"
email = "darcy.ab.jones@gmail.com, evan.krell@tamucc.edu"
blurb = (
    "{program}\n"
    "{:<}"
    )
short_blurb = (
    "Formats the json output of gi2tax.py as FASTA headers compatable with Munch's Statistical Assignmnet Package (SAP): github.com/kaspermunch/sap"
    )
license = (
    '{program} {version}\n'
    '{short_blurb}\n\n'
    'Copyright (C) {date},  {author}'
    '\n\n'
    'This program is free software: you can redistribute it and/or modify '
    'it under the terms of the GNU General Public License as published by '
    'the Free Software Foundation, either version 3 of the License, or '
    '(at your option) any later version.'
    '\n\n'
    'This program is distributed in the hope that it will be useful, '
    'but WITHOUT ANY WARRANTY; without even the implied warranty of '
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the '
    'GNU General Public License for more details.'
    '\n\n'
    'You should have received a copy of the GNU General Public License '
    'along with this program. If not, see <http://www.gnu.org/licenses/>.'
    ).format(**locals())


# Taxonomic ranks accepted by SAP
accepted_taxon_ranks = ['superkingdom', 'kingdom', 'subkingdom',
             'superphylum', 'phylum', 'subphylum',
             'superclass', 'class', 'subclass', 'infraclass',
             'superorder', 'order', 'suborder', 'infraorder', 'parvorder',
             'superfamily', 'family', 'subfamily',
             'supertribe', 'tribe', 'subtribe',
             'supergenus', 'genus', 'subgenus',
             'species group', 'species subgroup',
             'species', 'subspecies', 'varietas',
             'otu']

def main(in_file, out_file):
    
    output = None
    if out_file is not None:
        output = open(out_file, 'w')
    
    for line in in_file:
        try:
            record = json.loads (line)
            tax_path = record.pop ('tax_path')
            gi = str(record.pop ('gi'))
            id = record.pop ('id')

            rank_str = ">" + gi + " ; "
            for node in tax_path:
                if node['rank'] in accepted_taxon_ranks:
                    rank_str = rank_str + node['rank'].replace(",", "_")  + ": " + (node['rank_name']) + ", "

            rank_str = rank_str[0:-2]
            rank_str = rank_str + " ; " + id
            rank_str = re.sub (r'<.*>', "", rank_str)
            if output is not None:
                output.write(rank_str)
            else: # STDOUT
                sys.stdout.write(rank_str + '\n')
        except Exception:
            pass

    if output is not None:
        output.close()

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser ()
    arg_parser.add_argument(
        "-i", "--infile",
        dest='in_file',
        default=sys.stdin,
        type=argparse.FileType('r'),
        help=(
            "Path to input file containing taxonomy "
            "output by gi2tax.py. "
            "Enter '-' for stdin (default)."
        )
    )
    arg_parser.add_argument(
        "-o", "--outfile",
        dest="out_file",
        default=None,
        help=( 
            "Path to write output to. "
            "Default is stdout."
            )
        )
    args = arg_parser.parse_args ()
    main (**args.__dict__)

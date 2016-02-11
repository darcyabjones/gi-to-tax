#! /usr/bin/env python3

from __future__ import print_function
import sys
import argparse
import re
from collections import defaultdict
import os
import json
import sqlite3
import random
import hashlib

from Bio import SeqIO, Phylo, Entrez
from Bio.Blast import NCBIXML
from Bio.Phylo.PhyloXML import Phylogeny
from Bio.Phylo.PhyloXML import Clade
from Bio.Phylo.PhyloXML import Taxonomy
from Bio.Phylo.PhyloXML import Id
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.BaseTree import Clade as Clade2

program = "tax2format"
version = "0.1.0"
author = "Darcy Jones"
date = "14 January 2016"
email = "darcy.ab.jones@gmail.com"
short_blurb = (
    "Remove duplicate sequences from a sequence file."
    )
license = (
    '{program}-{version}\n'
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
    )

license = license.format(**locals())


"################################# Globals ##################################"

ORDER = [
    "code",
    "id",
    "gi",
    "name",
    "taxid",
    "description",
    "subspecies",
    "species",
    "subgenus",
    "genus",
    "subfamily",
    "family",
    "superfamily",
    "suborder",
    "order",
    "superorder",
    "subclass",
    "class",
    "superclass",
    "subphylum",
    "phylum",
    "superphylum",
    "subkingdom",
    "kingdom",
    "superkingdom",
    "no_rank",
    ]

"################################# Functions ################################"


def phyloxml_tree_generator(tree):
    children = phyloxml_tree_generator_recurse(tree, taxid=1)
    return Phylogeny(name='gi2tax - Common tree', root=children)


def phyloxml_tree_generator_recurse(tree, taxid=1):
    """ . """
    children = list()
    for child in tree[taxid]['children']:
        children.append(phyloxml_tree_generator_recurse(tree, child))
    if len(children) == 0:
        children = None

    if taxid == 1:
        rank_name = None
        node_id = 1
        tax_data = None
    else:
        if 'rank' in tree[taxid]:
            rank = tree[taxid]['rank']
            if rank not in Phylo.PhyloXML.Taxonomy.ok_rank:
                rank = "other"
        else:
            rank = None

        if children is None:
            rank_name = taxid
        else:
            rank_name = None
            node_id = taxid
    this_node = Clade(name=rank_name, clades=children, branch_length=0.1)
    if taxid != 1:
        this_node.taxonomy = Taxonomy(
            scientific_name=tree[taxid]['rank_name'],
            rank=rank,
            id=Id(taxid, provider='ncbi_taxonomy'),
            )
    return this_node


def newick_tree_generator(tree):
    children = newick_tree_generator_recurse(tree, taxid=1)
    return Tree(name='gi2tax - Common tree', root=children)


def newick_tree_generator_recurse(tree, taxid=1):
    """ . """
    children = list()

    for child in tree[taxid]['children']:
        children.append(newick_tree_generator_recurse(tree, child))

    # Drops nodes with only one child. eg genus with only one species.
    if len(children) == 1:
        return children[0]
    else:
        if taxid == 1:
            rank_name = None
        else:
            if len(children) == 0:
                children = None
                rank_name = str(tree[taxid]['rank_name'])
                # rank_name = str(taxid)  # might be more useful?
                # Remove newick illegal characters from name
                # One or more whitespace characters
                regex = re.compile("[\s,]+")
                rank_name = re.sub(regex, '_', rank_name)
                rank_name.replace('(', '[')
                rank_name.replace(')', ']')
            else:
                rank_name = None

        this_node = Clade2(name=rank_name, clades=children, branch_length=1.0)
        return this_node


def build_tree(record, tree):
    for node in record['tax_path']:
        if node['taxid'] in tree:
            if 'rank' not in tree[node['taxid']]:
                tree[node['taxid']]['rank'] = node['rank']

            if 'rank_name' not in tree[node['taxid']]:
                tree[node['taxid']]['rank_name'] = node['rank_name']

            if 'children' not in tree[node['taxid']]:
                tree[node['taxid']]['children'] = set()

        else:
            tree[node['taxid']] = {
                'rank': node['rank'],
                'rank_name': node['rank_name'],
                'children': set()
                }

        if node['parent_taxid'] in tree:
            tree[node['parent_taxid']]['children'].add(node['taxid'])

        else:
            tree[node['parent_taxid']] = {'children': {node['taxid']}}

    if 'gi' in record:
        if record['taxid'] in tree:
            tree[record['taxid']]['gi'] = record['gi']

        else:
            tree[record['taxid']] = {'gi': record['gi']}

    else:
        if record['taxid'] in tree:
            tree[record['taxid']]['id'] = record['id']

        else:
            tree[record['taxid']] = {'id': record['id']}

    return tree


def table_writer(record, sep="\t", order=ORDER):
    """ . """

    template = sep.join(["{{{}}}".format(c) for c in order])
    for col in order:
        if col not in record:
            record[col] = ""

    return template.format(**record) + "\n"


"################################### Main ###################################"


def main(
        in_file,
        out_file,
        out_format
        ):
    # Only required for tree formats
    tree = dict()

    if out_format == "tsv":
        sep = "\t"

    else:
        sep = ","

    if out_format in {'tsv', 'csv'}:
        columns = dict(zip(ORDER, ORDER))
        out_file.write(table_writer(columns, sep))

    for line in in_file:
        record = json.loads(line)

        if out_format in ('newick', 'phyloxml'):
            build_tree(record, tree)

        else:
            tax_path = record.pop('tax_path')

            no_rank = []
            for node in tax_path:
                rank = node['rank']
                rank_name = node['rank_name']

                if rank != 'no rank':
                    record['{}_taxid'.format(rank)] = node['taxid']
                    record[rank] = rank_name

                else:
                    no_rank.append(rank_name)
                    record['no_rank'] = ';'.join(no_rank)

            out_file.write(table_writer(record, sep))

    if out_format == 'newick':
        common_tree = newick_tree_generator(tree)
    elif out_format == 'phyloxml':
        common_tree = phyloxml_tree_generator(tree)

    if out_format in {'newick', 'phyloxml'}:
        Phylo.write(common_tree, out_file, out_format)

    return


"############################# Command line args ############################"


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=license,
        epilog=(
            "Example usage:\n"
            "$ %(prog)s -i my_tax.json -f tsv -o my_tax.tsv"
            )
        )
    arg_parser.add_argument(
        "-i", "--infile",
        dest='in_file',
        default=sys.stdin,
        type=argparse.FileType('r'),
        help=(
            "Path to input file containing gi's. "
            "Enter '-' for stdin (default)."
            )
        )
    arg_parser.add_argument(
        "-f", "--format",
        dest='out_format',
        default="tsv",
        choices=['tsv', 'csv', 'phyloxml', 'newick'],
        help=(
            ),
        )
    arg_parser.add_argument(
        "-o", "--outfile",
        dest="out_file",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help=(
            "Path to write output to. "
            "Default is stdout."
            )
        )
    arg_parser.add_argument(
        "--version",
        action="version",
        version='%(prog)s {}'.format(version),
        )

    args = arg_parser.parse_args()
    main(**args.__dict__)

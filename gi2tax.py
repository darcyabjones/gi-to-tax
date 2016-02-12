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

program = "gi2tax"
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

GI_REGEX = r"gi[\|\_]?[\s]*(\d+)\S*"


"################################# Classes ##################################"


class Coder(object):

    def __init__(
            self,
            length=6,
            alphabet=(
                "abcdefghijklmnopqrstuvwxyz"
                "123456789"
                "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                ),
            ):
        self.length = length
        self.alphabet = alphabet
        self.n = len(self.alphabet)
        return

    def __getitem__(self, key):
        return self.index(key)

    def index(self, pattern):
        if len(pattern) == 0:
            return 0

        if isinstance(pattern, str):
            pattern = list(pattern)

        character = pattern.pop(-1)
        return self.n * self.index(pattern) + self.alphabet.index(character)

    def pattern(self, index, k=None):
        if k is None:
            k = self.length

        assert k > 0
        if k == 1:
            return self.alphabet[index]

        prefix_index = index // self.n
        r = index % self.n
        return self.pattern(prefix_index, k - 1) + self.alphabet[r]


"################################# Functions ################################"


def printer(string, clear=None):
    """ Prints running information to screen.

    If clear is specified (int) it will clear that line before printing.
    """
    if isinstance(clear, int):
        sys.stdout.write('\033[{}A\033[K'.format(clear))
    print(string)
    if isinstance(clear, int):
        sys.stdout.write('\033[{}B'.format(clear))
    return


def nodes_table(connection, cursor, handle):
    """ . """

    def table_generator(handle):
        """ returns an iterable generator to write nodes with. """
        for line in handle:
            sline = line.rstrip('\t|\n').split('\t|\t')
            yield tuple(sline[:5])
        return

    # Delete table if it already in the database.
    cursor.execute("drop table if exists nodes")
    cursor.execute((
        "create table nodes\n"
        "(tax_id integer primary key, parent_tax_id integer, "
        "rank text, embl_code text, division_id text)\n"
    ))
    cursor.executemany(
        ("insert into nodes(tax_id, parent_tax_id, rank, embl_code, "
         "division_id) values (?, ?, ?, ?, ?)"),
        table_generator(handle)
    )
    connection.commit()
    return


def names_table(connection, cursor, handle):
    """ . """

    def table_generator(handle):
        """ returns an iterable generator to write nodes with. """
        for line in handle:
            sline = line.rstrip('\t|\n').split('\t|\t')
            tax_id = sline[0]
            name_txt = sline[1]
            unique_name = sline[2]
            name_class = sline[3]
            if name_class == 'scientific name':
                yield (tax_id, name_txt, unique_name, name_class)
        return

    # Delete table if it already in the database.
    cursor.execute("drop table if exists names")
    cursor.execute((
        "create table names\n"
        "(tax_id integer primary key, name_txt text, "
        "unique_name text, name_class text)\n"
        ))

    cursor.executemany((
            "insert into names(tax_id, name_txt, unique_name, name_class) "
            "values (?, ?, ?, ?)"
        ),
        table_generator(handle)
        )
    connection.commit()
    return


def get_gi_taxid_index(handle, index_freq=100000):
    """
    Loops through all lines in the gi_ti file.
    If the lines gi is in the list that we want to know about it passes the
    ti to findTaxPath function. If all of the gi's have been encountered it
    will stop the loop to avoid wasting time.
    """

    index_counter = 0

    # Hold position of the last EOL
    last_position = 0

    handle.seek(0)

    # Note that we can't do 'for line in handle' because we need tell() method.
    for i, line in enumerate(iter(handle.readline, '')):
        gi = 0
        sline = line.rstrip('\n').split('\t')
        gi = int(sline[0])

        if i >= index_counter:
            # Stores end of line for gi.
            yield gi, last_position
            index_counter += index_freq

        last_position = handle.tell()


def get_gi_taxid(handle, gis, gi_taxid_index):
    """
    Loops through all lines in the gi_ti file.
    If the lines gi is in the list that we want to know about it add it to a
    dict. If all of the gi's have been encountered it will stop the loop to
    avoid wasting time.
    """

    len_gi_taxid_index = len(gi_taxid_index)
    last_position = 0  # position of the last EOL

    for i in range(0, len_gi_taxid_index):
        if len(gis) <= 0:
            break
        lower_bound_gi = gi_taxid_index[i][0]
        if i + 1 >= len_gi_taxid_index:
            upper_bound_gi = float('inf')
            upper_bound = float('inf')
        else:
            upper_bound_gi = gi_taxid_index[i + 1][0]
            upper_bound = gi_taxid_index[i + 1][1]
        gi_block = {gi for gi in gis if lower_bound_gi <= gi < upper_bound_gi}
        if len(gi_block) == 0:
            continue
        handle.seek(gi_taxid_index[i][1])
        last_position = gi_taxid_index[i][1]

        for line in iter(handle.readline, ''):
            if len(gi_block) > 0 and handle.tell() <= upper_bound:
                sline = line.rstrip('\n').split('\t')
                gi = int(sline[0])
                ti = int(sline[1])

                if gi in gis:
                    gi_block.discard(gi)
                    yield gi, ti
            else:
                break
    return


def get_gis_ncbixml(handle):
    """ . """
    blast_record = NCBIXML.parse(handle)
    gis = dict()
    for blast_query in blast_record:
        gi = gi_extract(blast_query.query)
        if len(gi) != 0:
            gis[gi] = {
                'gi': gi,
                'id': blast_query.query,
                'description': ''
                }
        for alignment in blast_query.alignments:
            gi = gi_extract(str(alignment.hit_id))
            if len(gi) != 0:
                gis[gi] = {
                    'gi': gi,
                    'id': alignment.hit_id,
                    'description': alignment.hit_def
                    }
            else:
                gi = gi_extract(str(alignment.hit_def))
                if len(gi) != 0:
                    gis[gi] = {
                        'gi': gi,
                        'id': alignment.hit_id,
                        'description': alignment.hit_def
                        }
    return gis


def gi_extract(gi_string):
    """ Takes string object, returns gid. """
    gis = dict()
    regex = re.compile(GI_REGEX, re.I)
    hits = regex.finditer(gi_string.strip())
    for match in hits:
        id_ = match.group(0)
        gi = int(match.group(1))
        gis[gi] = {'gi': gi, 'id': id_}
    return gis


def find_tax_path(cursor, ti, gi_tax_path):
    """
    Recursive function inputs ti, returns dict of all parent tis.
    """

    cursor.execute("""
        SELECT * FROM nodes WHERE tax_id = ?
        """, (ti,))
    """ We specified row_factory as sqlite3.Row earlier so we can access
    values by column name, ACE! """
    node_row = cursor.fetchone()
    rank = node_row['rank']  # eg Superkingdom, kingdom, genus etc.
    parent_taxid = node_row['parent_tax_id']
    cursor.execute("""
        SELECT * FROM names WHERE tax_id = ?
        """, (ti,))
    name_row = cursor.fetchone()
    """ If there is an unique name use that name rather than the duplicate
    Otherwise we can assume that the name_txt column is unique. """
    if name_row['unique_name'] != "":
        rank_name = name_row['unique_name']
    else:
        rank_name = name_row['name_txt']

    new_node_entry = {
        'taxid': ti,
        'parent_taxid': parent_taxid,
        'rank': rank,
        'rank_name': rank_name
        }
    gi_tax_path.append(new_node_entry)

    """ 1 is the lowest node possible.
    All taxonomic paths end up at 1 after superkingdom. """
    if parent_taxid != 1:
        gi_tax_path = find_tax_path(cursor, parent_taxid, gi_tax_path)

    return gi_tax_path


def jsonl_handler(handle):
    try:
        for i, line in enumerate(handle):
            yield i + 1, json.loads(line)
    except ValueError:
        handle.seek(0)
        return enumerate(json.load(handle))


def xsv_handler(handle, sep=r"\t"):
    columns = None
    for i, line in enumerate(handle):
        sline = line.rstrip().split(sep)
        if columns is None:
            columns = sline
            continue

        record = dict()
        for col, val in zip(columns, sline):
            try:
                record[col] = int(val)
            except ValueError:
                record[col] = val
        yield i + 1, record


"################################### Main ###################################"


def main(
        in_file,
        in_format,
        out_file,
        gi_regex=GI_REGEX,
        db_path='tax_db',
        db_type='both',
        self_file=None,
        email=None,
        quiet=False
        ):

    GI_REGEX = gi_regex

    # Set quiet to true if writing output to stdout.
    if out_file == sys.stdout:
        quiet = True

    if db_type == 'both':
        db_type = ['protein', 'nucleotide']
    else:
        db_type = [db_type]

    # Redefine printer if quiet is true.
    if quiet:
        def printer(string, clear=None):
            """ Dummy printer. """
            return
    else:
        def printer(string, clear=None):
            """ Prints running information to screen.

            If clear is specified (int) it will clear that line
            before printing.
            """
            if isinstance(clear, int):
                sys.stdout.write('\033[{}A\033[K'.format(clear))
            print(string)
            if isinstance(clear, int):
                sys.stdout.write('\033[{}B'.format(clear))
            return

    printer("gi2tax")
    printer("======\n")

    printer('Using parameters:')
    printer('- input file: {}'.format(in_file.name))
    printer('- input format: {}'.format(in_format))
    printer('- taxid input path: {}'.format(
        self_file.name if self_file is not None else "none"
        ))
    printer('- output file: {}'.format(out_file.name))
    printer('- database path: {}'.format(db_path))
    printer('- database contains: {} gi\'s\n'.format(' and '.join(db_type)))

    gis = dict()

    printer("Searching for gis using {} method...".format(in_format))

    if in_format == 'regex':
        for line in in_file:
            gis.update(gi_extract(line))
    elif in_format == 'blastxml':
        gis.update(get_gis_from_xml(in_file))
    elif in_format in {'fasta', 'clustal', 'genbank'}:
        sequences = SeqIO.parse(in_file, in_format)
        for seq in sequences:
            these_gis = list(gi_extract(seq.id))
            if len(these_gis) != 0:
                gi = these_gis[0]
                gi['description'] = seq.description
                gis[gi['gi']] = gi
    elif in_format in {'json', 'tsv'}:
        if in_format == 'json':
            handler = jsonl_handler(in_file)
        elif in_format == 'tsv':
            handler = xsv_handler(in_file)

        for i, record in handler:
            if 'gi' in record:
                gis[str(record['gi'])] = record
            elif 'id' in record and 'taxid' in record:
                gis["id:" + str(record['id'])] = record
            else:
                sys.stderr.write(
                    "Error reading line {} in {}. ".format(i, in_file.name) +
                    "Each record must have at least a gi, or an id and taxid."
                    )
                sys.exit()

    coder = Coder()
    for i, (gi, record) in enumerate(gis.items()):
        record['code'] = coder.pattern(i)
        if "name" not in record:
            record['name'] = ""
        if "description" not in record:
            record['description'] = ""

    len_gis = len(gis)

    printer("Searching for gis using {} method... Done".format(in_format), 1)
    printer("Found {} gis\n".format(len_gis, in_format))

    if os.path.isfile(os.path.join(db_path, "ncbi_tax.db")):
        # TODO: check that tables exist in db.
        printer('Found existing SQLite database\n')
    else:
        try:
            printer('Creating SQLite database... ')
            connection = sqlite3.connect(os.path.join(db_path, 'ncbi_tax.db'))
            nodes_handle = open(os.path.join(db_path, 'nodes.dmp'), 'r')
            names_handle = open(os.path.join(db_path, 'names.dmp'), 'r')
            cursor = connection.cursor()
            nodes_table(connection, cursor, nodes_handle)
            names_table(connection, cursor, names_handle)
            printer('Creating SQLite database... Done\n', 1)
        except KeyboardInterrupt:
            sys.stderr.write((
                'Error - user cancelled sqlite db write '
                'in gi2tax, cannot continue.\n'
                ))
            sys.exit()
        finally:
            nodes_handle.close()
            names_handle.close()
            connection.close()

    gi_taxid_dbs = {
        'protein': os.path.join(db_path, 'gi_taxid_prot.dmp'),
        'nucleotide': os.path.join(db_path, 'gi_taxid_nucl.dmp'),
        }

    no_tax_info = list()

    # Need to get GIs
    for db in db_type:
        try:
            handle = open(gi_taxid_dbs[db], 'r')
            connection = sqlite3.connect(os.path.join(db_path, 'ncbi_tax.db'))
            connection.row_factory = sqlite3.Row
            cursor = connection.cursor()

            if not os.path.isfile(gi_taxid_dbs[db] + '.index'):
                handle.seek(0, 2)
                handle_size = handle.tell()
                handle.seek(0)

                printer("Indexing {} gi_taxid file... ".format(db))
                index = list()
                for gi, position in get_gi_taxid_index(handle):
                    index.append((gi, position))
                    printer(
                        "Indexing {} gi_taxid file... {:>4.0%}".format(
                            db,
                            position / handle_size
                            ),
                        1
                        )
                with open(gi_taxid_dbs[db] + '.index', 'w') as json_handle:
                    json.dump(gi_taxid_index, json_handle)
                printer("Indexing {} gi_taxid file... Done".format(db), 1)

            printer("Searching {} gi_taxid file...".format(db))
            with open(gi_taxid_dbs[db] + '.index', 'r') as json_handle:
                gi_taxid_index = json.load(json_handle)

            remaining = len([gi for gi, r in gis.items() if 'taxid' not in r])
            found_count = len_gis - remaining
            for gi, ti in get_gi_taxid(handle, gis, gi_taxid_index):
                gis[gi]['taxid'] = ti
                found_count += 1
                printer("Searching {} gi_taxid file... {} of {}".format(
                    db,
                    found_count,
                    len_gis),
                    1
                    )
                if ti == 0:
                    no_tax_info.append(str(gi))
                else:
                    gis[gi]['tax_path'] = find_tax_path(
                        cursor,
                        ti,
                        list()
                        )
                out_file.write(json.dumps(gis[gi]) + "\n")

            printer("Searching {} gi_taxid file... Done".format(db), 1)

        finally:
            connection.close()
            handle.close()

        remaining = len([gi for gi, r in gis.items() if 'taxid' not in r])
        found_count = len_gis - remaining
        if email is None or remaining <= 0:
            continue

        printer("Searching Entrez {} database...".format(db))
        Entrez.email = email
        outdated = list()

        try:
            connection = sqlite3.connect(os.path.join(db_path, 'ncbi_tax.db'))
            connection.row_factory = sqlite3.Row
            cursor = connection.cursor()

            for gi, record in gis.items():
                if 'taxid' in record:
                    continue
                try:
                    # Raises error if the gi wasn't found in the db
                    entrez_record = Entrez.read(
                        Entrez.esummary(
                            db=db,
                            id=str(gi)
                            )
                        )
                except:
                    raise
                    # taxid not found

                if 'TaxId' in entrez_record[0]:
                    if 'Status' in entrez_record[0]:
                        if (entrez_record[0]['Status'].lower() != 'current' or
                                entrez_record[0]['Status'].lower() != 'live'):
                            outdated.append(gi)
                    record['taxid'] = entrez_record[0]['TaxId']
                    record['tax_path'] = find_tax_path(
                        cursor,
                        ti,
                        list()
                        )
                    out_file.write(json.dumps(record) + "\n")
                    found_count += 1
                    printer(
                        "Searching Entrez {} database... {} of {}".format(
                            db,
                            found_count,
                            len_gis),
                        1
                        )

        finally:
            connection.close()

        printer("Searching for gis using Entrez... Done", 2)

    not_found = [str(k) for k, v in gis.items() if 'taxid' not in v]
    if len(not_found) > 0:
        printer("Note - gi2tax could not find taxids for gis:")
        for i, j in zip(range(0, len(not_found)), range(3, len(not_found))):
            printer("\t" + ", ".join(not_found[i:j]))

    if len(no_tax_info) > 0:
        printer(
            "Note - the organism is unknown for the following gis "
            "(taxid == 0):"
            )
        for i, j in zip(
                range(0, len(no_tax_info)),
                range(3, len(no_tax_info))
                ):
            printer("\t" + ", ".join(no_tax_info[i:j]))


"############################# Command line args ############################"


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=license,
        epilog=(
            "Example usage:\n"
            "$ %(prog)s -i my_fasta.fna -f regex -o my_fasta.faa "
            "-d ./tax_dbs -t protein -e your@email.com\n"
            ""
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
        dest='in_format',
        default="regex",
        choices=['regex', 'blastxml', 'fasta', 'clustal', 'json', 'tsv'],
        help=(
            "The format of the files that you are finding gi's from. "
            "By default gi2tax uses a regular expression pattern."
            ),
        )
    arg_parser.add_argument(
        "-r", "--regex",
        dest='gi_regex',
        default=GI_REGEX,
        type=str,
        help=(
            "A custom regular expression to capture gis. "
            "By default, gi2tax uses the expression {} ".format(GI_REGEX) +
            "for all formats. "
            "See documentation for more information."
            ),
        )
    arg_parser.add_argument(
        "-o", "--outfile",
        dest="out_file",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help=(
            "Path to write output to. "
            "Enter '-' for stdout (default)."
            )
        )
    arg_parser.add_argument(
        "-d", "--dbpath",
        dest="db_path",
        default="",
        help=(
            "Path to the NCBI taxonomy files (a directory). "
            "Default is current working directory."
            )
        )
    arg_parser.add_argument(
        "-t", "--dbtype",
        dest="db_type",
        default='both',
        help=(
            "Specifies which gi to taxid databases (local and Entrez) gi2tax "
            "should search for your gis in. Searching through both will "
            "take longer than going through one. "
            "Default is 'both'."
            ),
        choices=['protein', 'nucleotide', 'both']
        )
    arg_parser.add_argument(
        "-e", "--email",
        default=None,
        help=(
            "Stores an email address to search for gis not found in the "
            "gi_taxid files with ncbi entrez query system."
            )
        )
    arg_parser.add_argument(
        "-q", "--quiet",
        default=False,
        action='store_true',
        help=(
            "Boolean toggle to suppress running feedback. "
            "Always True when output is to stdout."
            )
        )
    arg_parser.add_argument(
        "--version",
        action="version",
        version='%(prog)s {}'.format(version),
        )

    args = arg_parser.parse_args()
    main(**args.__dict__)

#!/usr/local/bin/env python

#Version 1. Darcy Jones, January 2014.
#Contact Darcy Jones, darcy.ab.jones@gmail.com

    #----------------------------------- LICENSE ---------------------------------------#
    #    gi2tax - finds taxonomic information for a given set of genbank identifiers    #
    #    Copyright (C) 2013  Darcy Jones                                                #
    #                                                                                   #
    #    This program is free software: you can redistribute it and/or modify           #
    #    it under the terms of the GNU General Public License as published by           #
    #    the Free Software Foundation, either version 3 of the License, or              #
    #    (at your option) any later version.                                            #
    #                                                                                   #
    #    This program is distributed in the hope that it will be useful,                #
    #    but WITHOUT ANY WARRANTY; without even the implied warranty of                 #
    #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  #
    #    GNU General Public License for more details.                                   #
    #                                                                                   #
    #    You should have received a copy of the GNU General Public License              #
    #    along with this program.  If not, see <http://www.gnu.org/licenses/>.          #
    #                                                                                   #
    #-----------------------------------------------------------------------------------#

###Import modules

import sys;
import argparse;
from ftplib import FTP;
import os;
import time;
import pickle;
import json;
import calendar;
import sqlite3;
import zipfile;
import re;
from Bio import SeqIO, Phylo;
from Bio.Blast import NCBIXML; #for generic xml handling
from Bio.Phylo.PhyloXML import Phylogeny, Clade, Taxonomy, Id;
from Bio.Phylo.BaseTree import Tree;
from Bio.Phylo.BaseTree import Clade as Clade2



###Function definitions.


def debug(string, gap_=False):
    pass;

def findDbPathContents():
    # Find what files are available to us.
    # Returns a dictionary keyed by file name, with values = a list of relative paths to that file.
    debug('### funct findDbPathContents ###', gap_=True);

    db_path_files_= dict(zip(os.listdir(db_path), [[os.path.join(db_path, file_)] for file_ in os.listdir(db_path)])); # Gets list of files already in the db directory.

    if os.path.isdir(os.path.join(db_path, 'taxdmp')):
        for file_ in os.listdir(os.path.join(db_path, 'taxdmp')):
            if file_ in db_path_files_:
                db_path_files_[file_].append(os.path.join(db_path, 'taxdmp', file_));
            else:
                db_path_files_[file_] = [os.path.join(db_path, 'taxdmp', file_)];
    for file_ in db_path_files_:
        db_path_files_[file_].sort(key=lambda file_: os.stat(file_).st_mtime); # Sort in place by mtime, small to big

    debug("All files in db_path: {}".format(db_path_files_.values()));
    return db_path_files_;

def updateArchives():
    if not quiet:
        print('##### Checking local tax archives  #####');
        print('Connecting to ncbi ftp.');
    debug('### funct updateArchives ###');
    ftp_=connectTaxFTP(); # Connect to the ncbi ftp server
    if not quiet:
        print('Finding which archive files to download.');
    returned_=checkTaxFTPData(ftp_); # Determine which files you need to download or update.
    to_get_=returned_[0];
    to_keep_=returned_[1];
    success_=True;
    if len(to_get_)>0:
        if not quiet:
            print('Found {} archives to be updated.'.format(len(to_get_)));
            print('\n###### Downloading files from ftp ######');
        got_=getTaxFTPData(ftp_, to_get_); # Download files that checkTaxFTPData said that you need.
        if len(to_get_) != len(got_): # If we didn't get all of the files we didnt succeed.
            print('could not successfully download all files.');
            print('Could not retrieve {}.'.format(set(to_get_)-set(got_)));
            success_=False;
        elif not quiet:
            print('successfully downloaded archive files.');
    elif not quiet:
        print('All archive files are current and complete.');
    debug('updateArchives successful?: {}'.format(success_))
    return to_keep_;

def connectTaxFTP():
    # Connects to the taxonomic ftp server and moves to the right directory.
    # Check is here because ftplib return values. Later we might check for success/failure before proceeding.
    ftp_=FTP(tax_ftp); # Open new connection
    check=ftp_.login(); # Login as anonymous user
    debug(check);
    check=ftp_.cwd(tax_ftp_path); # change current working directory to taxonomy directory.
    debug(check);
    check=ftp_.sendcmd('type i'); # Change connection to binary mode (so that i can read filesizes easily and download in blocks.)
    debug(check);
    return ftp_; # return ftp handle.

def checkTaxFTPData(ftp_):
    # Checks whether tax files need to be downloaded/updated.
    # Note that because of possible time differences between the server and the local computer, this method is not always accurate in comparing dates. Later versions might store modification dates in a 'versions' file.
    debug('### funct checkTaxFTPData ###', gap_=True);
    # Find what files are available to us.
    if os.path.isfile(os.path.join(db_path,'VERSIONS')):
        with open(os.path.join(db_path,'VERSIONS'), 'r') as versions_:
            versions_dict_=pickle.load(versions_);
    else:
        versions_dict_=dict();
 
    debug('### Versions file ###');
    for d in versions_dict_:
        debug('\tVersions file: {}={}.'.format(d, versions_dict_[d]));

    to_get_=set(); # Set of files to download
    to_keep_=dict(); # Dictionary keyed by preferred archive, with value = a list of paths to useable files. 

    for req_file_ in req_files: # for group in req_files, eg loops through gi_taxid_prot, taxdump ... in dictionary
        debug('Checking versions of: {}'.format(req_file_), gap_=True);

        # Find what which set of files from the required group are present.
        present_=[];

        if req_file_ in db_path_files:
            present_.append(req_file_);
        elif req_files[req_file_] in db_path_files:
            present_.append(req_files[req_file_]);

        debug("hits for {}: {}".format(req_file_, ' '.join(present_)));

        if len(present_) == 0: # if there are no archives or files for the data.
            to_get_.add(req_files[req_file_]); # Download the preferred data archive.
            debug('No hits for req file(s) {}, will download containing archive: {}.'.format(req_file_, req_files[req_file_]));
        else:
            present_paths_=[];
            for file_ in present_:
                present_paths_.extend(db_path_files[file_]);
            present_paths_.sort(key=lambda file_: os.stat(file_).st_mtime, reverse=True); # Sort in place by mtime, big to small
            debug('File hits ordered by mtime: {}'.format(present_paths_));

            returned_=checkTaxFTPCompare(ftp_, present_paths_, versions_dict_);
            if returned_ == None:
                to_get_.add(req_files[req_file_]);
            else:
                to_keep_[req_file_] = returned_;

    return (to_get_, to_keep_);

def checkTaxFTPCompare(ftp_, present_paths_, versions_dict_):
    # Loop through all hits untill we get a hit that is up to date and the (right size, to come later) so that get_ ==false. if get_ is still true after we go through all files the preferred archive is added to the download list. If an appropriate file is found, the path to that file is added to a dictionary.
    debug('### funct checkTaxFTPCompare ###', gap_=True);

    preferred_archive_dict={
    'gi_taxid_prot.dmp':'gi_taxid_prot.zip',
    'gi_taxid_prot.zip':'gi_taxid_prot.zip',
    'gi_taxid_prot.dmp.gz':'gi_taxid_prot.zip',
    'gi_taxid_nucl.dmp':'gi_taxid_nucl.zip',
    'gi_taxid_nucl.zip':'gi_taxid_nucl.zip',
    'gi_taxid_nucl.dmp.gz':'gi_taxid_nucl.zip',
    'categories.dmp':'taxcat.zip',
    'taxcat.zip':'taxcat.zip',
    'taxcat.tar.gz':'taxcat.zip',
    'names.dmp':'taxdmp.zip',
    'nodes.dmp':'taxdmp.zip',
    'taxdmp.zip':'taxdmp.zip',
    'taxdump.tar.gz':'taxdmp.zip'}

    containing_archive_=preferred_archive_dict[os.path.split(present_paths_[0])[1]];
    ftp_file_mtime_= calendar.timegm(time.strptime(ftp_.sendcmd('MDTM {}'.format(containing_archive_))[4:],'%Y%m%d%H%M%S')); # http://alexharvey.eu/code/python/get-a-files-last-modified-datetime-using-python/, see also http://docs.python.org/3/library/time.html                      
    ok_path_list_=[];
    delete_path_set_=set();

    for path_ in present_paths_:
        local_file_ = path_;
        local_file_mtime_ = os.stat(local_file_).st_mtime;

        debug('Checking version of {}'.format(os.path.split(local_file_)[1]));
        debug('Mod times for {}: local {}, ftp {}'.format(local_file_, local_file_mtime_, ftp_file_mtime_));

        if ftp_file_mtime_ > (local_file_mtime_+3600): # if the ftp version is newer than the local version. NB 3600 == 1 hour, it is added to account for any daylight savings difference.
            file_true_=False; 
            debug('test = {} > {} + 3600 = True'.format(ftp_file_mtime_, local_file_mtime_));
            delete_path_set_.add(local_file_);

        else:
            debug('test = {} > {} + 3600 = False'.format(ftp_file_mtime_, local_file_mtime_));
            file_true_=localArchiveQC(local_file_, ftp_, versions_dict_); # Checks if the archive is ok. If not then it is added to the list to update. Returns true false,
            if file_true_:
                ok_path_list_.append(local_file_);
                get_=False;
            else:
                delete_path_set_.add(local_file_);

    if len(ok_path_list_) >1:
        if '.dmp' in ok_path_list_[1] and '.dmp' not in ok_path_list_[0] and os.stat(ok_path_list_[1]).st_mtime >= os.stat(ok_path_list_[1]).st_mtime:
            best_hit_=ok_path_list_[1];
        else:
            best_hit_=ok_path_list_[0];
        debug('Found complete set of required files {}'.format(to_keep_[containing_archive_]));
    elif len(ok_path_list_)==1:
        best_hit_=ok_path_list_[0];
    else:
        debug('No up to date sets found in {}, will download preferred archive {}'.format(present_paths_, containing_archive_));
        best_hit_=None;
    for file_ in delete_path_set_: # If the files are outdated or corrupted we delete them to keep it clean
        try:
            os.remove(file_);
            debug('Deleted old file {}'.format(file_[0]));
        except:
            debug('Failed to delete old file {}'.format(file_[0]));
            debug(sys.exc_info());

    return best_hit_

def getTaxFTPData(ftp_, to_get_):
    # Retrieves files from tax_ftp that are new or are not in the directory.
    #loop through to_get, download each file. add to set got on success.

    debug('### funct getTaxFTPData ###', gap_=True);
    got_=set();
    if os.path.isfile(os.path.join(db_path,'VERSIONS')):
        with open(os.path.join(db_path,'VERSIONS'), 'r') as versions_:
            versions_dict_=pickle.load(versions_);
    else:
        versions_dict_=dict();

    debug('### Versions file ###');
    for d in versions_dict_:
        debug('\tVersions file: {}={}.'.format(d, versions_dict_[d]));

    for file_ in to_get_:
        try:
            debug('Attempting to download: {}'.format(file_), gap_=True);
            ftp_file_size_=int(ftp_.size(file_));
            debug('The downloaded file should be {} bytes'.format(ftp_file_size_));
            ftp_file_mtime_= calendar.timegm(time.strptime(ftp_.sendcmd('MDTM {}'.format(file_))[4:],'%Y%m%d%H%M%S')); # http://alexharvey.eu/code/python/get-a-files-last-modified-datetime-using-python/, see also http://docs.python.org/3/library/time.html
            debug('The downloaded file was last modified {} seconds since epoch'.format(ftp_file_mtime_));
            if not quiet:
                print('Retrieving file: {}, size: {} bytes'.format(file_, str(ftp_file_size_)));
            tries=0;
            success=False;
            while success==False and tries < 10:
                try:
                    debug('Download attempt {}'.format(tries));
                    file_path_ = os.path.join(db_path, file_)
                    check=ftp_.retrbinary('RETR {}'.format(file_), open(file_path_, 'w').write);
                    debug(check);
                    if os.stat(file_path_).st_size == ftp_file_size_: # If we have successfully downloaded the file, break the loop and add the information to the versions pickle dict.
                        success=True;
                        versions_dict_[file_.split('.')[0]]={'path':file_path_, 'mtime':ftp_file_mtime_, 'size':ftp_file_size_};
                        got_.add(file_path_); # Add the file to the list of files successfully downloaded.

                    else:
                        tries+=1;
                except KeyboardInterrupt:
                    sys.stderr.write('Keyboard interrupt: user cancelled download in gi2tax\n');
                    sys.exit();
                except:
                    tries+=1;
                    debug(sys.exc_info());
        except KeyboardInterrupt:
            sys.stderr.write('Keyboard interrupt: user cancelled download in gi2tax\n');
            sys.exit();
        except:
            debug(sys.exc_info());

    with open(os.path.join(db_path,'VERSIONS'), 'w') as versions_:
        pickle.dump(versions_dict_, versions_, protocol=2); # Write the updated versions dict to the Versions file. NB i'm using protocol 2 for backwards compatibilty with python 2.*

    return got_; # Return a set of files successfully retrieved.

def localArchiveQC(file_, ftp_, versions_dict_): 
    # Checks that the local copy of the archive is the size that it should be. If it's not we assume that it is corrupted.
    debug('### funct localArchiveQC ###', gap_=True);
    file_name_=os.path.split(file_)[1].split('.')[0];

    local_file_size_=os.stat(file_).st_size; # find the size of the local archive file.
    if file_name_ in versions_dict_: # if there is size information in the versions pickle dict use that
        ftp_file_size_=versions_dict_[file_name_]['size'];
    elif file_.split('.')[-1] != 'dmp': # otherwise connect to the ftp server and retrieve size info from there if the file is an archive.
        try:
            ftp_file_size_=ftp_.size(file_);
        except:
            debug(sys.exc_info());
            ftp_file_size_=0;
    if local_file_size_ != ftp_file_size_: # if the file sizes are different we assume that it is corrupted or incomplete and download it again.
        debug('test: {} != {} = True'.format(local_file_size_, ftp_file_size_));
        return False;
    else:
        debug('test: {} != {} = False'.format(local_file_size_, ftp_file_size_));
        return True; # ie it is ok


def archiveFileExtractor(archive_path_, member_): # Archive path is the archive that you want to use. member_ is the file within the archive that you want.
    #   Takes an archive file in .zip, .tar.gz, .gz, .tar.Z and returns an open file, saved in a temporary directory
    #   Note that member really should always be specified.
    debug('### funct archiveFileExtractor ###', gap_=True);

    if not quiet:
        print('\nExtracting {} from archive {}...'.format(member_, archive_path_));
    if zipfile.is_zipfile(archive_path_):
        archive_=zipfile.ZipFile(archive_path_);
        extracted_file_ = archive_.extract(member_, db_path); # Note that this returns a file into the temp folder.
        if not quiet:
            sys.stdout.write("\033[1A\033[K"); # Moves one line up, and clears to the end of the line.
            print('Extracting {} from archive {}...  Done'.format(member_, archive_path_));
        return open(extracted_file_, 'rU');

    else:
        try:
            if os.path.isfile(archive_path_):
                if not quiet:
                    sys.stdout.write("\033[1A\033[K"); # Moves one line up, and clears to the end of the line.
                    print('Extracting {} from archive {}...  Done'.format(member_, archive_path_));
                return open(archive_path_);
        except:
            debug(sys.exc_info());
            print('The required archive, {} containing {} could not be found.'.format(archive_path_, member_));
            print('Please run the script again in update mode, or extract the files yourself to the db_path');
            sys.exit();

def updateDB(archive_files_):
    if not quiet:
        print('\n#####  Generating SQLite database  #####');
    debug('### funct updateDB ###');
    connection_= sqlite3.connect(os.path.join(db_path, 'ncbi_tax.db'));
    cursor_= connection_.cursor();

    try:
        if not quiet:
            print('Generating table: nodes.');
        extracted_file_ = archiveFileExtractor(archive_files_['nodes.dmp'], 'nodes.dmp');
        success_= nodesTable(connection_, cursor_, extracted_file_);
        extracted_file_.close()
        if not quiet:
            print('Generating table: names.');
        extracted_file_ = archiveFileExtractor(archive_files_['names.dmp'], 'names.dmp');
        success_= namesTable(connection_, cursor_, extracted_file_);
        extracted_file_.close();

    except KeyboardInterrupt:
        sys.stderr.write('Keyboard Interrupt: user cancelled sqlite db write in gi2tax, cannot continue.\n');
        sys.exit();
    except:
        debug(sys.exc_info());
        success_=False;
    finally:
        try: # We dont want these extracted files hangin around. Delete them and save some space!
            if isinstance(extracted_file_, file):
                extracted_file_.close();
        except:
            debug(sys.exc_info()); # nb errno # 2 is no such directory
    return success_;

def nodesTable(connection_, cursor_, tempfile_):
    debug('### funct nodesTable ###', gap_=True);
    success_=True;
    cursor_.execute("drop table if exists nodes"); #  if the table is already in the database, ie if we are updating the database, delete the previous table.
    cursor_.execute('''
        create table nodes
        (tax_id integer primary key, parent_tax_id integer, rank text, embl_code text, division_id text)
    ''') # Create a new table with the columns, ti (taxoonomic identifier), parent nodes taxonomic id, and the taxonomic rank (Kingdom, phylum, class etc)
    debug('\tGenerated new table: nodes');
    if not quiet:
        print('\tWriting taxonomic nodes to database (Be patient).');

    cursor_.executemany("insert into nodes(tax_id, parent_tax_id, rank, embl_code, division_id) values (?, ?, ?, ?, ?)", nodesTableGenerator(tempfile_));
    connection_.commit();

    if not quiet:
        print('\tSuccessfully created nodes table in SQLite database');

    return success_;

def nodesTableGenerator(file_): # returns an iterable generator to write nodes with.
    for line_ in file_:
        line_list_=line_.rstrip('\t|\n').split('\t|\t');
        tax_id=line_list_[0];
        parent_tax_id=line_list_[1];
        rank=line_list_[2];
        embl_code=line_list_[3];
        division_id=line_list_[4];
        yield (tax_id, parent_tax_id, rank, embl_code, division_id);

def namesTable(connection_, cursor_, tempfile_):
    debug('### funct namesTable ###', gap_=True);
    success_=True;

    cursor_.execute("drop table if exists names"); #  if the table is already in the database, ie if we are updating the database, delete the previous table.
    cursor_.execute('''
        create table names
        (tax_id integer primary key, name_txt text, unique_name text, name_class text)
    ''') # Create a new table with the columns, ti (taxoonomic identifier), parent nodes taxonomic id, and the taxonomic rank (Kingdom, phylum, class etc)
    debug('\tGenerated new table: names');

    if not quiet:
        print('\tAssigning names for each node in database.');

    cursor_.executemany("insert into names(tax_id, name_txt, unique_name, name_class) values (?, ?, ?, ?)", namesTableGenerator(tempfile_));
    connection_.commit();

    if not quiet:
        print('\tSuccessfully created names table in SQLite database');
    return success_;

def namesTableGenerator(file_): # returns an iterable generator to write nodes with.
    for line_ in file_:
        line_list_=line_.rstrip('\t|\n').split('\t|\t');
        tax_id=line_list_[0];
        name_txt=line_list_[1];
        unique_name=line_list_[2];
        name_class=line_list_[3];
        if name_class == 'scientific name':
            yield (tax_id, name_txt, unique_name, name_class);

def findArchiveFiles(db_path_files_, archive_files_):
    for file_ in req_files:
        if file_ in db_path_files_ and file_ not in archive_files_:
            archive_files_[file_] = db_path_files_[file_][0];
        elif req_files[file_] in db_path_files_ and file_ not in archive_files_:
            archive_files_[file_] = db_path_files_[req_files[file_]][0];
    return archive_files_;


def giFinder(input_path_, input_format_):
    # In > input_file containing gi's to find info for. 
    # Out > set of gi's
    gis_=set();
    if input_path_ == None:
        return {}
    try:
        if input_path_ == sys.stdin or input_path_ == '-':
            in_file_ = sys.stdin;
        else:
            in_file_ = open(input_path_, 'rU');

        if input_format_ == 'regex':
            regex=re.compile('gi[\|\_]?[\s]*(\d+)', re.I); # matches gi|[any amount of white space][numbers], gi_[any amount of white space][numbers], gi[any amount of white space][numbers],  gi_[numbers][any amount of white space], gi[any amount of white space][numbers]. (Case insensitive)
            gis_=set(regex.findall(in_file_.read()));
        elif input_format_ == 'blastXML':
            gis_= getGisFromXML(in_file_);
        elif input_format_ == 'plain':
            regex=re.compile("(\d+)"); # Matches any set of numbers
            gis_=set(regex.findall(in_file_.read()));
        elif input_format_ in {'fasta', 'clustal', 'genbank'}:
            sequences_=SeqIO.parse(in_file_, input_format_);
            regex=re.compile('gi[\|\_]?[\s]*(\d+)', re.I); # matches gi|[any amount of white space][numbers], gi_[any amount of white space][numbers], gi[any amount of white space][numbers],  gi_[numbers][any amount of white space], gi[any amount of white space][numbers]. (Case insensitive)
            for seq_ in sequences_:
                gis_|=set(regex.findall(seq_.id));
                # gis_|=set(regex.findall(seq_.name); # Maybe later
                # gis_|=set(regex.findall(seq_.description);
    except:
        debug(sys.exc_info());
        raise;
    finally:
        if in_file_ != sys.stdin:
            in_file_.close();

    if len(gis_) == 0:
        debug("No gi\'s found in {} using {} search method. Try again.".format(input_path_, input_format_));
        gis_={};

    return {int(gi) for gi in gis_};

def getGisFromXML(in_file_):
    blast_record_ = NCBIXML.parse(in_file_);
    gis_=set();
    for blast_query_ in blast_record_:
        gi_=gi_extract(blast_query_.query)
        if len(gi_) != 0:
            gis_|=set(gi_);
        for alignment_ in blast_query_.alignments:    
            for hsp_ in alignment_.hsps:
                gi_ = gi_extract(str(alignment_.title));
                #hit_evalue = float(hsp.expect);
                if len(gi_)!=0:
                    gis_|=set(gi_);
    return gis_;

def gi_extract(gi_string): #takes string object, returns gid
    regex = re.compile("gi[\|\_]?[\s]*(\d+)", re.I);
    hit_id = regex.findall(gi_string.strip());
    #if len(hit_id)!=0: #handles the event that the query has no gid, ie it may be an unpublished sequence.
    return hit_id;


def giTaxidContainer(gi_ti_file_, archive_files_, gis_, index_=False):
    debug('### funct giTaxidContainer ###', gap_=True);
    connection_= sqlite3.connect(os.path.join(db_path, 'ncbi_tax.db'));
    connection_.row_factory = sqlite3.Row;
    cursor_= connection_.cursor();
    try:
        extracted_file_ = archiveFileExtractor(archive_files_[gi_ti_file_], gi_ti_file_);
        if index_:
            gis_ = indexGiTaxidFile(cursor_, extracted_file_, gi_ti_file_, gis_);
        else:
            gis_ = giTaxidFile(cursor_, extracted_file_, gi_ti_file_, gis_);
    except:
        debug('Index Failed');
        debug(sys.exc_info());
        if debug:
            raise;
    finally:
        extracted_file_.close();
        debug("Closed {}".format(gi_ti_file_))
    return gis_;

def indexGiTaxidFile(cursor_, extracted_file_, gi_ti_file_, gis_):
    # Loops through all lines in the gi_ti file, if the lines gi is in the list that we want to know about it passes the ti to findTaxPath function. If all of the gi's have been encountered it will stop the loop to avoid wasting time. Will return a set of gi's that weren't found.
    debug('### funct indexGiTaxidFile ###', gap_=True);
    if not quiet:
        print('\tSearching for gi\'s in {} and finding taxonomic paths.'.format(gi_ti_file_));

    gi_dict_={};
    # gi_dict will have the structure. {gi:{'species':'speciesName', 'species_ti':tax_id, 'genus':'genusName', 'genus_ti':tax_id ... 'superKingdom':'Superkingdom_name' etc}, anotherGI:{}};
    gi_taxid_index_=[];
    len_gi_taxid_index_=0;
    max_len_gi_taxid_index_=200000; # Makes max size of the file approx 2 mb. Memory size probably slightly bigger.

    i=0;
    index_freq=500000;
    index_counter=0;

    last_position=0; #position of the last EOL

    found_count_=0;
    len_gis_=len(gis_);
    if not quiet:
        print('\tFound 0 of {} gi\'s'.format(len_gis_));
        print('\tSearched {} lines in {}'.format(i, gi_ti_file_));
    for line_ in iter(extracted_file_.readline, ''):
        record_= line_.rstrip('\n').split('\t');
        gi_ = int(record_[0]);
        ti_ = int(record_[1]);
        if gi_ in gis_:
            gi_dict_[gi_]={'gi':gi_,'taxid':ti_, 'tax_path':[]};
            gi_dict_[gi_]['tax_path']=findTaxPath(cursor_, ti_, gi_dict_[gi_]['tax_path']);
            output_handler(gi_dict_[gi_]);
            gis_.discard(gi_);
            if not quiet:
                sys.stdout.write("\033[2A\033[K"); # Moves two lines up, and clears to the end of the line.
                print('\tFound {} of {} gi\'s'.format(found_count_, len_gis_));
                sys.stdout.write("\033[2B"); # Moves two lines down
            if len_gi_taxid_index_<max_len_gi_taxid_index_: # Adds index for found gis, impoves researching same gis
                gi_taxid_index_.append((gi_, last_position)); # Stores end of line for gi.
                len_gi_taxid_index_+=1;
            found_count_+=1;
        if i>=index_counter:
            gi_taxid_index_.append((gi_, last_position)); # Stores end of line for gi.
            index_counter+=index_freq;
            len_gi_taxid_index_+=1
            if not quiet:
                sys.stdout.write("\033[1A\033[K"); # Moves one line up, and clears to the end of the line.
                print('\tSearched {} lines in {}'.format(i, gi_ti_file_));
        last_position=extracted_file_.tell()
        i += 1;
    if not quiet:
        sys.stdout.write("\033[1A\033[K"); # Moves one line up, and clears to the end of the line.
        print('\tCompleted writing {} gi and ti data'.format(gi_ti_file_));

    with open(os.path.join(db_path,'{}.index'.format(os.path.splitext(gi_ti_file_)[0])), 'w') as index_:
        pickle.dump(gi_taxid_index_, index_, protocol=2); # Write the updated versions dict to the Versions file. NB i'm using protocol 2 for backwards compatibilty with python 2.*

    return gis_; # Any problems ? This set should have len = 0. if len >0 the program wasn't able to find the gi in the gi_ti file given (perhaps nuc vs protein or bad gi call);

def giTaxidFile(cursor_, extracted_file_, gi_ti_file_, gis_):
    # Loops through all lines in the gi_ti file, if the lines gi is in the list that we want to know about it passes the ti to findTaxPath function. If all of the gi's have been encountered it will stop the loop to avoid wasting time. Will return a set of gi's that weren't found.
    debug('### funct giTaxidFile ###', gap_=True);
    if not quiet:
        print('\tSearching for gi\'s in {} and finding taxonomic paths.'.format(gi_ti_file_));

    with open(os.path.join(db_path,'{}.index'.format(os.path.splitext(gi_ti_file_)[0])), 'r') as index_:
        gi_taxid_index_ = pickle.load(index_); # Load the index. NB i'm using protocol 2 for backwards compatibilty with python 2.*

    new_gi_taxid_index_=[]; # We write this one as we go to add new indices for searched gis
    len_new_gi_taxid_index_=0;

    gi_dict_={};
    # gi_dict will have the structure. {gi:{'species':'speciesName', 'species_ti':tax_id, 'genus':'genusName', 'genus_ti':tax_id ... 'superKingdom':'Superkingdom_name' etc}, anotherGI:{}};

    len_gi_taxid_index_=len(gi_taxid_index_);
    max_len_gi_taxid_index_=200000; # Makes max size of the file approx 2 mb. Memory size probably slightly bigger.
    debug('There are {} blocks in index file'.format(len_gi_taxid_index_));
    last_position=0; #position of the last EOL

    found_count_=0;
    len_gis_=len(gis_);
    if not quiet:
        print('\tFound 0 of {} gi\'s'.format(len_gis_));

    i=0;
    while i<len_gi_taxid_index_:
        if len(gis_) >0:
            lower_bound_gi_=gi_taxid_index_[i][0]
            if i+1 >= len_gi_taxid_index_:
                upper_bound_gi_= float('inf');
                upper_bound_ = float('inf');
            else:
                upper_bound_gi_=gi_taxid_index_[i+1][0];
                upper_bound_= gi_taxid_index_[i+1][1];
            gi_block_={gi_ for gi_ in gis_ if lower_bound_gi_<= gi_ <upper_bound_gi_};
            if len(gi_block_) >0:
                extracted_file_.seek(gi_taxid_index_[i][1]);
                last_position=gi_taxid_index_[i][1];
                for line_ in iter(extracted_file_.readline, ''):
                    if len(gi_block_)>0 and extracted_file_.tell()<=upper_bound_:
                        _record_= line_.rstrip('\n').split('\t');
                        gi_ = int(_record_[0]);
                        if gi_ in gis_:
                            ti_ = int(_record_[1]);
                            gi_dict_[gi_]={'gi':gi_,'taxid':ti_, 'tax_path':[]};
                            gi_dict_[gi_]['tax_path']=findTaxPath(cursor_, ti_, gi_dict_[gi_]['tax_path']);
                            output_handler(gi_dict_[gi_]);
                            gis_.discard(gi_);
                            gi_block_.discard(gi_);
                            if not quiet:
                                sys.stdout.write("\033[1A\033[K"); # Moves two lines up, and clears to the end of the line.
                                print('\tFound {} of {} gi\'s'.format(found_count_, len_gis_));
                            if len_gi_taxid_index_+len_new_gi_taxid_index_ < max_len_gi_taxid_index_ and gi_ != lower_bound_gi_: # Adds index for found gis, if not already in the list. impoves researching same gis Looking for best way to add this.
                                new_gi_taxid_index_.append((gi_, last_position)); # Stores end of line for previous gi. (beginning of line for this one).
                                len_new_gi_taxid_index_+=1;
                            found_count_+=1;
                    else:
                        break;
                    last_position=extracted_file_.tell()
        else:
            break;
        i+=1;

    if not quiet:
        print('\tCompleted writing {} gi and ti data'.format(gi_ti_file_));

    with open(os.path.join(db_path,'{}.index'.format(os.path.splitext(gi_ti_file_)[0])), 'w') as index_:
        pickle.dump(sorted(gi_taxid_index_+new_gi_taxid_index_, key=lambda tup: tup[0]), index_, protocol=2); # Write the updated indexto the index file. NB i'm using protocol 2 for backwards compatibilty with python 2.*

    return gis_; # Any problems ? This set should have len = 0. if len >0 the program wasn't able to find the gi in the gi_ti file given (perhaps nuc vs protein or bad gi call);

def findTaxPath(cursor_, ti_, gi_tax_path_): #Recursive function inputs ti_, returns dict of all parent tis;
    debug("### funct findTaxPath ###", gap_=True);
    cursor_.execute("""
        SELECT * FROM nodes WHERE tax_id = ?
        """, (ti_,));
    node_row_=cursor_.fetchone(); # We specified row_factory as sqlite3.Row earlier so we can access values by column name, ACE!
    rank=node_row_['rank']; # eg Superkingdom, kingdom, genus etc.
    parent_taxid=node_row_['parent_tax_id'];

    cursor_.execute("""
        SELECT * FROM names WHERE tax_id = ?
        """, (ti_,));
    name_row_=cursor_.fetchone();
    
    if name_row_['unique_name'] != "": # If there is an unique name specified use that name rather than the duplicate
        rank_name=name_row_['unique_name'];
    else: # if there wasn't an unique name specified we can assume that the name_txt column is unique.
        rank_name=name_row_['name_txt'];

    debug('{} is a {} called {}'.format(ti_, rank, rank_name));

    new_node_entry_={'taxid':ti_, 'parent_taxid':parent_taxid, 'rank':rank, 'rank_name':rank_name};
    gi_tax_path_.append(new_node_entry_);
    if parent_taxid != 1: # 1 is the lowest node possible. All taxonomic paths end up at 1 after superkingdom.
        gi_tax_path_=findTaxPath(cursor_, parent_taxid, gi_tax_path_);
    
    return gi_tax_path_;


def output_handler(record_):
    debug('funct output_handler');
    if 'newick' in output_format or 'phyloXML' in output_format:
        treeFormatHandler(record_);

    if 'SQLite' in output_format or 'pickle' in output_format or 'json' in output_format or 'tab' in output_format:
        row_=columnFormatHandler(record_);

    if 'SQLite' in output_format:
        output_cursor.execute("""
        insert into taxonomy(gi,taxid,subspecies,species,subgenus,genus,subfamily,family,superfamily,suborder,order_,superorder,subclass,class,superclass,subphylum,phylum,superphylum,subkingdom,kingdom,superkingdom,no_rank)
         values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, rowWriter(row_, 'SQLite'));

    if 'pickle' in output_format or 'json' in output_format:
        output_dict[record_['gi']] = row_;

    if 'tab' in output_format:
        out_tab_file.write(rowWriter(row_, 'tab'));

def treeFormatHandler(record_):
    debug('funct treeFormatHandler', gap_=True);
    for node_ in record_['tax_path']:
        if node_['taxid'] in tree_dict:
            if 'rank' not in tree_dict[node_['taxid']]:
                tree_dict[node_['taxid']]['rank']=node_['rank'];
            if 'rank_name' not in tree_dict[node_['taxid']]:
                tree_dict[node_['taxid']]['rank_name']=node_['rank_name'];
            if 'children' not in tree_dict[node_['taxid']]:
                tree_dict[node_['taxid']]['children'] = set();
        else:
            tree_dict[node_['taxid']]={'rank':node_['rank'], 'rank_name':node_['rank_name'], 'children':set()};
        if node_['parent_taxid'] in tree_dict:
            tree_dict[node_['parent_taxid']]['children'].add(node_['taxid']);
        else:
            tree_dict[node_['parent_taxid']]={'children':{node_['taxid']}};
    if record_['taxid'] in tree_dict:
        tree_dict[record_['taxid']]['gi']=record_['gi'];
    else:
        tree_dict[record_['taxid']]={'gi':record_['gi']};

def columnFormatHandler(record_):
    debug('funct columnFormatHandler')
    row_={'gi':record_['gi'], 'taxid':record_['taxid']};
    no_rank=[];
    for node_ in record_['tax_path']:
        rank=node_['rank'];
        rank_name=node_['rank_name'];
        if rank != 'no rank':
            row_['{}_taxid'.format(rank)]=node_['taxid'];
            row_[rank]=rank_name;
        else:
            no_rank.append(rank_name);
    row_['no_rank']=';'.join(no_rank);
    return row_;

def rowWriter(row_, type_):
    debug('funct rowWriter');
    debug(row_);
    if type_=='tab':
        out_row_=[str(row_['gi']), str(row_['taxid'])];
    else:
        out_row_=[row_['gi'], row_['taxid']];
    column_headings_=['subspecies','species','subgenus','genus','subfamily','family','superfamily','suborder','order','superorder','subclass','class','superclass','subphylum','phylum','superphylum','subkingdom','kingdom','superkingdom','no_rank'];

    for rank_ in column_headings_:
        if rank_ in row_:
            out_row_.append(row_[rank_]);
        else:
            out_row_.append('.');

    if type_ == 'SQLite':
        return tuple(out_row_);
    elif type_ == 'tab':
        return '{}\n'.format('\t'.join(out_row_));


def phyloxmlTreeGenerator():
    debug('funct phyloxmlTreeGenerator', gap_=True);
    children=phyloxmlTreeGeneratorRecurse(taxid=1);
    
    return Phylogeny(name='gi2tax - Common tree', root= children);

def phyloxmlTreeGeneratorRecurse(taxid=1):
    children_=[];
    for child_ in tree_dict[taxid]['children']:
        children_.append(phyloxmlTreeGeneratorRecurse(child_));
    if len(children_)==0:
        children_=None;

    if taxid == 1:
        rank_name_=None;
        node_id_=1, 
        tax_data_=None

    else:
        if 'rank' in tree_dict[taxid]:
            rank_=tree_dict[taxid]['rank']
            if rank_=='no_rank' or rank_=='no rank' or rank_=='superkingdom':
                rank_=None;
        else:
            rank_=None;
        
        if 'gi' in tree_dict[taxid]:
            rank_name_='gi|{}'.format(tree_dict[taxid]['gi']);
            node_id_=tree_dict[taxid]['gi'];
        else:
            rank_name_=None;
            node_id_=taxid;
    this_node_ = Clade(name=rank_name_, clades=children_);
    if taxid!=1:
        this_node_.taxonomy=Taxonomy(scientific_name=tree_dict[taxid]['rank_name'], rank=rank_, id=Id(taxid, provider='ncbi_taxonomy'));
    return this_node_;

def newickTreeGenerator():
    debug('funct newickTreeGenerator', gap_=True);
    children=newickTreeGeneratorRecurse(taxid=1);
    return Tree(name='gi2tax - Common tree', root= children);

def newickTreeGeneratorRecurse(taxid=1):
    children_=[];

    for child_ in tree_dict[taxid]['children']:
        children_.append(newickTreeGeneratorRecurse(child_));

    if len(children_)==1:
        return children_[0];
    else:
        if taxid == 1:
            rank_name_=None; 
        else:
            if len(children_)==0:
                children_=None;
            if 'gi' in tree_dict[taxid]:
                rank_name_='{}|gi|{}'.format(tree_dict[taxid]['rank_name'], tree_dict[taxid]['gi']);
            else:
                rank_name_=tree_dict[taxid]['rank_name']

            # Remove newick illegal characters from name
            regex=re.compile("[\s,]+"); # One or more whitespace characters
            rank_name_= re.sub(regex, '_', rank_name_);
            rank_name_.replace('(', '[');
            rank_name_.replace(')', ']');

        this_node_ = Clade2(name=rank_name_, clades=children_);
        return this_node_;

###Code

def main(input_path, input_format, output_path, output_format, db_path='./tax_db', db_type='both', update=False, update_db=False, reindex=False, quiet=False):
    ### Setup taxonomic db
    ## Check if files are in local directory or are up to date if not > download them

    if output_path == '-' or output_path == sys.stdout:
        quiet=True;
    if input_path == '-' or input_path == sys.stdin:
        input_path=sys.stdin;


    if not quiet:
        print('####################### Begin gi2tax #######################\n');
        print('Using parameters...');
        print('\tinput path: {}'.format(input_path));
        print('\tinput format: {}'.format(input_format));
        print('\toutput path: {}'.format(output_path));
        print('\toutput format: {}'.format(' '.join(output_format)));
        print('\tdatabase path: {}'.format(db_path));
        if db_type=='protein':  type_='protein gi\'s';
        elif db_type=='nucleotide': type_='nucleotide gi\'s';
        else:   type_='protein and nucleotide gi\'s';
        print('\tdatabase contains: {}'.format(type_));
        print('\tupdate: {}'.format(str(update)));

    if not os.path.isdir(db_path):
        os.makedirs(db_path);

    ###Variable Definitions
    global tax_ftp, tax_ftp_path, req_files, db_path_files, contents;
    tax_ftp='ftp.ncbi.nlm.nih.gov';
    tax_ftp_path='pub/taxonomy';
    req_files={};

    if db_type=='protein' or db_type=='both':
        req_files['gi_taxid_prot.dmp']='gi_taxid_prot.zip';
    if db_type=='nucleotide' or db_type=='both':
        req_files['gi_taxid_nucl.dmp']='gi_taxid_nucl.zip';
    if not os.path.isfile(os.path.join(db_path, 'ncbi_tax.db')) or update or update_db:
        req_files['names.dmp']='taxdmp.zip';
        req_files['nodes.dmp']='taxdmp.zip';
        #req_files['categories.dmp']='taxcat.zip';

    contents={
        'gi_taxid_nucl.zip':['gi_taxid_nucl.dmp'],
        'gi_taxid_prot.zip':['gi_taxid_prot.dmp'],
        'taxcat.zip':['categories.dmp'],
        'taxdmp.zip':['names.dmp', 'nodes.dmp'],
    }; # rest of taxdmp.zip = 'citations.dmp','delnodes.dmp','division.dmp', 'gc.dmp', 'gencode.dmp', 'merged.dmp', 

    db_path_files = findDbPathContents();

    archive_files={};
    if update or ((not os.path.isfile(os.path.join(db_path, 'ncbi_tax.db')) or update_db) and ('names.dmp' not in db_path_files or 'nodes.dmp' not in db_path_files)  and 'taxdmp.zip' not in db_path_files):
        if not quiet:
            print('\n#########  Updating archives   #########');
        archive_files=updateArchives();
        debug("Updated Archives: {}".format(archive_files))

    if not os.path.isfile(os.path.join(db_path, 'ncbi_tax.db')) or update or update_db:
        if len(archive_files) < len(req_files):
            db_path_files = findDbPathContents();
            archive_files = findArchiveFiles(db_path_files, archive_files);
            debug("Updated Archives: {}".format(archive_files));
        success=updateDB(archive_files);

    gis=giFinder(input_path, input_format)
    debug('Found gi\'s: {}'.format(gis));

    if len(archive_files) < len(req_files):
        db_path_files = findDbPathContents();
        archive_files = findArchiveFiles(db_path_files, archive_files);
        debug("Updated Archives: {}".format(archive_files));

    try:
        if 'newick' in output_format or 'phyloXML' in output_format:
            global tree_dict;
            tree_dict={};
        if 'SQLite' in output_format:
            global output_connection, output_cursor;
            if output_path=='-' or output_path==sys.stdout:
                outfile_='gi2tax_output'
            else:
                outfile_=output_path;
            output_connection= sqlite3.connect('{}.db'.format(os.path.splitext(outfile_)[0]));
            output_cursor= output_connection.cursor();
            output_cursor.execute("drop table if exists taxonomy"); #  if the table is already in the database, ie if we are updating the database, delete the previous table.
            output_cursor.execute('''
                create table taxonomy
                (gi integer primary key, taxid integer, subspecies text, species text, subgenus text, genus text, subfamily text, family text, superfamily text, suborder text, order_ text, superorder text, subclass text, class text, superclass text, subphylum text, phylum text, superphylum text, subkingdom text, kingdom text, superkingdom text, no_rank text)
            ''');

        if 'pickle' in output_format or 'json' in output_format:
            global output_dict;
            output_dict={};
        if 'tab' in output_format:
            global out_tab_file;
            if output_path == sys.stdout or output_path == '-':
                out_tab_file=sys.stdout;
            else:
                out_tab_file=open('{}.tab'.format(os.path.splitext(output_path)[0]), 'w');
            column_headings_=['gi', 'taxid', 'subspecies','species','subgenus','genus','subfamily','family','superfamily','suborder','order','superorder','subclass','class','superclass','subphylum','phylum','superphylum','subkingdom','kingdom','superkingdom','no_rank'];
            out_tab_file.write('{}\n'.format('\t'.join(column_headings_)));

        if not quiet:
            print('\n##### Finding Taxonomic information ####')

        if (db_type == 'protein' or db_type == 'both') and (not os.path.isfile(os.path.join(db_path, 'gi_taxid_prot.index')) or reindex):
            gis=giTaxidContainer('gi_taxid_prot.dmp', archive_files, gis, index_=True);
        elif db_type == 'protein' or db_type == 'both':
            gis=giTaxidContainer('gi_taxid_prot.dmp', archive_files, gis, index_=False);
        if (db_type == 'nucleotide' or db_type == 'both') and (not os.path.isfile(os.path.join(db_path, 'gi_taxid_nucl.index')) or reindex):
            gis=giTaxidContainer('gi_taxid_nucl.dmp', archive_files, gis, index_=True);
        elif db_type == 'nucleotide' or db_type == 'both':
            gis=giTaxidContainer('gi_taxid_nucl.dmp', archive_files, gis, index_=False);

        if not quiet:
            print('\nCompleted finding taxonomic information.');
        if len(gis) != 0:
            if not quiet:
                print('Could not find taxonomic information for gis:');
                gi_list=list(gis);
                i=0;
                j=1;
                sys.stdout.write('\t');
                while i<len(gi_list):
                    sys.stdout.write('{} '.format(gi_list[i]));
                    if j == 3:
                        sys.stdout.write('\n\t');
                        j=1;
                    else:   j+=1;
                    i+=1;
                sys.stdout.write('\n');

        if 'newick' in output_format:
            newick_common_tree=newickTreeGenerator();
        if 'phyloXML' in output_format:
            phyloxml_common_tree=phyloxmlTreeGenerator();

    except:
        debug('Error in finding taxonomic path and writing files.');
        debug(sys.exc_info());
        if debug:
            raise;

    finally:
        if not quiet:
            print('\n######  Writing Selected Output   ######');
        if 'tab' in output_format:
            if isinstance(out_tab_file, file) and out_tab_file != sys.stdout:
                out_tab_file.close();
            if not quiet:   print('\tFinished writing tab file');
        if 'SQLite' in output_format:
            output_connection.commit();
            output_connection.close();
            if not quiet:   print('\tFinished writing SQLite database');
        if 'pickle' in output_format:
            if output_path == '-' or output_path == sys.stdout:
                with open('gi2tax_output.pickle', 'w') as outfile_:
                    pickle.dump(output_dict, outfile_, protocol=2); # Write the updated versions dict to the Versions file. NB i'm using protocol 2 for backwards compatibilty with python 2.*
                if not quiet:   print('\tFinished writing pickle file to: gi2tax_output.pickle');
            else:
                with open('{}.pickle'.format(os.path.splitext(output_path)[0]), 'w') as outfile_:
                    pickle.dump(output_dict, outfile_, protocol=2); # Write the updated versions dict to the Versions file. NB i'm using protocol 2 for backwards compatibilty with python 2.*
                if not quiet:   print('\tFinished writing pickle file to: {}'.format('{}.pickle'.format(os.path.splitext(output_path)[0])));
        if 'json' in output_format:
            if output_path == '-' or output_path == sys.stdout:
                with open('gi2tax_output.json', 'w') as outfile_:
                    json.dump(output_dict, outfile_);
                if not quiet:   print('\tFinished writing json file to: gi2tax_output.json');
            else:
                with open('{}.json'.format(os.path.splitext(output_path)[0]), 'w') as outfile_:
                    json.dump(output_dict, outfile_);
                if not quiet:   print('\tFinished writing json file to: {}'.format('{}.json'.format(os.path.splitext(output_path)[0])));
        if 'newick' in output_format:
            if output_path == '-' or output_path == sys.stdout:
                Phylo.write(newick_common_tree, sys.stdout, 'newick');
                if not quiet:   print('\tFinished writing newick file to: stdout');
            else:
                with open('{}.nwk'.format(os.path.splitext(output_path)[0]), 'w') as outfile_:
                    Phylo.write(newick_common_tree, outfile_, 'newick');
            if not quiet:   print('\tFinished writing newick file to: {}'.format('{}.nwk'.format(os.path.splitext(output_path)[0])));
        if 'phyloXML' in output_format:
            if output_path == '-' or output_path == sys.stdout:
                Phylo.write(phyloxml_common_tree, sys.stdout, 'phyloxml');
                if not quiet:   print('\tFinished writing phyloXML file to: stdout');
            else:
                with open('{}.xml'.format(os.path.splitext(output_path)[0]), 'w') as outfile_:
                    Phylo.write(phyloxml_common_tree, outfile_, 'phyloxml');
                if not quiet:   print('\tFinished writing phyloXML file to: {}'.format('{}.xml'.format(os.path.splitext(output_path)[0])));

        if not quiet:
            print('\n####################### End gi2tax #######################\n');

if __name__== '__main__':
    ###Argument Handling
    arg_parser=argparse.ArgumentParser(description='description');
    arg_parser.add_argument("-f", '--input_path', default=None, help="Path to input file(s) or directory, enter '-' for stdin (default). You may indicate multiple files/directories or combinations by separating with a space eg. -f one.faa two/ etc.");
    arg_parser.add_argument("-g", '--input_format', default="regex", help="The format of the files that you are finding gi's from. Default = Use a regular expression pattern", choices=['regex', 'blastXML', 'plain', 'fasta', 'clustal']);    
    arg_parser.add_argument("-o", '--output_path', default=None, help="Path to output file, enter '-' for stdout. Default for tab output is stdout. Default for sqlite, json, and pickle creates an 'gi2tax_output.*' file in working direcory.");
    arg_parser.add_argument("-p", '--output_format', nargs='*',default=["tab"], help="Format(s) to output the taxonomic information. You may specify multiple output formats by separating with a space. Default = tab delimited file.", choices=['tab', 'pickle', 'json', 'newick', 'phyloXML', 'SQLite']);
    arg_parser.add_argument("-d", '--db_path', default="./tax_db/", help="Path to SQLite taxonomic db (this should be a directory/folder). If no db exists at path one will be created. Default= ./tax_db/");
    arg_parser.add_argument("-t", '--db_type', default='both', help="Which set of gi's do you want to map to the taxonomy database. default='both'.", choices=['protein', 'nucleotide', 'both']);
    arg_parser.add_argument("-u", '--update', default=False, action='store_true', help="Boolean toggle. Option to update both the archives and the taxonomic db.");
    arg_parser.add_argument("-v", '--update_db', default=False, action='store_true', help="Boolean toggle. Option to update the sqlite taxonomic db.");
    arg_parser.add_argument("-r", '--reindex', default=False, action='store_true', help="Boolean toggle. Option to update the gi_taxid file index system.");
    arg_parser.add_argument("-q", '--quiet', default=False, action='store_true', help="Boolean toggle. Suppress running feedback. ");
    arg_parser.add_argument('--debug', default=False, action='store_true', help="Boolean toggle. Option to give extensive running feedback to monitor the programs\' progression");
    args = arg_parser.parse_args();

    input_path=args.input_path;
    output_path=args.output_path;
    input_format=args.input_format;
    output_format=args.output_format;
    db_path=args.db_path;
    update=args.update;
    update_db=args.update_db;
    quiet=args.quiet;
    db_type=args.db_type;
    debug=args.debug;
    reindex=args.reindex;

    if debug:
        def debug(string, gap_=False):
            if gap_:
                print('');
            print('DEBUG\t{}'.format(string));
    else:
        def debug(string, gap_=False):
            pass;

    main(input_path, input_format, output_path, output_format, db_path, db_type, update, update_db, reindex,quiet);

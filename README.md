gi-to-tax
=========

what: Takes a group of genbank identifiers (gi's) and finds a full taxonomic path for it from ncbi's taxonomy database.
how: Constructs an SQLite table containing all names and nodes of taxomonic ids. Then Searches through gis to find taxonomic identifiers (using indexing where appropriate) and a recursive function to find the full taxonomic lineage.


Example
-------------------------

	python gi2tax.py -f myBlastResult.xml -g blastXML -t protein -o myBlastResultTaxonomy -p tab newick -u

Downloads all required files from the ncbi taxonomy ftp server, constructs an SQLite database and indexes and searches the protein gi_taxid file for the gi's present in myBlastResult.xml (stores this information in ./tax_db/ (Default)). Writes a tab delimited file of the information and a newick tree to the file myBlastResultTaxonomy.*

	python gi2tax.py -f sequenceFile.fa -g fasta -d my_db_path/ -o myFastaTaxonomy

Searches the for taxonomic information for gi's found in sequenceFile.fa in both the nucleotide and protein set of gi's. Outputs information as a tab delimited file (default) to the path myFastaTaxonomy.tab (.tab added automatically). Required files will not be downloaded unless the taxonomic db doesn't already exist and the required files are not in my_db_path. An SQLite db will not constructed unless one does not already exist in my_db_path. The protein and nucleotide gi_taxid files will not be indexed unless an index file (*.index) does not already exist.

	python gi2tax.py -vrq

Constructs a new SQLite database with the required files present in ./tax_db/ (default for -d). Reindexes both the nucleotide and protein gi_taxid files. Does not attempt to download any new data from ncbi, which protects your local copies from deletion. Does not attempt to find gis and their paths and outputs no information. Does not print any running feedback (quiet, -q). This is useful if you had old copies of the taxonomy files and you just wanted to create a new database for later searching. But could easily search for some information by specifying -f -g -o and -p.

	python gi2tax.py -f myBlastResult.xml -g blastXML -o myBlastResultTaxonomy -p tab phyloXML newick SQLite pickle json -i myTaxids.tab -e myemail@example.com

Searches the for taxonomic information for gi's found in myBlastResult.xml in both the nucleotide and protein set of gi's. gi2tax will attempt to retrieve taxids for gis not found in the local nucleotide or protein database using the entrez api. User defined taxa found in myTaxids.tab will be included in the output. Outputs information as a tab delimited file, SQlite db, Newick and phyloXML trees and as pickle and json data wrappers to myBlastREsultTaxonomy.* (extensions added automatically). Required files will not be downloaded unless the taxonomic db doesn't already exist and the required files are not in my_db_path. An SQLite db will not constructed unless one does not already exist in my_db_path. The protein and nucleotide gi_taxid files will not be indexed unless an index file (*.index) does not already exist.


Options
-------------------------

usage: python gi2tax.py [-f input_path] [-g input_format] [-i input_self_path] [-o output_path] [-p output_format] [-d db_path] [-t db_type] [-u update] [-v update_db] [-r reindex] [-e email] [-q quiet] [--debug]

<table>

	<tr>
	<th>Arg</th><th>Name</th><th>Description</th>
	</tr>
	
	<tr>
	<td>-f</td><td>input_path</td><td>Path to input file, enter '-' for stdin. Default = None (no gis will be searched for in db).</td>
	</tr>
	
	<tr>
	<td>-g</td><td>input_format</td><td>Format of the input file/method to search the file with. Choose one of {'regex', 'blastXML', 'plain', 'fasta'}. Default (regex) uses a regular expression to find gis (matches gi|XXXXXXXX etc.). plain searches for any continuous set of integers. blastXML and fasta search for gis in the name fields with a regex.</td>
	</tr>

	<tr>
	<td>-i</td><td>input_self_path</td><td>Path to tab file containing additional taxids to include. Default is to not include additional taxids. Instructions on how to format your tab file is included in this readme (below)</td>
	</tr>

	<tr>
	<td>-o</td><td>output_path</td><td>Path to output file. Enter '-' for stdout (default). Appropriate extensions will be added automatically based on output format</td>
	</tr>

	<tr>
	<td>-p</td><td>output_format</td><td>Format(s) to output the taxonomic information. Choose from {'tab', 'pickle', 'json', 'newick', 'phyloXML', 'SQLite'}. You may specify multiple output formats by separating with a space. Default = tab delimited file.</td>
	</tr>

	<tr>
		<td>-d</td><td>db_path</td><td>Path to write database information to/search for required files and existing SQLite databases. This should be a directory. If the specified directory doesn't exist it will be created.</td>
	</tr>

	<tr>
		<td>-t</td><td>db_type</td><td>Determines the set of gis to download and search/index. Choose one of {'nucleotide', 'protein', 'both'}. Default = both.</td>
	</tr>

	<tr>
		<td>-u</td><td>update</td><td>Checks the NCBI ftp server for more recent versions of the datasets. Downloads newer versions and DELETES OLD VERSIONS (be careful, the ncbi does not publicly keep backlogs of old files). Also Reconstructs the SQLite3 database with the updated set of files.</td>
	</tr>

	<tr>
		<td>-v</td><td>update_db</td><td>Reconstructs the SQLite db using present local files only. Will always use the newest files available in db_path. If not all required files are found the program will exit with a message telling you whats missing. Useful for working with older datasets as it will preserve your data.</td>
	</tr>
	
	<tr>
	 <td>-v</td><td>update_db</td><td>Reconstructs the SQLite db using present local files only. Will always use the newest files available in db_path. If not all required files are found the program will exit with a message telling you whats missing. Useful for working with older datasets as it will preserve your data.</td>
	</tr>
	
	<tr>
	 <td>-r</td><td>reindex</td><td>Generates a new index for the gi_taxid files (as specified in db_type). You can still search for gi's with this but the process will take longer that with the indexed file. Useful if you manually added a new set of gi_taxid files and didn't delete the *.index files or if you are starting with a new set of gi's and you've maxed out the index size limit (200000 file positions)  from previous searches.</td>
	</tr>
	
	<tr>
		<td>-e</td><td>email</td><td>If email is provided, gi2tax will attempt to find taxids for gis that were not found in the gi_taxid file using the entrez query system. It will only search the databases specified in db_type. The NCBI requires all entrez queries to be accompanied by an email. If they feel that you are overusing the api, they will contact you.</td>
	</tr>

	<tr>
		<td>-q</td><td>quiet</td><td>Runs the program in quiet mode, with no running feedback</td>
	</tr>
	
	<tr>
		<td>-h</td><td>help</td><td>Print help message and exit</td>
	</tr>

	<tr>
		<td>--debug</td><td>debug</td><td>Prints lots of running feedback. Shows output of most function, shows details of handled exceptions and helps you keep track of where the program is at. Useful for debugging. Users shouldn't have to use this. If you are having a problem with the program, please send me an email: darcy.ab.jones@gmail.com</td>
	</tr>

</table>

User specified taxids
-------------------------

The tab delimited user specified taxids file is a plain text file with columns delimited by a single tab and rows as different lines. You can view/edit/save these tab delimited files in most spreadsheet programs.

The first row of the tab file must contain column headings (all lower case). Two columns are required, an id column containing an UNIQUE identifier for the record, and a taxid which corresponds to the organism of interest (doesn't have to be unique). If either one of these columns is missing, the tab file cannot be used. If any record/row is missing information for either of these columns, that record will be excluded from later analysis. If a separate column containing gis is specified and a gi value exists for a row, this will used as instead of id as the unique id, so in the case where an id values is missing for a row but a gi value is specified the record will be included.

Optional columns might include: code, gi, name, description etc. Where possible these row values will be conserved for later output. Row values for these columns can be left empty/unspecified by writing a period (.) where the value would normally be. In the case where a gi has been specified and that gi is already present in the input_path (-f) set of gis, the gi2tax will attempt to merge information from those two records and use the user specified taxid rather than searching for it in the gi_taxid files. 

The purpose of including your own taxids is intended to be used for sequences that do not have a gi but that you still want to include in the output. However, one might use it to avoid having to make entrez queries by including outdated gis with their corresponding taxid or you might just provide the tab output of a previous gi2tax search as user specified taxids in which case the taxids for gis provided may not have to found again. In the latter case, the difference in speed is minimal between searching an indexed gi_taxid file and using the taxid from a user specified tab file, but it might be useful for combining the output of two searched.

Example tab file:

	id taxid gi\n
	gene1 3702  .\n
	gi|20338420|gb|AY093582.1|  192845  20338420\n
	. 3702  17529069\n

The first line contains column headings. For the first record (second line) gi2tax will use gene1 as the unique id and find taxonomic information corresponding to taxid 3702 (Arabidopsis thaliana). The second record contains both an id and a gi, gi2tax will use the gi as the unique identifier and find taxonomic information for taxid 192845 (Erucastrum strigosum). While the third record is missing an id, because it has a value for gi this record is still valid and gi2taxid will use gi as the unique identifier and find taxonomic information for taxid 3702.
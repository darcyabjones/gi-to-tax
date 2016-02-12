gi-to-tax
=========

Takes a group of genbank identifiers (gi's) and finds a full taxonomic path for it from NCBI's taxonomy database.
To make the most of UNIX pipes and avoid rerunning the time consuming steps, gi-to-tax is split into two parts: gi2tax and tax2format.
gi2tax constructs an SQLite table containing all names and nodes of taxomonic ids, then Searches through gis to find taxonomic identifiers (using indexing where appropriate) and a recursive function to find the full taxonomic lineage, with output as a [line delimited JSON file](http://jsonlines.org/).
tax2format takes the JSON output of gi2tax and can produce output in other convenient formats such as [tsv](https://en.wikipedia.org/wiki/Tab-separated_values), [csv](https://en.wikipedia.org/wiki/Comma-separated_values), [newick tree](https://en.wikipedia.org/wiki/Newick_format), or [phyloXML tree](http://www.phyloxml.org/).


Setup
-----

To use gi2tax, you'll need to have the (uncompressed) NCBI taxonomy database somewhere that you can access them.
Specifically gi2tax needs the `names.dmp`, `nodes.dmp`, and one or both of the `gi_taxid_nucl.dmp` or `gi_taxid_prot.dmp` files (depending on the options used).

In Linux (and probably Macintosh) you can easily download, check, and decompress the files as follows:

	wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.gz.md5
	md5sum -c gi_taxid_nucl.dmp.gz.md5
	gunzip gi_taxid_nucl.dmp.gz

	wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz.md5
	md5sum -c gi_taxid_prot.dmp.gz.md5
	gunzip gi_taxid_prot.dmp.gz

	wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
	wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
	md5sum -c taxdump.tar.gz.md5
	tar -zxf taxdump.tar.gz

All of the required files need to be in the same directory, which you will point gi2tax to at runtime.
When needed, gi2tax will create an SQLite database of the `names.dmp` and `nodes.dmp` files (from the taxdump archive) , and a JSON formatted index of each of the `gi_taxid_nucl.dmp` or `gi_taxid_prot.dmp` files, in the same directory as they were found.
This unfortunately means that you need write permission to that folder.
If this is a problem (eg when using compute clusters) I suggest sym-linking the required taxonomy files to somewhere that you can write.

Please note that because gis are not permanent, they are sometimes not found in the `gi_taxid` files because they have been replaced.
If possible, use NCBI taxonomic information from the same time that your homology search was performed or sequences were retrieved.
Because this is not always possible, gi2tax provides a way to find taxids for remaining gis (including depreciated ones) using the Entrez api.
You will however, need to provide an email address (a requirement for Entrez API use) and have an active web connection.


Examples
-------

	gi2tax.py -i my_blast_result.xml -g blastxml -t protein -o my_blast_result_taxonomy.json

Searches the protein gi_taxid file for the gis present in `my_blast_result.xml`.

	gi2tax.py -i sequence_file.fa -g fasta -d my_db_path -o sequence_file_taxonomy.json

Searches the for taxonomic information for gis found in `sequence_file.fa` in both the nucleotide and protein sets of gis.

	gi2tax.py -i my_taxids.json -g json -e myemail@example.com > my_taxids_taxonomy.json

Finds taxonomic information for a set of gis or taxids from a line delimited json file.
This could contain records that are not associated with a gi, but have a known taxid (eg. from a draft genome sequence; see below for a detailed explanation on how to use this feature).
gi2tax will attempt to retrieve taxids for gis not found in the local nucleotide or protein database using the NCBI Entrez api.

	tax2format.py -i sequence_file_taxonomy.json -f tsv -o sequence_file_taxonomy.tsv

Converts line-delimited JSON output from gi2tax into a spreadsheet friendly tab separated values format.

User specified taxids
-------------------------

Two special input formats are provided to allow users to get information for records without gis.
Users can provide a line delimited JSON file (or a JSON file with a root array) or a TSV file as input.
Each record must contain either a 'gi', or an 'id' and a 'taxid'.
The TSV format must have a header row specifying the order of the fields.

The following files are examples of valid input with the same information.

```
{"id": "my_sequence", "taxid": 5206}
{"gi": 111111}
{"gi": 1202, "taxid": 999}
```

```
[
	{"id": "my_sequence", "taxid": 5206},
	{"gi": 111111},
	{"gi": 1202, "taxid": 999},
]
```

```
id	gi	taxid
my_sequence		5206
	111111
	1202	999
```

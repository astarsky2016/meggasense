# meggasense
Pipeline for metagenome reads functional assembly and annotation.
DISCLAIMER:
This pipeline was not intended as a stand alone application!
Scripts and programs commited were created by SemGen Ltd. as part of "Amylomics" FP7 funded project.
Pipeline was tested on OS Linux Debian version 8.6
Prerequisites: installed Newbler assebler, BLAST suite of programs, HMMER v > 3.0, PfamScan and glimmer-mg
Analysis data stored in MongoDB, sqlite, MySQL and python ZODB OO database (stores sequence reads and all steps of subsequent analyses)
Requires more than 1 CPU (relies on multiprocessing).
Pipeline relies on KEGG for it's annotations. Since KEGG stopped providing it's data for free, only data from obsolete KEGG v 58 are provided!

Download ftp for KEGG data v 58 used in meggasense pipeline:




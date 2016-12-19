# meggasense
Pipeline for metagenome reads functional assembly and annotation.
DISCLAIMER:
This pipeline was not intended as a stand alone application!
Scripts and programs commited were created by SemGen Ltd. as part of "Amylomics" FP7 funded project (http://cordis.europa.eu/result/rcn/157198_en.html).
Pipeline was intended for annotation of metagenome WGS on Roche 454 platform (.sff raw data) with special emphasis on carbohydrate acting enzymes (CAZymes).
Pipeline was tested on OS Linux Debian version 8.6!

Prerequisites: installed Newbler assebler, BLAST suite of programs, HMMER v > 3.0, PfamScan and glimmer-mg. There is a series of easily obtainable python modules such as pyparsing, zodb, sqlite etc. Third party software such as glimmer-mg and PfamScan.pl (in Perl) was heavily modified and cannot be distributed as such. We can share the modified source upon direct request. In this case please contact us directly.

Contact: 
http://semgen.info, contact@semgen.info

Analysis data is being stored in MongoDB, sqlite, MySQL and python ZODB OO database (stores sequence reads and all steps of subsequent analyses).
Once set up, pipeline requires more than 1 CPU (relies on multiprocessing).
Pipeline relies on KEGG for it's annotations. Since KEGG stopped providing it's data for free, only data from obsolete KEGG v 58 are provided!

Download link for KEGG data v 58 used in meggasense pipeline:

http://bioserv.pbf.hr/kegg_58_hmm.tar.gz

Download link for processed CAZy data:

http://bioserv.pbf.hr/cazy.tar.gz






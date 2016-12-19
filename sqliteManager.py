'''
Created on Aug 30, 2012

@author: astar
'''
import apsw
class sqliteManager(object):
    def __init__(self, db):
        self.memcon=apsw.Connection(":memory:")
        self.disc = db
        self.disc_connection=apsw.Connection(db)
        # Copy into memory
        with self.memcon.backup("main", self.disc_connection, "main") as backup:
            backup.step() # copy whole database in one go
    def alterTable(self, columns):
        for column in columns:
            self.memcon.cursor().execute('ALTER TABLE hits ADD COLUMN %s INTEGER' % column)
    def close(self):
        self.memcon.close()
    def disc_close(self):
        self.disc_connection.close()
    def getContigCount(self, ko, library):
        ids = [str(row[0]) for row in self.memcon.cursor().execute("SELECT sequences.id FROM sequences, hits WHERE hits.keggid = '%s' AND sequences.id = hits.seqid AND sequences.assembly_type = 'assembled_read' AND sequences.metaid = '%s'" % (ko, library))]
        contigs = set()
        [contigs.add(str(row[0])) for row in self.memcon.cursor().execute("SELECT contigid FROM assembly WHERE readid IN (%s) AND kegg = '%s'" % ((", ".join(ids)), ko))]
        return len(contigs)
    def getReadCount(self, ko, library):
        row = self.memcon.cursor().execute("SELECT COUNT (sequences.id) FROM sequences, hits WHERE hits.keggid = '%s' AND sequences.id = hits.seqid AND sequences.metaid = '%s'" % (ko, library)).next()
        return row[0]
    def addColumnValues(self, contig, library, updates):
        for row in self.memcon.cursor().execute("SELECT id FROM sequences WHERE seqid = '%s' AND metaid = '%s'" % (contig, library)):
            assert len(row) == 1, "more than 1 entry %s in library %s!!!" % (contig, library)
            seq = row[0]
            values = tuple(updates) + (seq,)
            self.memcon.cursor().execute("UPDATE hits SET start = %i, stop = %i WHERE keggid = '%s' AND seqid = %i" % values)
    def query(self, queryStr):
        return self.memcon.cursor().execute(queryStr)
    def backupToDisc(self):
        destination = apsw.Connection(self.disc)
        # Copy onto disc
        with destination.backup("main", self.memcon, "main") as backup:
            backup.step() # copy whole database in one go
        print "backed up to: %s" % self.disc
        self.memcon.close()
        destination.close()
                                                           
    

if __name__ == "__main__":
    # There will be no disk accesses for this query
    count = 0
    x = sqliteManager("/media/Transcend/TAGS/Jura_meta/reads/metagenomes.sqlite")
    print x.getContigCount('K00186', '927-2E_test_reads')
    print x.getReadCount('K00186', '927-2E_test_reads')
    x.memcon.close()

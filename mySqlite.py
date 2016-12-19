'''
Created on Dec 2, 2011

@author: astar
'''
import sqlite3, sys, logging
logging.basicConfig(level=logging.DEBUG)
class MySqlite(object):
    def __init__(self, db):
        try:
            self.conn = sqlite3.connect(db)
        except sqlite3.OperationalError: # Can't locate database file
            sys.exit(1)
        self.cursor = self.conn.cursor()
     
    def createTable(self, cmd):
        self.cursor.execute(cmd)
        self.conn.commit()
    
    def execute(self, cmd, records):
        try:
            self.cursor.execute(cmd, records)
        except Exception, ex:
            logging.exception("problem with SQL execution!")
            print cmd
            for rec in records:
                print rec      
            sys.exit(1)
        
    def commit(self): 
        self.conn.commit()
    
    def getRecords(self, cmd, ids):
        self.cursor.execute(cmd, ids)
        return(self.cursor.fetchall())
    
    def closeHandle(self):
        '''Closes the connection to the database'''
        self.conn.commit() # Make sure all changes are saved
        self.conn.close()
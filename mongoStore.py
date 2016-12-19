'''
Created on Dec 2, 2011

@author: astar
'''
import sys
from pymongo import Connection
from pymongo.errors import ConnectionFailure
class MongoConnect(object):
    def __init__(self, db):
        """connect to MongoDB"""
        try:
            c = Connection(host="localhost", port=27017)
        except ConnectionFailure, exc:
            sys.stderr.write("Could not connect to MongoDB: %s" % exc)
            sys.exit(1)
        #get a database handle
        self.dbh = c[db]
        assert self.dbh.connection == c
        print "Successfully set up a MongoDB database handle"
if __name__ == "__main__":
    conn = MongoConnect("myTestDB")
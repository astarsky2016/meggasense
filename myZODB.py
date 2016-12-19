'''
Created on Dec 2, 2011

@author: astar
'''
from ZODB import FileStorage, DB
import transaction

class MyZODB(object):
    count = 0
    def __init__(self, path):
        MyZODB.count += 1
        self.storage = FileStorage.FileStorage(path)
        self.db = DB(self.storage)
        self.connection = self.db.open()
        self.dbroot = self.connection.root()
        if MyZODB.count == 100:
            print "Performing routine packing of ZODB database: %s please have patience" % path
            self.db.pack()
            MyZODB.count = 0

    def close(self):
        self.connection.close()
        self.db.close()
        self.storage.close()
        
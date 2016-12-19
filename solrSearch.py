'''
Created on Feb 23, 2012

@author: astar
'''
from urllib2 import *
import simplejson
conn = urlopen('http://localhost:1234/solr/core0/select/?q=protein&wt=json')
rsp = simplejson.load(conn)
print rsp

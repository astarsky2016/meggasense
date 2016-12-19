'''
Created on Oct 13, 2013

@author: astar
'''
from bs4 import BeautifulSoup
from Bio import SeqIO
import MySQLdb as mdb
import urllib2, collections, time, shelve, sys
from Bio import Entrez
Entrez.email = "astar@pbf.hr"
def dd():
    return collections.defaultdict(list)
data = {}
org_list = set()
taxon = shelve.open("taxon.db")
taxDict = collections.defaultdict(set)
for key, val in taxon.items():
    for org in val:
        if org not in ['Eukaryotae', 'Eucarya', 'Eucaryotae', 'Eukarya', 'eukaryotes', 'eucaryotes', 'biota', 'cellular organisms', 'all', 'root']:
            taxDict[org].add(key)
taxon.close()
print "taxa loaded!"
source = "http://www.cazy.org/"
data["GH"] = "http://www.cazy.org/Glycoside-Hydrolases.html"
data["GT"] = "http://www.cazy.org/GlycosylTransferases.html"
data["PL"] = "http://www.cazy.org/Polysaccharide-Lyases.html"
data["CE"] = "http://www.cazy.org/Carbohydrate-Esterases.html"
data["AA"] = "http://www.cazy.org/Auxiliary-Activities.html"
data["CBM"] = "http://www.cazy.org/Carbohydrate-Binding-Modules.html"
con = mdb.connect('localhost', 'cazyuser', 'lubbert', 'cazydb')
con.set_character_set('utf8')
# missing_from_nr = []
missed = collections.defaultdict(list)
gbKeys = collections.defaultdict(dd)
nr_keys = {}
CAZyFasta = collections.defaultdict(set)
def is_coherent(seq):
    for pos, el in enumerate(range(seq[0], seq[-1]+1)):
        if seq[pos] != el:
            return(pos)   
    return(len(seq) + 1)    
def batchEntrez(connection):
    cursor = connection.cursor() 
    print "missing", len(missed)
    tryAgain = []
    orgs = set()
#     out_handle = open("missing_entries.fasta", "w")
    missing_list = missed.keys()
#     batch = 0
    for missing_acc in missing_list:
    #for i in xrange(0, len(missing_list), 100):
        downloaded = False
        while not downloaded:
            try:   
                #net_handle = Entrez.efetch(db="protein", id=", ".join(missing_list[i:i+100]), rettype="fasta", retmode="text")
                net_handle = Entrez.efetch(db="protein", id=missing_acc, rettype="fasta", retmode="text")
                data = net_handle.read()
                store = []
                for el in missed[missing_acc]:
                    store.append(list(el[:4]) + [data] + list(el[4:]))
                missed[missing_acc] = store
                #out_handle.write(data)
#                 batch += 100 
#                 print batch
#                 out_handle.write("\n")  
                net_handle.close()
                downloaded = True
            except Exception,e: 
                print str(e)
                print "skipping"
                tryAgain.append(missing_acc)
                break
    times = 0
    while tryAgain and times < 100:
        times += 1
        for acc in tryAgain:
            try:
                net_handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
                data = net_handle.read()
                missed[acc] = missed[acc][:4] + [data] + missed[acc][4:]
#                 out_handle.write(data)
#                 out_handle.write("\n")
                net_handle.close()
                tryAgain.remove(acc)
            except:
                print acc, "will try again!"
   # for value in missed.values():
    for missing_acc in missing_list:
        for entry in missed[missing_acc]:
            if len(entry) == 7:
                cursor.execute('INSERT INTO Families (acc, family_id, class_id, orgn_name, fasta, ec_number, structure) VALUES (%s, %s, %s, %s, %s, %s, %s)', entry)
            else:
                print "kurslus"
                print missed[missing_acc]
            if entry[3] not in org_list:
                org_list.add(entry[3])
                orgs.add(entry[3])
    print orgs
    getTaxon(list(orgs), cursor)
#     out_handle.close()
#     for seq in SeqIO.parse(open("missing_entries.fasta"), "fasta"):
#         missing_from_nr.append((seq.description, seq.format("fasta")))
def clearHash():
    global databases
    databases = {"emb": [], "dbj": [], "gb": [], "ref": [], "tpg": [], "tpe": [], "tpd": [], "pdb": [], "gi": []}
def hasTheRightParent(tag):
    return tag.parent.name == "select" and tag.parent["name"] == "cazypage"
def getAllAcc(bowl, key, element):
    for tag in bowl.find_all("tr"):
        accFlag = False
        pdb = False
        acc = None
        ec = None
        taxn = None
        for child in tag.descendants:
            if child.name == "b":
                link = child.parent.get('href')
                if link:
                    if link.find("entrez") != -1:
                        acc = child.string
                        accFlag = True
                    if link.find("cazy.org") != -1:
                        taxn = child.string.rstrip().lstrip()
            if child.name == "a":
                link = child.get('href')
                if link:
                    if link.find("ec=") != -1:
                        ec = child.string
                    if link.find("structureId") != -1:
                        pdb = True
                    if link.find("Taxonomy/Browser") != -1:
                        taxn = child.string.rstrip().lstrip()
        if accFlag: 
            gbKeys[key][element].append((acc, ec, pdb, taxn))
def findAcc(tags):
    for db in databases.keys():
        for pos in [i + 1 for i, x in enumerate(tags) if x == db]:
            databases[db].append(tags[pos])
# def taxoDict():
#     if os.path.isfile("gi2tax.db"):
#         return
#     else:
#         gi2tax = shelve.open("gi2tax.db")
#         for line in open("gi_taxid_prot.dmp"):
#             data = [el.lstrip().rstrip() for el in line.split()]
#             gi2tax[data[0]] = data[1]
#         gi2tax.close()    
def getTaxon(batch_orgn, cursor):
    print "entering get taxon!"
    btch = []
    if batch_orgn:
        for org in batch_orgn: 
            taxId = "-1"
            bacteria = False
            archaea = False
            eukaryota = False
            unclassified = False
            name = org
            names = list(taxDict.get(org, ["-1"]))
            if len(names) != 1:
                print "###########################"
                print names
                print org
                taxId = ",".join(names)
                print "###########################"
            else:
                taxId = names[0]
                if taxId in taxDict.get("Archaea", []): archaea = True
                if taxId in taxDict.get("Bacteria", []): bacteria = True
                if taxId in taxDict.get("Eukaryota", []): eukaryota = True
                if taxId in taxDict.get("Unclassified", []): unclassified = True 
            btch.append((taxId, name, bacteria, archaea, eukaryota, unclassified))
    cursor.executemany('INSERT INTO Organisms (tax_id, name, bacteria, archaea, eukaryota, unclassified) values (%s, %s, %s, %s, %s, %s)', 
                   btch)   
    print "exiting get taxon!"
    return(1)
            
def initiateCazy(connection):  
#     taxon = shelve.open("gi2tax.db")
    cursor = connection.cursor()         
    batch_families = []
    batch_organisms = []
    limit = 5000
    for family in gbKeys.keys():
        for cazy in gbKeys[family].keys():
            if gbKeys[family][cazy]:
                accs = [el[0] for el in gbKeys[family][cazy]]
                ec_s = [el[1] for el in gbKeys[family][cazy]]
                structures = [el[2] for el in gbKeys[family][cazy]]
                taxonomy = [el[3] for el in gbKeys[family][cazy]]
                accToall = {}
                for pos, acc in enumerate(accs):
                    accToall[acc] = (ec_s[pos], structures[pos], taxonomy[pos])
                cursor.execute('SELECT id, fasta FROM nr WHERE id IN (' + ", ".join(["'" + acc + "'" for acc in accs]) + ')')
                data = cursor.fetchall()
                fastas = [el[1]for el in data]
                nr_accs = [el[0] for el in data]
                #gis = [el[2] for el in data]
#                 taxid = collections.defaultdict(int)
#                 for taxa in [taxon.get(el, None) for el in gis]:
#                     taxid[taxa] += 1
#                 representativeTax = sorted(taxid, reverse=True)[0]
                for pos, fasta in enumerate(fastas):
                    faa = fasta
                    if fasta.find(">") == -1: # CREATE TABLE Families family_id VARCHAR(25) PRIMARY KEY, class_id VARCHAR(25), tax_id VARCHAR(25), fasta VARCHAR(25), ec_number VARCHAR(25), structure BOOLEAN, gis text'
                        cursor.execute('SELECT fasta FROM nr WHERE id = ' + "'" + fasta + "'")
                        faa = cursor.fetchall()[0][0]
#                         cursor.execute('INSERT INTO Families (acc, family_id, class_id, tax_id, fasta, ec_number, structure) VALUES (%s, %s, %s, %s, %s, %s, %s)', 
#                                        (nr_accs[pos], cazy, family, accToall[nr_accs[pos]][2], fasta, accToall[nr_accs[pos]][0], accToall[nr_accs[pos]][1]))
#                         #CAZyFasta["_".join((family, cazy)) + ".faa"].add(fasta)
#                     else:
#                         
#                         cursor.execute('SELECT fasta FROM nr WHERE id = ' + "'" + fasta + "'")
#                         faa = cursor.fetchall()[0][0] acc VARCHAR(25), family_id VARCHAR(25), class_id VARCHAR(25), tax_id VARCHAR(100), fasta TEXT, ec_number VARCHAR(25), structure BOOLEAN'
                    batch_families.append((nr_accs[pos], cazy, family, accToall[nr_accs[pos]][2], faa, accToall[nr_accs[pos]][0], accToall[nr_accs[pos]][1]))
                    #if len(batch_families) == limit:
                    #    cursor.executemany('INSERT INTO Families (acc, family_id, class_id, orgn_name, fasta, ec_number, structure) values (%s, %s, %s, %s, %s, %s, %s)', 
                    #                   batch_families)
                        #batch_families = []
                        #print "batch families inserted!"
                        #CAZyFasta["_".join((family, cazy)) + ".faa"].add(faa)
                    if accToall[nr_accs[pos]][2] not in org_list:
                        org_list.add(accToall[nr_accs[pos]][2])
                   #     batch_organisms.append(accToall[nr_accs[pos]][2])
                   #     if len(batch_organisms) == limit:
                   #         status = getTaxon(batch_organisms, cursor)
                   #         batch_organisms = []
                   #         print "batch organisms inserted"
                for missing_acc in set(accs).difference(set(nr_accs)):
                    missed[missing_acc].append((missing_acc, cazy, family, accToall[missing_acc][2], accToall[missing_acc][0], accToall[missing_acc][1]))
    #if batch_families: cursor.executemany('INSERT INTO Families (acc, family_id, class_id, orgn_name, fasta, ec_number, structure) values (%s, %s, %s, %s, %s, %s, %s)',
    #                                   batch_families)
    #if batch_organisms: 
    #    status = getTaxon(batch_organisms, cursor)
            #missed.update(set(accs).difference(set(nr_accs)))
# def batchUpdate():
#     pass
#     for family in gbKeys.keys():
#         for cazy in gbKeys[family].keys():
#             if gbKeys[family][cazy]:
#                 for acc in gbKeys[family][cazy]:
#                     for batchRetrieved in missing_from_nr:
#                         if batchRetrieved[0].find(acc) != -1:
#                             CAZyFasta["_".join((family, cazy)) + ".faa"].add(batchRetrieved[1])
# def finalize():
#     for family, fastas in CAZyFasta.items():
#         if fastas:
#             f = open(family, "w")
#             for fasta in fastas:
#                 f.write(fasta + "\n")
#             f.close()
def create_nr(connection):
    cur = connection.cursor()
    cur.execute("DROP TABLE nr") 
    cur.execute("SHOW TABLES")
    if "nr" not in [table[0] for table in cur.fetchall()]:
        print "loading cazy database with nr data..."
        cur.execute('CREATE TABLE nr (id VARCHAR(20), fasta MEDIUMTEXT, gis TEXT)')
        clearHash()
        batchCount = 0
        batch = []
        for seq in SeqIO.parse(open("nr"), "fasta"):
            batchCount += 1
            tags = seq.description.replace("\x01", "|").split("|")
            findAcc(tags)
            first = True
            pointer = None
            for accList in databases.values():
                if accList:  
                    gis = ",".join(databases.get("gi", None))
                    for acc in accList:
                        if first:
                            batch.append((acc, seq.format("fasta"), gis))
                            pointer = acc
                            first = False
                        else:
                            batch.append((acc, pointer, None))
            if batchCount % 50000 == 0: 
                cur.executemany('INSERT INTO nr (id, fasta, gis) values (%s, %s, %s)', batch) 
                print "batch %i commited!" % batchCount
                batch = []
            clearHash()
        if batch:
            cur.executemany('INSERT INTO nr (id, fasta, gis) values (%s, %s, %s)', batch) 
        print "nr data loaded"
    cur.execute('CREATE INDEX idx ON nr (id)')    
def create_cazy(connection):
    cur = connection.cursor()
    cur.execute("SHOW TABLES")
    if "Cazy_classes" not in [table[0] for table in cur.fetchall()]:
        cur.execute('CREATE TABLE Cazy_classes (class_id VARCHAR(25) PRIMARY KEY, date_created DATE)')
        cur.execute('CREATE TABLE Families (acc VARCHAR(25), family_id VARCHAR(25), class_id VARCHAR(25), orgn_name TEXT, fasta MEDIUMTEXT, ec_number VARCHAR(25), structure BOOLEAN)')
        cur.execute('CREATE TABLE Organisms (tax_id VARCHAR(100), name TEXT, bacteria BOOLEAN, archaea BOOLEAN, eukaryota BOOLEAN, unclassified BOOLEAN)')
# def index_cazy(cur):
#     cur.execute()
if __name__ == "__main__":
    with con:
        cur = con.cursor()
        cur.execute('SET NAMES utf8')
        cur.execute('SET CHARACTER SET utf8')
        cur.execute('SET character_set_connection=utf8')
        cur.execute("SELECT VERSION()")
        ver = cur.fetchone()
        print "Database version : %s " % ver
        #create_nr(con)
        #cur.execute("DROP TABLE Cazy_classes")
        #cur.execute("DROP TABLE Families")
        #cur.execute("DROP TABLE Organisms")
        #create_cazy(con)
        for key, item in data.items():
            #cur.execute('INSERT INTO Cazy_classes (class_id, date_created) VALUES (%s, %s)', (key, time.strftime("%Y-%m-%d")))
            response = urllib2.urlopen(item)
            soup = BeautifulSoup(response)
            for cazy_type in soup.find_all(hasTheRightParent):
                portion = BeautifulSoup(urllib2.urlopen(cazy_type["value"]))
                for link in portion.find_all('a'):
                    if link.string == "All":
                        spoon = BeautifulSoup(urllib2.urlopen(link.get("href")))
                        pages = [(1, link.get("href"))]
                        nextPage = True
                        while nextPage:
                            appending = []
                            for page in spoon.find_all('a'):
                                if page.get("href"):
                                    if page.get("href").find("#pagination") > 0 and page['rel'][0] == 'nofollow': 
                                        if int(page.string) > pages[-1][0]:
                                            appending.append((int(page.string), source + page.get("href")))
                            if appending:
                                appending.sort(key = lambda tup: tup[0])
                                diff = is_coherent([el[0] for el in appending])
                                pages += appending[0:diff]
                                spoon = BeautifulSoup(urllib2.urlopen(pages[-1][1]))
                            else:
                                nextPage = False
                        for page in pages:
                            getAllAcc(BeautifulSoup(urllib2.urlopen(page[1])), key, cazy_type.string)
                        print cazy_type.string                         
    #taxoDict()
        initiateCazy(con)
        print "CAZy initiated!"
        batchEntrez(con)
        print "Entrez batch loaded!" 
#     print "Missing %i sequences fetched" % len(missing_from_nr)
#     batchUpdate()
        print "CAZy updated!"
    #finalize()
        print "finished"
    con.commit()     
    con.close()
         
            
    
    
        
                
            

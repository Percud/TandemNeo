import gzip
from Bio import SeqIO
import os
import re
import pandas as pd
import json

cwd = os.path.dirname(os.getcwd())

if 'database.json' in os.listdir(cwd):
    js = json.load(open(cwd + '/database.json'))
    jskeys = js.keys()

class database():   
    
    def suffixes():
        """return suffixes used by ensembl to distinguish 
        between species with format. 
        e.g.: {'ENSG00': 'Homo_sapiens}"""

        suffixes = {}
        path = cwd + '/fa/'
        for f in os.listdir(path):
            if f.endswith('.fa.gz'):
                s = f.split('.fa')[0].capitalize()
                file = gzip.open(path + f).readline()
                l = file.split()[0][1:7]
                l = l.decode("utf-8")
                suffixes.update({s: l})
                suffixes.update({l: s})
                
        suffixes.update({'ENSG00': 'Homo_sapiens'})

        return suffixes 
    
        #database.suffixes()
    
    def whichdup(acc):
        """return wich kind of duplication 
        for a given ensembl protein ID"""

        def read(kind):
            path = cwd + '/orthologues/' + kind + '.csv'
            return pd.read_table(path, sep=';')

        t = read('Tandem').isin([acc]).any()
        d = read('Divergent').isin([acc]).any()
        c = read('Convergent').isin([acc]).any()

        if any(t):
            return 'Tandem'
        elif any(d):
            return 'Divergent'
        elif any(c):
            return 'Convergent'
        else:
            return None

        #database.whichdup('ENSGALP00000072350')
        
    def dataframe(acc):
        """return a dataframe containing 
        a given ensembl protein ID"""

        kind = database.whichdup(acc)
        path = cwd + '/orthologues/' + kind + '.csv'
        df = pd.read_table(path, sep=';')
        df = df.set_index(['Classes', 'Orders', 'Species'])

        return df
        
        #database.dataframe('ENSGALP00000072350')

    def pairs(acc):
        """return the orthogroups pair identifier 
        for a given ensembl protein ID"""

        df = database.dataframe(acc)
        p = df.columns[df.isin([acc]).any()][0]
        p_num = re.split('A|B', p)[0]

        if 'B' in p:
            pairs = [p_num + 'A', p]
        else:
            pairs = [p, p_num + 'B']

        return [p, pairs]

        #database.pairs('ENSGALP00000072350')

    def orthodf(k):
        """return the given duplication 
        kind orthologues dataframe"""

        path = cwd + '/orthologues/' + k + '.csv'
        df = pd.read_table(path, sep=';')
        return df
    
        # database.orthodf('tandem')

    def acklist(s, k):
        """return an accessions list of 
        a given specie and kind of duplications"""

        df = database.orthodf(k)
        df = df[df['Species'] == s]
        l = df.dropna(axis=1).values.tolist()[0][3:]
        return l
    
       # database.acklist('Homo_sapiens', 'tandem')

    def aclist(s):
        """return an accession list of
        a given specie of any kind of duplications"""

        d = ['tandem', 'convergent', 'divergent']
        l = [database.acklist(s, d[0]) + 
             database.acklist(s, d[1]) +
             database.acklist(s, d[2])]
        
        return l[0]
        
        # database.aclist('Homo_sapiens')

    def faline(q, file):
        """return the FASTA format strings for 
        a given ensembl protein ID in a given file,
        usefull for speed up parsing processes"""
        
        if not file:
            s = database.suffixes().get(q[0:6])
            path = cwd + '/fa/' + s + '.fa.gz'
            handle = gzip.open(path, 'rt')
            file = list(SeqIO.parse(handle, "fasta"))
        
        alist = [l.id.split('.')[0] for l in file]
        i = alist.index(q)

        return file[i]
        
        #database.faline('ENSGALP00000072350', file.fa)

    def forquery(q, file):
        """return a dictionary containing: protein ID, gene ID, 
        transcript ID, symbol, product, sequence, sequence lenght, 
        orthogroup for a given ensembl protein ID,
        can give a file for speed up the process 
        or let the system find it"""
        
        if not file:
            s = database.suffixes().get(q[0:6])
            path = cwd + '/fa/' + s + '.fa.gz'
            handle = gzip.open(path, 'rt')
            file = list(SeqIO.parse(handle, "fasta"))

        l = database.faline(q, file)
        x = re.split('[.:\s]', l.description)
        pid = x[0]
        gid = x[x.index('gene')+1]
        tid = x[x.index('transcript')+1]
        try:
            sym = x[x.index('gene_symbol')+1]
        except:
            sym = None
        try:
            fr = x.index('description')+1
            to = x.index('[Source')-1
            pro = ' '.join(x[fr:to])
        except:
            pro = None
        try:
            ort = database.pairs(pid)[0]
        except:
            ort = None

        main = {
            'protein': pid,
            'gene': gid,
            'transcript': tid,
            'symbol': sym,
            'product': pro,
            'sequence': l.seq,
            'lenght': len(l.seq),
            'orthogroup': ort
        }

        return main
    
        #database.forquery('ENSGALP00000072350', file.fa (or None))
        
    def info(l):
        """parsing of a given FASTA entry"""

        x = re.split('[.:\s]', l.description)

        try:
            fr = x.index('description')+1
            to = x.index('[Source')-1
            pro = ' '.join(x[fr:to])
        except:
            pro = None

        try:
            sym = x[x.index('gene_symbol')+1]
        except:
            sym = None
        
        db = ({x[0]:
                {
                    'transcript id': x[x.index('transcript')+1],
                    'gene id': x[x.index('gene')+1],
                    'symbol': sym,
                    'product': pro,
                    'sequence': str(l.seq),
                    'sequence lenght': len(l.seq)
                }
        })
        
        return db
        
        #database.info(FASTA entry)
        
class dbinfo:
    """create an object class 
    to speedup future parsing of an entry"""

    def __init__(self, pid):
        
        self.pid = pid
        self.gid = js.get(pid)['gene id']
        self.sym = js.get(pid)['symbol']
        self.tid = js.get(pid)['transcript id']
        self.pro = js.get(pid)['product']
        self.seq = js.get(pid)['sequence']
        self.sle = js.get(pid)['sequence lenght']
        #self.pair = db.pairs(pid)[0]
        #self.let = db.pairs(pid)[0][-1]
        #self.num = re.split('A|B', db.pairs(pid)[0])[0]
        #self.spe = db.suffixes().get(pid[0:6])
        
        # a = 'ENSP0918398'
        # ainfo = dbinfo(a)
        # ainfo.gid

    def ids_to(specie, kind, ID_kind):
        """IDs conversion"""
    
        ortholist = database.acklist(specie, kind)
        
        if ID_kind == 'protein':
            ortholist = ortholist
        elif ID_kind == 'gene':
            ortholist = [js.get(l)['gene id'] for l in ortholist if l in js.keys()] 
        elif ID_kind == 'transcript':
            ortholist = [js.get(l)['transcript id'] for l in ortholist if l in js.keys()]
        
        return ortholist

        # dbinfo.ids_to('homo_sapiens', 'convergent', 'gene')
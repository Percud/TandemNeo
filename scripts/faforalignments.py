import os
import pandas as pd
import json
from database import database as db
from database import dbinfo
import re

cwd = os.path.dirname(os.getcwd())

class faforalignments:
    
    def pairslist(kind):
        """return a list containing all indexes 
        of orthologues table for a given 
        kind of duplication"""
        
        df = db.orthodf(kind)
        cols = df.columns.tolist()[3:]
        pairs = [list(a) for a in 
                 zip([l for l in cols 
                      if 'A' in l], 
                     [l for l in cols 
                      if 'B' in l])]

        return pairs
    
        # alignemnts.pairslist('tandem')
        
    def printfa(df, p, suff, json, d):
        """print the FASTA file for a 
        given orthogroup"""

        for x in [0, 1]:
            ogroup = p[x]
            number = p[x][0:-1]
            letter = p[x][-1]

            accessions = df[p[x]].values.tolist()
            for a in accessions:

                try:
                    pid    = a
                    gid    = json.get(a)['gene id']
                    seqlen = json.get(a)['sequence lenght']
                    sym    = json.get(a)['symbol']
                    pro    = json.get(a)['product']
                    seq    = json.get(a)['sequence']
                    spe    = suff.get(a[0:6])
                    num    = number
                    tag    = letter

                    path = cwd + '/alignments/' + d + '/'
                    f = open(path + num + '.fa', 'a')

                    print('>', end='', file=f)

                    print(*[
                        pid, 
                        gid, 
                        seqlen, 
                        sym, 
                        pro, 
                        spe, 
                        ogroup,
                        num, 
                        tag
                        ], sep='\t', file=f)

                    print(seq, file=f)
                except:
                    pass
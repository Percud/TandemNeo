# /usr/bin/env python

import requests, sys, os, pandas as pd, numpy as np, re, gzip, itertools, json, csv, shutil, xmltodict, urllib.request
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast.Applications import NcbimakeblastdbCommandline as makeblastdb
from multiprocessing import Pool
from ftplib import FTP
from Bio import SeqIO
import natsort as ns
from natsort import natsorted
from Bio import Align, AlignIO, SeqIO, Phylo
from Bio.Align.Applications import ClustalOmegaCommandline as clustalo
from Bio.Align import substitution_matrices, MultipleSeqAlignment as MSA, AlignInfo
from Bio.PDB import PDBParser as PDB, PDBIO
from itertools import combinations
from random import sample

cwd = os.getcwd()

# DUPLICATIONS

for file in os.listdir(cwd + '/blast_queries'):
    if file.endswith('blast.txt'):

        # apro il file blast e seleziono solo le coppie che non siano se stesse contro se stesse e che siano best reciprocal hit
        
        x = pd.read_table(cwd + '/blast_queries/' + file, comment='#', header=None, names=['Query', 'Subject', 'Identity', 'A.lenght', 'Mismatches', 'Gap', 'Q.start', 'Q.end', 'S.start', 'S.end', 'Evalue', 'Bits'])
        x = x[x['Identity'] != 100].drop_duplicates('Query')
        x['Query'] = x['Query'].apply(lambda x: x.split('.')[0])
        x['Subject'] = x['Subject'].apply(lambda x: x.split('.')[0])
        columnab = [[line[0], line[1]] for line in x.values.tolist()]
        columnabbrh = [line for line in [[line[0], line[1]] for line in x.values.tolist()] if [line[1], line[0]] in columnab]
        xbrhdf = pd.merge(x, pd.DataFrame(columnabbrh), left_on=['Query', 'Subject'], right_on=[0, 1]).drop(columns=0)
        xy = pd.merge(pd.read_table(cwd + '/main/' + file.split('_blast')[0] + '.tsv'), xbrhdf, right_on='Query', left_on='Protein_id', how='outer').drop(columns=[1, 'Query'])
        xy = xy.values.tolist()

        # scrivo tre liste con le duplicazioni in tandem che rispettano le seguenti condizioni: strand uguali, se strand '+' allora lo stop del primo deve essere minore dello start del secondo, se '-' allora lo start del primo deve essere minore dello stop del secondo, se invece gli strand sono diversi e se lo strand del primo è '+' allora lo stop del primo deve essere minore dello stop del secondo mentre se lo strand del primo è '-' allora lo start del primo deve essere minore dello strand del secondo
        
        tandem, convergent, divergent = [], [], []
        for x in range(len(xy)-1):
            strand, start, stop = xy[x][2], xy[x][4], xy[x][5]
            strand2, start2, stop2 = xy[x+1][2], xy[x+1][4], xy[x+1][5]
            query, subject = xy[x][3], xy[x+1][8]
            if query == subject:
                if strand == strand2:
                    if (strand == '+' and stop < start2) or (strand == '-' and start < stop2):
                        tandem.append([xy[x][1], xy[x][3], xy[x+1][1], xy[x+1][3], xy[x][9], xy[x][17]])
                if strand != strand2:
                    if strand == '+' and stop < stop2:
                        convergent.append([xy[x][1], xy[x][3], xy[x+1][1], xy[x+1][3], xy[x][9], xy[x][17]])
                    if strand == '-' and start < start2:
                        divergent.append([xy[x][1], xy[x][3], xy[x+1][1], xy[x+1][3], xy[x][9], xy[x][17]])
               
        # scrivo tre file.tsv per le coppie: in tandem, divergenti, convergenti
        for line in tandem:
            print(*line, sep='\t', file=open(cwd + '/duplications/' + file.split('_blast')[0] + '_tandem.tsv', 'a'))
        for line in divergent:
            print(*line, sep='\t', file=open(cwd + '/duplications/' + file.split('_blast')[0] + '_divergent.tsv', 'a'))
        for line in convergent:
            print(*line, sep='\t', file=open(cwd + '/duplications/' + file.split('_blast')[0] + '_convergent.tsv', 'a'))

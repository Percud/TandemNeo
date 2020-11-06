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

# FASTA FOR ALIGNMENTS

# apro il database scritto con lo script numero 6
database = json.load(open('database.json'))

# creo le tre cartelle (una per ogni classe di duplicazione genica) che conterranno i FASTA per gli allineamenti
if not os.path.exists(cwd + '/alignments/tandem_alignments') and not os.path.exists(cwd + '/alignments/convergent_alignments') and not os.path.exists(cwd + '/alignments/divergent_alignments'):
    os.mkdir(cwd + '/alignments/tandem_alignments'), os.mkdir(cwd + '/alignments/convergent_alignments'), os.mkdir(cwd + '/alignments/divergent_alignments')
    
# scrivo un file FASTA per ciascun ortogruppo contenente tutte le sequenze
def fastaforaligments(kind):
    data = pd.read_csv(cwd + '/orthologues/' + kind + '_orthologues.csv', sep=';').set_index(['Classes', 'Orders', 'Species'])
    pairs_ab = list(zip([line for line in data.columns.tolist() if 'A' in line], [line for line in data.columns.tolist() if 'B' in line]))

    for x in range(len(pairs_ab)):
        for f in data[list(pairs_ab[x])].values.tolist():

            def getfastas(symbol, product, lenght, gene, specie, tandem, pid, fasta):
                    print('>' + pid + '\t' + gene + '\t' + lenght + '\t' + symbol + '\t' + product + '\t' + specie + '\t' + tandem + '\n' + fasta, file=open(cwd + '/alignments/' + kind + '_alignments' + '/' + re.split('A|B', tandem)[0] + '.fa', 'a'))            

            if not 'nan' in str(f[0]): 
                getfastas(database.get(f[0])[8], database.get(f[0])[7], 
                          str(len(database.get(f[0])[0])), database.get(f[0])[1], 
                          database.get(f[0])[9], data[list(pairs_ab[x])].columns.tolist()[0], 
                          f[0], database.get(f[0])[0])
                
            if not 'nan' in str(f[1]):
                getfastas(database.get(f[1])[8], database.get(f[1])[7], 
                          str(len(database.get(f[1])[0])), database.get(f[1])[1], 
                          database.get(f[1])[9], data[list(pairs_ab[x])].columns.tolist()[1], 
                          f[1], database.get(f[1])[0])
            
fastaforaligments('tandem')
fastaforaligments('convergent')
fastaforaligments('divergent')
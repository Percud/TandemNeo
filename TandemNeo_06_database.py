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

# DATABASE

# apro le tre tabelle di ortologia delle tre classi di duplicazioni
tandem_orthologs = pd.read_table(cwd + '/orthologues/' + 'tandem_orthologues.csv', sep=';').set_index(['Classes', 'Orders'])
convergent_orthologs = pd.read_table(cwd + '/orthologues/' + 'convergent_orthologues.csv', sep=';').set_index(['Classes', 'Orders'])
divergent_orthologs = pd.read_table(cwd + '/orthologues/' + 'divergent_orthologues.csv', sep=';').set_index(['Classes', 'Orders'])
orthologs = tandem_orthologs.values.tolist() + divergent_orthologs.values.tolist() + convergent_orthologs.values.tolist()

# aggiungo ad un dizionario le seguenti informazioni: FASTA, protein_ID, gene_ID, transcript_ID, cromosoma, strand, coordinate di start e stop, simbolo genico, prodotto, specie, lunghezza della sequenza
data = {}
for x in range(len(orthologs)):
    specie = orthologs[x][0] 
    
    main = {}
    for fasta in SeqIO.parse(gzip.open(cwd + '/fa/' + specie + '.fa.gz', 'rt'), 'fasta'):

        splitted = re.split('\s', fasta.description)
        pid = splitted[0].split('.')[0]
        gid = splitted[3].split('gene:')[1].split('.')[0]
        tra = splitted[4].split('transcript:')[1].split('.')[0]
        chrom = splitted[2].split(':')[2]
        strand = splitted[2].split(':')[5]
        if strand == '1':
            sta, sto = splitted[2].split(':')[3], splitted[2].split(':')[4]
        elif strand == '-1':
            sta, sto = splitted[2].split(':')[4], splitted[2].split(':')[3]
        try:
            symbol = splitted[7].split('gene_symbol:')[1]
        except:
            symbol = 'no_data'
        try:
            if '[Source' in fasta.description:
                product = re.split('description:', fasta.description)[1].split(' [Source')[0]
            else:
                product = re.split('description:', fasta.description)[1]
        except:
            product = 'no_data'

        main.update({pid: [str(fasta.seq), gid, chrom, strand, int(sta), int(sto), len(fasta.seq), product.strip(), symbol, specie]})

    for line in orthologs[x]:
        if not 'nan' in str(line):
            data.update({line: main.get(line)})
            
# salvo il dizionario in formato json
json.dump(data, open('database.json', 'w')) 

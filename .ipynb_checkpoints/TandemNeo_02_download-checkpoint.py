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


# DOWNLOAD

# apro la lista delle specie
species = [line.split()[2] + '_' + line.split()[3] for line in open(species_file_path)]

# annoto i codici degli assemblies e delle release
assemblies = json.load(open('assemblies.json'))
release = [line['release'] for line in assemblies['species'] if line['name'] == 'homo_sapiens']

def getassembly(specie):
    """return assembly identifier for a given species name"""
    
    for assembly in assemblies['species']:
        if specie.lower() == assembly['name']:
            return assembly['assembly']
        
def download(specie):
    """download faa from ensembl ftp"""
    
    ftp = FTP('ftp.ensembl.org')
    ftp.login()
    ftp.cwd('/pub/current_fasta/' + specie.lower() + '/pep/')
    ftp.retrbinary('RETR ' + specie + '.' + getassembly(specie) + '.pep.all.fa.gz', open(cwd + '/fa/' + specie + '.fa.gz', 'wb').write)
    if specie in queries:
        ftp.cwd('/pub/current_gtf/' + specie.lower() + '/')
        ftp.retrbinary("RETR " + specie + '.' + getassembly(specie) + '.' + str(release[0]) + '.gtf.gz', open(cwd + '/gtf/' + specie + '.gtf.gz', 'wb').write)
    ftp.quit()
        
# attivo la funzione download in multiprocessing
if __name__ == '__main__':
    pool = Pool(processes=20)
    pool.map(download, [specie for specie in species])

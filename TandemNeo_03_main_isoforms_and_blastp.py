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

# MAIN ISOFORMS AND BLASTP

species = [line.split()[2] + '_' + line.split()[3] for line in open('species_list.txt')]
queries = ['Homo_sapiens', 'Gallus_gallus', 'Danio_rerio']

assemblies = json.load(open('assemblies.json'))
from TandemNeo_02_download import getassembly

# per ogni query
for line in queries:
    
    # scarico da appris un dataframe contenente le isoforme principali per la data specie
    appris = pd.read_table(urllib.request.urlopen("http://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/" + line.lower() + "/" + getassembly(line.lower()) + "/appris_data.principal.txt"), header=None)
    appris = appris[appris[4].str.contains('PRINCIPAL')].sort_values([1, 4]).drop_duplicates(1)
    for line2 in appris.values.tolist():
        print(*line2, sep='\t', file=open(cwd + '/appris/' + line + '_principal_appris.txt', 'a'))
    main_isoforms = appris[2].values.tolist()

    # scremo i fasta per le isoforme principali e scrivo un file che le contiene
    for fasta in SeqIO.parse(gzip.open(cwd + '/fa/' + line + '.fa.gz', 'rt'), 'fasta'):
        if 'gene_biotype:protein_coding' in fasta.description and re.split('\s', fasta.description)[4].split('transcript:')[1].split('.')[0] in main_isoforms:
            print('>' + fasta.id + '\n' + fasta.seq, file=open(cwd + '/main/' + line + '_main.fa', 'a'))

    # faccio lo stesso per i gtf
    gtf = pd.read_table(cwd + '/gtf/' + line + '.gtf.gz', compression='gzip', comment='#', header=None)
    gtf = gtf.join(gtf[8].str.split('"', expand=True).add_prefix('sec'))           
    gtf = gtf[gtf['sec5'].isin(main_isoforms)]

    # scrito un file.tsv con le diverse informazioni inerenti agli accession e soprattutto alle coordinate 
    main = []
    for line3 in gtf.values.tolist():
        if 'CDS' in line3[2]:
            chromosome, strand, gene_id, transcript_id, gene_name = line3[0], line3[6], line3[10], line3[14], line3[20]
            protein_id = [line for line in line3 if str(line).startswith(gene_id.split('0')[0][:-1] + 'P')]
        elif 'start' in line3[2]:
            start = line3[3]
        elif 'stop' in line3[2]:
            stop = line3[4]
            main.append([chromosome, gene_id, strand, protein_id[0], start, stop, transcript_id, gene_name])

    tsv = pd.DataFrame(main, columns=['Chromosome', 'Gene_id', 'Strand', 'Protein_id', 'Start', 'Stop', 'Transcript_id', 'Gene_name'])
    tsv['Chromosome'] = pd.Categorical(tsv['Chromosome'], ordered=True, categories=ns.natsorted(tsv['Chromosome'].unique()))
    tsv['Mean'] = tsv[['Start', 'Stop']].mean(axis=1)
    tsv = tsv.sort_values(['Chromosome', 'Mean'], ascending=(True, True)).drop(columns='Mean').to_csv((cwd + '/main/' + line + '.tsv'), index=False, sep='\t')
    
    # blast intra-specie
    makeblastdb(dbtype='prot', input_file=(cwd + '/main/' + line + '_main.fa'))()
    blastp(query=(cwd + '/main/' + line + '_main.fa'), db=(cwd + '/main/' + line + '_main.fa'), num_threads=num_threads, evalue=evalue, max_hsps=max_hsps, outfmt=outfmt, max_target_seqs=max_target_seqs, out=cwd + '/blast_queries/' + line + '_blast.txt')()

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
import statistics as st

cwd = os.getcwd()

# ALIGNMENTS

def conservation(msa): 
    """ assegna uno score ad ogni posizione dell'allineamento con sum of pairs """
    
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load(matrix)
    numbers, columns, scores = [], [], []
    for n in range(msa.get_alignment_length()):
        columns.append((msa[:,n]))
        numbers.append(n+1)
    for c in columns:
        c = list(c)
        if 'X' in c:
            c.remove('X') # rimuove le X dalle sequenze
        pairs = list(combinations(c,2))
        score = []
        try:
            for p in pairs:
                if not '-' in p: # assegna valore 0 per i gap
                    score.append(aligner.score(p[0],p[1]))
            scores.append(sum(score)/len(pairs))
        except:
            pass
    return (list(zip(numbers,columns,scores)))

def getalignmentjson(folder):
    """return coverage, number of sequences, alignment lenght and positions for a given aligned fasta file"""

    alignments_log = {}

    # define filename    
    number_of_files = max([int(file.split('.')[0]) for file in os.listdir(cwd + '/alignments/' + folder) if file.endswith('.fa')])
    for name in range(number_of_files):
        name = str(name)
        fa = cwd + '/alignments/' + folder + '/' + name + '.fa'
        fasta = cwd + '/alignments/' + folder + '/' +  name + '.fasta'

        len_A, len_B = [], []
        for f in list(SeqIO.parse(fa, 'fasta')):
            if 'A' in f.description.split('\t')[-1]:
                len_A.append(int(f.description.split('\t')[2]))
            elif 'B' in f.description.split('\t')[-1]:        
                len_B.append(int(f.description.split('\t')[2]))

        coverage = round(min(st.mean(len_A),st.mean(len_B))/max(st.mean(len_A),st.mean(len_B))*100,2)
        if coverage > min_coverage:

            clustalo(infile = fa, outfile = fasta, force = False)()
            # sampling sequences from alignments
            alignment = AlignIO.read(fasta, 'fasta')

            dif_score, A, B, gap = [], [], [], []
            for a in alignment:
                if 'A' in a.description.split('\t')[-1]:
                    A.append(a) 
                elif 'B' in a.description.split('\t')[-1]:
                    B.append(a)
            A, B = sample(A, min(len(A),len(B))), sample(B, min(len(A),len(B)))
            aln_A, aln_B, aln = MSA(A), MSA(B), MSA(A+B)
            cons_A, cons_B = AlignInfo.SummaryInfo(aln_A).dumb_consensus(), AlignInfo.SummaryInfo(aln_B).dumb_consensus()

            differences = list(zip(conservation(aln), conservation(aln_A), conservation(aln_B)))   
            hit, res_A, res_B = [], [], []
            for d in differences:
                dif_score.append([d[0][0], d[1][1], d[2][1], round(min(d[1][2], d[2][2])-d[0][2],2)])
                hit.append(d[0][0])
                res_A.append(cons_A[d[0][0]-1])
                res_B.append(cons_B[d[0][0]-1])

            pos_A, pos_B, pos_AB = [], [], []
            for n in range(len(alignment)):
                for res in zip(res_A,hit,res_B):
                    if ensp in alignment[n].name: # first 5 letters of reference specie's protein_ID (i.e. ENSP0)
                        try:
                            pos_AB.append([res[0], res[2]])
                            if 'B' in alignment[n].description.split('\t')[-1][-1]:
                                pos_B.append(res[1]-list(alignment[n][0:res[1]-1]).count('-'))
                            else:
                                pos_B.append(None)
                            if 'A' in alignment[n].description.split('\t')[-1][-1]:
                                pos_A.append(res[1]-list(alignment[n][0:res[1]-1]).count('-'))
                            else:
                                pos_A.append(None)
                        except:
                            pass

            x = pd.concat([pd.Series(pos_A).dropna().reset_index().astype(np.int64), pd.Series(pos_B).dropna().reset_index().astype(np.int64)], axis=1).drop(columns='index').astype(str)
            y = zip(x.values.tolist(), pos_AB)
            ratio = round(len(hit)/np.mean(min(len_A,len_B))*100, 2)

            positions = {}
            for p in list(zip(dif_score, y)):
                key = p[0][0]
                column_A = p[0][1]
                column_B = p[0][2]
                score = p[0][3]
                homo_position_A = p[1][0][0]
                homo_position_B = p[1][0][1]
                consensus_A = p[1][1][0]
                consensus_B = p[1][1][1]
                position = {'Column A': column_A, 'Column B': column_B, 'Score': score, 
                            'Homo position A': homo_position_A, 'Homo position B': homo_position_B, 
                            'Consensus A': consensus_A, 'Consensus B': consensus_B}
                positions.update({key: position})

            accessions = []
            if any(ensp in string for string in [fasta.id for fasta in SeqIO.parse
                                                     (cwd + '/alignments/' + folder + '/' + name + '.fa', 'fasta')
                                                     if 'A' in fasta.description.split('\t')[6]]) == False:
                accessions.append(None)
            if any(ensp in string for string in [fasta.id for fasta in SeqIO.parse
                                                     (cwd + '/alignments/' + folder + '/' + name + '.fa', 'fasta')
                                                     if 'B' in fasta.description.split('\t')[6]]) == False:
                accessions.append(None)
            for fasta in SeqIO.parse(cwd + '/alignments/' + folder + '/' + name + '.fa', 'fasta'):
                if ensp in fasta.id and 'A' in fasta.description.split('\t')[6]:
                    accessions.append(fasta.id) 
                if ensp in fasta.id and 'B' in fasta.description.split('\t')[6]:
                    accessions.append(fasta.id) 


            alignments_log.update({name: {'Accessions': accessions, 'Coverage': coverage, 'Number of sequences': len(list(SeqIO.parse(fa, 'fasta'))), 'Alignment lenght': aln.get_alignment_length(), 'Ratio': ratio, 'Positions': positions, 'Scores': [v['Score'] for k,v in positions.items()]}})


    json.dump(alignments_log, open(cwd + '/alignments/' + folder + '/alignments.json', 'w'))
    

getalignmentjson('tandem_alignments')
getalignmentjson('convergent_alignments')
getalignmentjson('divergent_alignments')
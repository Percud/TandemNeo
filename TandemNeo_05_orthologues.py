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

# COMPARA ORTHOLOGUES

def getorthoss(line):
    """return a column of orthologues for a given query"""
    
    def getorthoss_letter(query, count, letter, specie, types):

        try:

            # scarico in formato json la lista degli ortologhi secondo COMPARA
            r = requests.get("https://rest.ensembl.org/homology/id/" + query + "?type=orthologues;aligned=0;cigar_line=0;compara=vertebrates", headers={"Content-Type": "application/json"})
            ortho = r.json()

            # elaboro il json per scrivere un file testuale, comodo per le fasi successive, contenente tutti una lista di ortologhi per ogni dato accession e le relative informazioni sulla percentuale d'identità, accession genici, proteici e sulle specie in cui sono stati trovati
            for x in range(len(ortho['data'][0]['homologies'])):
                target = ortho['data'][0]['homologies'][x]['target']
                per, pid, gid, spe = target['perc_id'], target['protein_id'], target['id'], target['species']
                source = ortho['data'][0]['homologies'][x]['source']
                per2, pid2, gid2, spe2 = source['perc_id'], source['protein_id'], source['id'], source['species']
                print(pid, gid, spe.capitalize(), count, letter, per, specie, str(count) + letter, file=open(cwd + '/orthologues/' + specie + '_' + types + '_orthologues.txt', 'a'))
            print(pid2, gid2, spe2.capitalize(), count, letter, 100, specie, str(count) + letter, file=open(cwd + '/orthologues/' + specie + '_' + types + '_orthologues.txt', 'a'))

        except:
            print(query, 'no_data', 'no_data', count, letter, 0, specie, str(count) + letter, file=open(cwd + '/orthologues/' + specie + '_' + types + '_orthologues.txt', 'a'))
    
    getorthoss_letter(line.split()[0], counter.get(line.split()[0]), 'A', line.split()[6].capitalize(), line.split()[7])
    getorthoss_letter(line.split()[2], counter.get(line.split()[2]), 'B', line.split()[6].capitalize(), line.split()[7])
       
# assegnazione di un progressivo per ogni coppia: ad esempio la prima coppia di accession sarà salvata come 1A per il partner 1 della coppia e 1B per il partner 2 della coppia   
counter = {} 
for file in os.listdir(cwd + '/duplications'):
    species = file.split('_')[0] + '_' + file.split('_')[1]
    types = file.split('_')[2].split('.')[0]
        
    for n, line in enumerate(open(cwd + '/duplications/' + file)):
        counter.update({line.split()[0]: n}), counter.update({line.split()[2]: n})
        
    # attivo la funzione per scaricare gli ortologhi in multiprocessing, rispettando le tempistiche dei server ENSEMBL
    if __name__ == '__main__':
        pool = Pool(processes=4)
        pool.map(getorthoss, [line + '\t' + species + '\t' + types for line in open(cwd + '/duplications/' + file)])

def orthologues_df(kind):

    # creo un dizionario con la tassonomia a partire dalla lista delle specie (specie, ordine, classe)
    classes, orders = {}, {}
    for line in open('species_list.txt'):
        classes.update({line.split()[2] + '_' + line.split()[3]: line.split()[0]})
        orders.update({line.split()[2] + '_' + line.split()[3]: line.split()[1]})

    # apro le tre tabelle
    h = pd.read_table(cwd + '/orthologues/Homo_sapiens_' + kind + '_orthologues.txt', sep='\s', engine='python', header=None)
    d = pd.read_table(cwd + '/orthologues/Danio_rerio_' + kind + '_orthologues.txt', sep='\s', engine='python', header=None)
    g = pd.read_table(cwd + '/orthologues/Gallus_gallus_' + kind + '_orthologues.txt', sep='\s', engine='python', header=None)

    # aggiungo lettere corrispondenti alle specie queries
    h[7] = 'H.' + h[7]
    d[7] = 'D.' + d[7]
    g[7] = 'G.' + g[7]

    # concateno le tre tabelle una dietro l'altra
    nobrh = pd.concat([h,d,g])

    # aggiungo voci della tassonomia
    nobrh[8] = nobrh[2].apply(lambda x: classes.get(x))
    nobrh[9] = nobrh[2].apply(lambda x: orders.get(x))
    
    
    
    ref_classes = ['Sauropsida', 'Mammalia', 'Actinopteri']
    
    
    
    nobrh = nobrh[nobrh[8].isin(ref_classes)]
    nobrh[10] = nobrh[7].apply(lambda x: re.split('\.', x)[0])
    nobrh = nobrh.dropna(subset=[8])

    # numeriamo le hit di compara per ogni query per ogni specie (one-to-many)
    # se con gene di gallus da 20 ortologhi, numera questi ortologhi da 1 a 20
    nobrhsort = nobrh.sort_values([10, 3, 4, 2, 5], ascending=(False, True, True, True, False))
    nobrhsort[11] = nobrhsort.groupby([2, 7]).cumcount()+1

    # drop duplicati stessa percentuale d'identità per ogni query per ogni specie per il sorting
    nobrhsort = nobrhsort.drop_duplicates([2, 5, 7])

    # sorting per percentuale d'identità globale, drop geni duplicati tenendo il primo (one-to-one), drop stessa specie stesso gruppo 
    # "best-hit"
    bh = nobrhsort.sort_values(5, ascending=False).drop_duplicates(1).drop_duplicates([2,7])

    # tiene solo le prime hit di one-to-one
    # come nei casi in cui da ortologo come prima hit quando in realtà era il terzo
    # "reciprocal"
    brh = bh[bh[11] == 1]
    brh = brh.sort_values([10, 3, 4], ascending=(False, True, True))
    brh[20] = brh.groupby(7, sort=False)[7].transform("count")

    counts = brh[[7, 20]].drop_duplicates(7).values.tolist()
    couple_high_then3 = []
    for x in range(len(counts)-1):
        if counts[x][1] and counts[x+1][1] > 3 and re.split('A|B', counts[x][0])[0].split('.')[1] == re.split('A|B', counts[x+1][0])[0].split('.')[1]:
            couple_high_then3.append(counts[x][0])
            couple_high_then3.append(counts[x+1][0])
    brhfull = brh[brh[7].isin(couple_high_then3)]
    brhfull = brhfull.reset_index().drop(columns=['index'])
    brhfull[10] = brhfull.groupby([3,6], sort=False, as_index=False).ngroup()
    brhfull[10] = brhfull[10].astype(str) + brhfull[4]

    # pivot table
    brhfullpivot = brhfull.pivot(index=2, columns=10, values=0)
    brhfullpivot = brhfullpivot.reset_index()
    brhfullpivot['Classes'] = brhfullpivot[2].apply(lambda x: classes.get(x))
    brhfullpivot['Orders'] = brhfullpivot[2].apply(lambda x: orders.get(x))
    brhfullpivot = brhfullpivot.rename(columns={2: 'Species'})
    brhfullpivot = brhfullpivot.set_index(['Classes', 'Orders', 'Species']).sort_index()
    brhfullpivot = brhfullpivot.reindex(natsorted(brhfullpivot.columns), axis=1)
    brhfullpivot.to_csv(cwd + '/orthologues/' + kind + '_orthologues.csv', sep=';')
    
# attivo tutto il procedimento per le tre classi di duplicazioni 
orthologues_df('tandem')
orthologues_df('divergent')
orthologues_df('convergent')

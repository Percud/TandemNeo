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
import configparser

cwd = os.getcwd()


# WORKING DIRECTORIES

if not os.path.exists(cwd + '/fa') and not os.path.exists(cwd + '/gtf') and not os.path.exists(cwd + '/main') and not os.path.exists(cwd + '/blast_queries') and not os.path.exists(cwd + '/duplications') and not os.path.exists(cwd + '/orthologues') and not os.path.exists(cwd + '/appris') and not os.path.exists(cwd + '/alignments') and not os.path.exists(cwd + '/uniprot'):
    os.mkdir(cwd + '/fa'), os.mkdir(cwd + '/gtf'), os.mkdir(cwd + '/main'), os.mkdir(cwd + '/blast_queries'), os.mkdir(cwd + '/duplications'), os.mkdir(cwd + '/appris'), os.mkdir(cwd + '/orthologues'), os.mkdir(cwd + '/alignments'), os.mkdir(cwd + '/uniprot')

# SPECIES LIST

config = configparser.ConfigParser()
config.optionxform = str
config.read('config.txt')
locals().update(dict(config.items('DEFAULT')))


if species_file == '':

    # scarico gli assembly da ENSEMBL al fine di recuperare le versioni dei fasta e dei gtf e per scrivere i percorsi delle cartelle cartelle dell'ftp
    assemblies = requests.get("http://rest.ensembl.org/info/species?", headers={"Content-Type": "application/json"}).json()
    json.dump(assemblies, open('assemblies.json', 'w'))

    # recupero da NCBI le informazioni sulla tassonomia e sul numero di pubblicazioni per una lista di specie
    def getaxonomy(specie):
        """return class, order, taxid, publications for a give specie name"""

        try:

            eutils = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
            taxid = xmltodict.parse(requests.get(eutils + 'esearch.fcgi?db=taxonomy&term=' + specie).content)['eSearchResult']['IdList']['Id']
            pubs = xmltodict.parse(requests.get(eutils + 'esearch.fcgi?db=pubmed&term=' + specie).content)['eSearchResult']['Count']
            taxons = xmltodict.parse(requests.get(eutils + 'efetch.fcgi?db=taxonomy&id=' + taxid + '&retmode=xml&rettype=full').content)['TaxaSet']['Taxon']['LineageEx']['Taxon']

            scientific_name = [taxons[x]['ScientificName'] for x in range(len(taxons)) if taxons[x]['ScientificName'] in ['Sauropsida', 'Mammalia', 'Actinopteri']]
            orders = [taxons[x]['ScientificName'] for x in range(len(taxons)) if taxons[x]['Rank'] == 'order']

            return scientific_name + orders + [specie.split('_', 1)[0]] + [specie.split('_', 1)[1]] + [taxid] + [pubs]

        except:
            pass

    # creo una lista delle specie disponibili e interessanti per il progetto
    available_species = []
    for assembly in assemblies['species']:
        taxonomy = getaxonomy(assembly['name'].capitalize())
        if not taxonomy == None and len(taxonomy) == 6:
            available_species.append(taxonomy)

    # filtro le specie: 5 per famiglia, genere unico
    species_list_df = pd.DataFrame(available_species)
    species_list_df[5] = species_list_df[5].astype(int)
    species_list_df_uniquegenre = species_list_df.sort_values(5, ascending=False).groupby(2, sort=False).head(1)
    species_list_df_cap = species_list_df_uniquegenre.groupby(1, sort=False).head(12).sort_values([0, 1])

    # scrivo il file con la lista delle specie e le informazioni correlate
    for line in species_list_df_cap.values.tolist():
        print(*line, sep=' ', file=open('species_list.txt', 'a'))

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
import urllib.parse
import urllib.request

cwd = os.getcwd()

def features(kind):

    # apro tabella ortologhi e creo lista degli accession di Homo da analizzare
    orthos = pd.read_table(cwd + '/orthologues/' + kind + '_orthologues.csv', sep=';').set_index(['Classes', 'Orders', 'Species'])
    homo = pd.DataFrame(orthos.loc[('Mammalia', 'Primates', 'Homo_sapiens')])
    homo_joined = ' '.join(map(str, [line[0] for line in homo.values.tolist()]))

    # preparo gli accession per la conversione Ensembl_ID --> Uniprot_ID
    params = {'from': 'ENSEMBL_PRO_ID', 'to': 'ACC', 'format': 'tab', 'query': homo_joined}
    url = 'https://www.uniprot.org/uploadlists/'
    req = urllib.request.Request(url, urllib.parse.urlencode(params).encode('utf-8'))
    with urllib.request.urlopen(req) as f:
        response = f.read()
        # 3 non li trova (su 1054)

    # scrivo un dizionario di conversione
    converted = {}
    for line in response.decode('utf-8').split('\n'):
        if not 'From' in line and not line.split() == []:
            converted.update({line.split()[1]: line.split()[0]})

    # recupero le features da uniprot in formato json
    def getfeatures(uniprot_id):
        """return the Uniprot complete dataset for a given Uniprot ID"""
        try:
            r = requests.get("https://www.ebi.ac.uk/proteins/api/proteins?size=-1&accession=" + uniprot_id, headers={ "Accept" : "application/json"})
            data = pd.json_normalize(r.json())
            return data
        except:
            return str(uniprot_id) + ' not_found'

    # concateno le features in formato tabulare e le scarico in locale
    appended_data = []
    for k, v in converted.items():
        data = getfeatures(k)
        appended_data.append(data)
    appended_data = pd.concat(appended_data)
    appended_data.to_csv(cwd + '/uniprot/' + kind + '_raw_features.csv')

    raw_features = pd.read_table(cwd + '/uniprot/' + kind + '_raw_features.csv', sep=',')
    raw_features = raw_features[raw_features['accession'].isin(list(converted.keys()))]

    # organizzo e filtro le features in base agli accession nel tabellone di ortologia, unendole in un unico dataframe
    dataframe_list = []
    for line in raw_features[['accession', 'gene', 'protein.recommendedName.ecNumber', 'protein.recommendedName.fullName.value', 'features']].iterrows():
        if not str(line[1]['features']) == 'nan' and 'description' in list(eval(line[1]['features'])[0].keys()):
            dataframe = pd.DataFrame(eval(line[1]['features']))[['type', 'category', 'description', 'begin', 'end']]
            dataframe['accession'] = line[1]['accession']
            if not str(line[1]['gene']) == 'nan':
                dataframe['gene'] = eval(line[1]['gene'])[0]['name']['value']
            else:
                dataframe['gene'] = 'NaN'
            if not str(line[1]['protein.recommendedName.ecNumber']) == 'nan':
                dataframe['EC_number'] = eval(line[1]['protein.recommendedName.ecNumber'])[0]['value']
            else:
                dataframe['EC_number'] = 'NaN'
            dataframe['full_protein_name'] = line[1]['protein.recommendedName.fullName.value']
            dataframe['ensembl_accession'] = converted.get(line[1]['accession'])
            dataframe_list.append(dataframe)

    df = pd.concat(dataframe_list)[['accession', 'ensembl_accession', 'gene', 'full_protein_name', 'type', 'category', 'description', 'EC_number', 'begin', 'end']]

    orthos = pd.read_table(cwd + '/orthologues/' + kind + '_orthologues.csv', sep=';')
    orthos = list(zip(orthos[orthos['Species'] == 'Homo_sapiens'].values.tolist()[0][3:], orthos[orthos['Species'] == 'Homo_sapiens'].columns.tolist()[3:]))
    orthos = pd.DataFrame(orthos, columns = ['ensembl_accession', 'orthogroup'])
    df = pd.merge(df, orthos, on='ensembl_accession')
    df = df[['accession', 'orthogroup', 'ensembl_accession', 'gene', 'full_protein_name', 'type', 'category', 'description', 'EC_number', 'begin', 'end']]

    df.to_csv(cwd + '/uniprot/' + kind + '_filtered_features.csv')

    # scrivo i range di posizioni che mi interessano (in questo caso -1 e +1 sia per begin che per end)
    df = pd.read_table(cwd + '/uniprot/' + kind + '_filtered_features.csv', sep=',').drop(columns='Unnamed: 0').reset_index()
    df['begin'] = df['begin'].apply(lambda x: None if '~' in str(x) else x)
    df['end'] = df['end'].apply(lambda x: None if '~' in str(x) else x)
    df['b+1'] = df['begin'].apply(lambda x: str(int(x)+1) if not x == None else x)
    df['b-1'] = df['begin'].apply(lambda x: str(int(x)-1) if not x == None else x)
    df['e+1'] = df['end'].apply(lambda x: str(int(x)+1) if not x == None else x)
    df['e-1'] = df['end'].apply(lambda x: str(int(x)-1) if not x == None else x)
    df['b+1'], df['b-1'] = df['b+1'].astype(float), df['b-1'].astype(float)
    df['e+1'], df['e-1'] = df['e+1'].astype(float), df['e-1'].astype(float)
    df['begin'], df['end'] = df['begin'].astype(float), df['end'].astype(float)

    uniprot_list = df[['ensembl_accession', 'index', 'begin', 'b+1', 'b-1', 'end', 'e+1', 'e-1']].values.tolist()

    positions = json.load(open(cwd + '/alignments/' + kind + '_alignments/alignments.json'))

    # confronto le posizioni trovate in uniprot con quelle trovate mediante allineamenti
    prova = []

    for group in positions:
        x = positions[group]
        acc_A, acc_B = x['Accessions'][0], x['Accessions'][1]

        for pos in x['Positions']:
            x2 = x['Positions'][pos]
            score = x['Positions'][pos]['Score']

            if score > 1:
                al_A, al_B = x2['Column A'], x2['Column B']
                pos_A, pos_B = float(x2['Homo position A']), float(x2['Homo position B'])
                res_A, res_B = x2['Consensus A'], x2['Consensus B']

                for line in uniprot_list:
                    accession = line[0]
                    index = line[1]
                    begin, begin_plus_one, begin_minus_one = line[2], line[3], line[4]
                    end, end_plus_one, end_minus_one = line[5], line[6], line[7]

                    if accession == acc_A:

                        if pos_A == begin or pos_A == begin_plus_one or pos_A == begin_minus_one or pos_A == end or pos_A == end_plus_one or pos_A == end_minus_one:
                            prova.append([acc_A, index, pos_A, res_A, res_B, al_A, al_B, score])
                        if pos_B == begin or pos_B == begin_plus_one or pos_B == begin_minus_one or pos_B == end or pos_B == end_plus_one or pos_B == end_minus_one:
                            prova.append([acc_A, index, pos_B, res_A, res_B, al_A, al_B, score])    

                    if accession == acc_B:
                        if pos_A == begin or pos_A == begin_plus_one or pos_A == begin_minus_one or pos_A == end or pos_A == end_plus_one or pos_A == end_minus_one:
                            prova.append([acc_B, index, pos_A, res_A, res_B, al_A, al_B, score])
                        if pos_B == begin or pos_B == begin_plus_one or pos_B == begin_minus_one or pos_B == end or pos_B == end_plus_one or pos_B == end_minus_one:
                            prova.append([acc_B, index, pos_B, res_A, res_B, al_A, al_B, score]) 

    positions_found = pd.DataFrame(prova, columns=['ensembl_accession', 'index', 'position_found', 'residue_A', 'residue_B', 'alignment_A', 'alignment_B', 'score'])

    final_df = pd.merge(df, positions_found, on=['index', 'ensembl_accession']).drop_duplicates()

    # recupero i full protein name di ensembl dal database (proteina in esame e partner)
    database = json.load(open('database.json'))
    orthos = pd.read_table(cwd + '/orthologues/' + kind + '_orthologues.csv', sep=';')
    orthos = list(zip(orthos[orthos['Species'] == 'Homo_sapiens'].values.tolist()[0][3:], orthos[orthos['Species'] == 'Homo_sapiens'].columns.tolist()[3:]))
    orthos = pd.DataFrame(orthos, columns = ['ensembl_accession', 'orthogroup'])
    final_df = pd.merge(final_df, orthos, on=['ensembl_accession', 'orthogroup'])
    orthos['full_name'] = orthos['ensembl_accession'].apply(lambda x: database.get(x)[7] if not database.get(x) == None else x)

    final_df['full_name_ensembl_partner'] = final_df['orthogroup'].apply(lambda x: database.get(
        orthos[orthos['orthogroup'] == re.split('A|B', x)[0] + {'A': 'B', 'B': 'A'}.get(x[-1])].values.tolist()[0][0])[7] if not database.get(
        orthos[orthos['orthogroup'] == re.split('A|B', x)[0] + {'A': 'B', 'B': 'A'}.get(x[-1])].values.tolist()[0][0]) == None else None)

    final_df['full_name_ensembl'] = final_df['orthogroup'].apply(lambda x: database.get(
        orthos[orthos['orthogroup'] == x].values.tolist()[0][0])[7] if not database.get(
        orthos[orthos['orthogroup'] == x].values.tolist()[0][0]) == None else None)

    # ordino la tabella e la salvo in locale
    final_df = final_df.drop(columns=['index', 'b+1', 'b-1', 'e+1', 'e-1']).drop_duplicates()[['orthogroup', 'accession', 'ensembl_accession', 'gene', 'full_protein_name', 'full_name_ensembl', 'full_name_ensembl_partner', 'type', 'category', 'description', 'begin', 'end', 'position_found', 'score', 'residue_A', 'residue_B', 'alignment_A', 'alignment_B', 'EC_number']]
    final_df = final_df[~final_df['category'].str.contains('STRUCTURAL')]
    final_df = final_df[~final_df['category'].str.contains('TOPOLOGY')]
    final_df.to_csv(cwd + '/uniprot/' + kind + '_features.csv')

features('Convergent')
features('Divergent')
features('Tandem')
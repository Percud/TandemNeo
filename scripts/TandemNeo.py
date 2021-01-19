#python3.7
#anaconda3

#conda create -n myenv3.7 python=3.7
#conda activate myenv3.7 
# conda install...

#clustalo            0.1.2
#biopython           1.78
#ipykernel           5.3.4
#ipython             7.19.0
#ipython-genutils    0.2.0
#json5               0.9.5
#jupyterlab          2.2.6
#matplotlib          3.3.2
#natsort             7.1.0
#numpy               1.19.2
#pandas              1.1.5
#pymol               2.4.1
#requests            2.25.1
#urllib3             1.26.2
#xmltodict           0.12.0

# Modules and congig import

import os, sys, gzip, json
import pandas as pd
from   multiprocessing import Pool, Manager
from   itertools       import product
from   Bio             import SeqIO

from   specieinfo      import assemblies
from   specieinfo      import specieinfo      as si
from   download        import download        as dl
from   mainisoforms    import mainisoforms    as mi
from   duplications    import duplications    as dup
from   orthology       import ortho           as ort 
from   faforalignments import faforalignments as ffal
from   database        import database        as db
from   database        import dbinfo
from   alignments      import alignments      as al
from   features        import features        as feat

cwd = os.path.dirname(os.getcwd())
sys.path.append(cwd)
from config import *

# 0. Working directories

os.mkdir(cwd + '/appris')
os.mkdir(cwd + '/alignments')
os.mkdir(cwd + '/alignments/tandem')
os.mkdir(cwd + '/alignments/divergent')
os.mkdir(cwd + '/alignments/convergent')
os.mkdir(cwd + '/blast_queries')
os.mkdir(cwd + '/duplications')
os.mkdir(cwd + '/fa')
os.mkdir(cwd + '/gtf')
os.mkdir(cwd + '/main')
os.mkdir(cwd + '/orthologues')
os.mkdir(cwd + '/features')

# 1. Species info

if not have_list:

    species = [a['name'] for a in assemblies['species']]       # collect all available species name in Ensembl database
    slist = [si(s).allinfo for s in species]                   # specieinfo class for each of those species name

    df = pd.DataFrame(slist, columns=[
                                    'Class',
                                    'Order',
                                    'Genus',
                                    'Specie',
                                    'Publications',
                                    'Taxid',
                                    'Assembly']).dropna()      # writing a dataframe with all collected informations
    df = df[df['Class'].isin(rf)]                              # filtering the dataframe for classes specified withing the config file
    df['Publications'] = df['Publications'].astype(int)
    df['Genus']        = df['Genus'].str.capitalize()          # capitalize genus value
    df = df.sort_values('Publications', ascending=False)       # sorting for the publications
    if mo:
        df = df.groupby('Order').head(mo)                      # limiting the max number of orders to collect (specified in config)
    if mg:
        df = df.groupby('Genus', sort=False).head(mg)          # same for the genres
    df = df.sort_values(['Class', 'Order'])                    # last sorting for the final dataframe
    df.to_csv(cwd + '/species/species_list.csv', index=False)

# 2. FASTAs and GTFs download

ls = pd.read_csv(cwd + '/species/species_list.csv')            # open species containing file
fa = [l[2] + '_' + l[3] for l in ls.values.tolist()]           # write a list containing all the specie names for the fasta files
gtf = rs                                                       # write a list containing all the specie names for the gtf files

pool = Pool(threads)                                           # activating the dl function in multiprocessing
pool.starmap(dl.downl, product(gtf, ['gtf']))                  # starmap multiprocessing accepts tuple with arguments [('Homo_sapiens', 'gtf'), ('Gallus_gallus', 'gtf'), etc...]
pool.starmap(dl.downl, product(fa, ['fa']))

# 3. Main isoforms and BLASTP

def mainblast(s, n, e, h, t):
    mi.tofa(s)                                                 # activating tofa function, write a principal isoforms containing .fa file
    mi.totsv(s)                                                # activating totsv function, write a principal isoforms, sorted over the genome, containing .tsv file
    mi.blast(s, n, e, h, t)                                    # performing an intraspecie blastP

blast_args = [num_threads,
              evalue,
              max_hsps,
              max_target_seqs]

pool = Pool(threads)                                           # activating the dl function in multiprocessing
args = [tuple([s] + blast_args) for s in rs]
pool.starmap(mainblast, args)                                  # starmap multiprocessing accepts tuple with arguments [('Homo_sapiens', 20, '10e-6', 1, 5), ('Gallus_gallus', 20, '10e-6', 1, 5), etc...]

# 4. Duplications

def tocsv(s, k):                                               # saving a .tsv file containing a list of duplications filtered for specie e duplication kind

    dups = dup.duplist(s, k)                                   # retrieving the accessions list for specie and duplication kind
    for l in dups:
        path = cwd + '/duplications/'
        file = path + s + '_' + k + '.tsv'
        print(*l, sep='\t', file=open(file, 'a'))              # writing duplicated pair in a tsv file

Pool(threads).starmap(tocsv, product(rs, dups))                # activating the tocsv function in multiprocessing

# 5. Orthology (COMPARA)

def forquery(ID, d):
    
    file = d + '_no_filter.csv'
    path = cwd + '/orthologues/' + file
    df   = ort.df(ID[0], True)                                 # per ogni ID fornito scarica la lista degli ortologhi in formato dataframe, se viene fornita una lista di specie viene filtrato per queste
    df['orthogroup'] = ID[1]                                   # assegna un numero corrispondende all'ortogruppo in base all'ordine di apparizione nel genoma dell'ID in esame
    df   = df.values.tolist()

    for l in df:
        print(*l, sep=';', file=open(path, 'a'))               # dopo la conversione in lista delle righe del dataframe, vengono stampati al momento i risultati in un file con tutte le informazioni non filtrate
        
def forspecie(rs, d):
        
    dups_df = ort.mergedups(rs, d)[[8, 10]]                    # unisce i dataframe contenenti i geni duplicati di tutte le specie di riferimento per ciascun tipo di duplicazione
    IDS = [l for l in dups_df.values.tolist()]                 # recupera tutti gli ID e i numeri degli ortogruppi assegnati per ogni tipo di duplicazione
    Pool(5).starmap(forquery, product(IDS, [d]))               # attiva la funzione (in multiprocessing) forquery per ogni ID presente nel dataframe nato dall'unione precedente

dups = ['convergent', 'divergent']
for d in dups:
    file = d + '_no_filter.csv'
    path = cwd + '/orthologues/' + file
    forspecie(rs, d)                                           # attiva la funzione for specie per ogni specie di riferimento e tipo di duplicazione
    ort.brh(path)                                              # brh interno nel dataframe di ortologia ottenuto, altri commenti presenti in orthology.py


# 6. Database

df    = db.orthodf('Tandem')                                   # opening orthologues dataframe for:
slist = df['Species'].values.tolist()                          # obtaining species list

manager = Manager()
data    = manager.dict()

def db_func(s):
 
    path   = cwd + '/fa/' + s + '.fa.gz'
    handle = gzip.open(path, 'rt')
    fastas = list(SeqIO.parse(handle, "fasta"))                # opening file.fa --> storing in a list

    orthos = db.aclist(s)                                      # accessions list from the orthologues tab
    f = [l for l in fastas 
         if l.id.split('.')[0] in orthos]                      # intersection between orthologues accession list and fastas

    for l in f:
        data.update(db.info(l))                                # activating database class info function
        
Pool(threads).map(db_func, 
    [s for s in slist])
        
json.dump(data.copy(), open(cwd + '/database.json', 'w'))      # dumping json database

# 7. FASTAs for alignments

for d in dups:                                                 # for each kind of duplication

    df    = db.orthodf(d)                                      # open and keep in memory the orthologues dataframe
    pairs = ffal.pairslist(d)                                  # write a list containing the orthologues dataframe column indexes corresponding to orthogroup pairs 
    suff  = db.suffixes()                                      # write a dictionary containing the species references based on accessions Ensembl coding {ENSP0: 'Homo_sapiens'}
    js    = json.load(open(cwd + '/database.json'))            # open the local database wrote in step 6
    
    for p in pairs:                                            # iterating over dataframe column indexes
        ffal.printfa(df, p, suff, js, d)                       # writing FASTA file (1.fa will contain the FASTA corresponding to columns 1A and 1B)

# 8. Alignments

for d in dups:
    
    folder = 'alignments/' + d
    ogroups = al.fanum(folder, '.fa')
    
    for o in ogroups:                                          # for each orthogroup
        
        fa   = 'alignments/' + d + '/' + str(o) + '.fa'        # FASTA file path
        faln = 'alignments/' + d + '/' + str(o) + '.fasta'     # Aligned FASTA file path
        al.clustifcov(fa, faln, th_coverage)                   # perform a clustalo if sequence lenghts coverage is higher than config threshold
        
    manager = Manager()
    logs = manager.dict()
    
    def log(ogroup, srefs, matrix):
        
        try:
            faln     = folder + '/' + str(ogroup) + '.fasta'   # Aligned FASTA file path
            falnfile = al.alignmentfile(faln)                  # aligned FASTA file
            if not ogroup in logs.keys():
                log  = al.log(falnfile, srefs, matrix)         # return a log containing alignments informations
                logs.update({ogroup: log})                     # storing those information inside a Manager dictionary for multiprocessing
        except:
            pass
    
    Pool(threads).starmap(                                     # activating log function in multiprocessing
        log, product(
            ogroups, [srefs], [matrix]))
        
    json.dump(logs.copy(), open(
        cwd + '/' + folder + '/' + d + '.json', 'w'))          # dumping those informations inside a json file

# 9. Features

for d in dups:

    ortholist = dbinfo.ids_to(                                 # IDs conversion from protein to gene 
        features_ref_specie, d, 'gene')
    converted = feat.convert_id(                               # IDs conversion from ENSEMBL to Uniprot
        'ENSEMBL_ID', 'ACC', ortholist)
    converted_IDS = [v for k,v in converted.items()            # storing in a list all uniprot IDS 
                     if not 'ENS' in v]

    manager = Manager()
    appended_data = manager.list()
    def getfeaturesparallel(k):
        data = feat.getfeatures(k)                             # using features class getfeatures function to retrieve a complete features dataset for each ID
        appended_data.append(data)
    Pool(threads).map(getfeaturesparallel, converted_IDS)      # activating getfeaturesparallel function in multiprocessing

    allfeatures = pd.concat(appended_data)                     # concatenated dataframe with features infos

    path = cwd + '/alignments/' + d + '/'
    aln = json.load(open(path + d + '.json'))                  # opening the json format alignments containing file
    alns = al.threshold_aln(
        aln, alignment_threshold)                              # set a threshold based on the alignment scores and return a dataframe
    alns = feat.add_genes_ids(alns, converted)                 # add more infos in the dataframe
    
    features = feat.intersect_alns_features(
        alns, allfeatures)                                     # return an intersection between alignments dataframe and allfeatures dataframe
    # --> feat.filter(allfeatures)
    features.to_csv(cwd + '/features/' + d + '.csv')


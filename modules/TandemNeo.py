import os, sys, gzip, json
import pandas as pd
import argparse
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

parser = argparse.ArgumentParser(prog='TandemNeo')

parser.add_argument('-t', '--threads', default=threads, help='Define available threads for multiprocessing processes')
parser.add_argument('-s', '--seed', default=seed, help='Define a seed for random sampling processes')

parser.add_argument('-l', '--have_species_list', default=have_list, help='Specify (True or False) if you already have a species list')
parser.add_argument('-ll', '--species_list', default=species_file_name, help='Species list follows the syntax: [[class, order, genus, specie, publications, taxid, assembly], [class2, order2, ...]]. Write None as seat holder if you miss any field')
parser.add_argument('-mg', '--max_genres', default=mg, help='Define max number of genres, None if you don\'t want a limit')
parser.add_argument('-mo', '--max_orders', default=mo, help='Define max number of orders, None if you don\'t want a limit')
parser.add_argument('-rf', '--ref_classes', default=rf, help='Define a reference classes list (python list format), None if you don\'t want to choose one')
parser.add_argument('-rs', '--ref_species', default=rs, help='Define a reference species list (python list format)')

parser.add_argument('-bt', '--blastp_threads', default=num_threads, help='Define blast params')
parser.add_argument('-ev', '--blastp_evalue', default=evalue, help='Define blast params')
parser.add_argument('-mh', '--blastp_max_hsps', default=max_hsps, help='Define blast params')
parser.add_argument('-mt', '--blastp_max_target_seqs', default=max_target_seqs, help='Define blast params')

parser.add_argument('-dc', '--duplication_classes', default=dups, help='Define duplication classes list (python list format)')

parser.add_argument('-m', '--matrix', default=matrix, help='Define matrix that will be used whithin alignment phases')
parser.add_argument('-ta', '--tags', default=srefs, help='Define reference specie ENSEMBL tags that will be used for extract alignment positions')
parser.add_argument('-tc', '--threshold_coverage', default=th_coverage, help='Define minimum accepted coverage between orthogroups sequence lenghts ')

parser.add_argument('-fr', '--features_ref_specie', default=features_ref_specie, help='Define reference specie for extract Uniprot features')
parser.add_argument('-pr', '--pos_range', default=pos_range, help='Define max distance between alignment residue positions and feature positions ')
parser.add_argument('-at', '--alignment_threshold', default=alignment_threshold, help='Define minumum accepted position score')

parser.add_argument('-cr', '--clust_ref_specie', default=clust_ref_specie, help='')
parser.add_argument('-ms', '--clust_min_samples', default=clust_min_samples, help='')
parser.add_argument('-ep', '--clust_eps', default=clust_eps, help='')

parser.add_argument('phase', nargs='*', action='append', help='Define a phase:')
parser.add_argument('species_list', default=None, nargs='?', help='Generate a species list file')
parser.add_argument('download', default=None, nargs='?', help='Download FASTA and GTF from Ensembl ftp servers')
parser.add_argument('main_isoforms', default=None, nargs='?', help='Looking for principal isoforms through Appris database')
parser.add_argument('duplications', default=None, nargs='?', help='Find duplicated genes along the genome')
parser.add_argument('orthology', default=None, nargs='?', help='Looking for orthology informations in Compara database')
parser.add_argument('database', default=None, nargs='?', help='Generate a database containing accessions found informations for speed up next processes')
parser.add_argument('fa_for_alignments', default=None, nargs='?', help='Generate FASTA file containing orthology columns to align')
parser.add_argument('alignments', default=None, nargs='?', help='Generate aligned files (FASTA format)')
parser.add_argument('features', default=None, nargs='?', help='Looking for features through the Uniprot database')
parser.add_argument('clustering', default=None, nargs='?', help='')

args = parser.parse_args()

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

if 'species_list' in args.phase[0] or args.phase[0] == []:

    if not args.have_list:

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
        df = df[df['Class'].isin(args.rf)]                         # filtering the dataframe for classes specified withing the config file
        df['Publications'] = df['Publications'].astype(int)
        df['Genus']        = df['Genus'].str.capitalize()          # capitalize genus value
        df = df.sort_values('Publications', ascending=False)       # sorting for the publications
        if mo:
            df = df.groupby('Order').head(args.mo)                 # limiting the max number of orders to collect (specified in config)
        if mg:
            df = df.groupby('Genus', sort=False).head(args.mg)     # same for the genres
        df = df.sort_values(['Class', 'Order'])                    # last sorting for the final dataframe
        df.to_csv(cwd + '/species/species_list.csv', index=False)

# 2. FASTAs and GTFs download

if 'download' in args.phase[0] or args.phase[0] == []:

    ls = pd.read_csv(cwd + '/species/species_list.csv')            # open species containing file
    fa = [l[2] + '_' + l[3] for l in ls.values.tolist()]           # write a list containing all the specie names for the fasta files
    gtf = args.rs                                                  # write a list containing all the specie names for the gtf files

    pool = Pool(args.threads)                                      # activating the dl function in multiprocessing
    pool.starmap(dl.downl, product(gtf, ['gtf']))                  # starmap multiprocessing accepts tuple with arguments [('Homo_sapiens', 'gtf'), ('Gallus_gallus', 'gtf'), etc...]
    pool.starmap(dl.downl, product(fa, ['fa']))

# 3. Main isoforms and BLASTP

if 'main_isoforms' in args.phase[0] or args.phase[0] == []:

    def mainblast(s, n, e, h, t):
        mi.tofa(s)                                                 # activating tofa function, write a principal isoforms containing .fa file
        mi.totsv(s)                                                # activating totsv function, write a principal isoforms, sorted over the genome, containing .tsv file
        mi.blast(s, n, e, h, t)                                    # performing an intraspecie blastP

    blast_args = [args.num_threads,
                  args.evalue,
                  args.max_hsps,
                  args.max_target_seqs]

    pool = Pool(args.threads)                                      # activating the dl function in multiprocessing
    args = [tuple([s] + blast_args) for s in args.rs]
    pool.starmap(mainblast, args)                                  # starmap multiprocessing accepts tuple with arguments [('Homo_sapiens', 20, '10e-6', 1, 5), ('Gallus_gallus', 20, '10e-6', 1, 5), etc...]

# 4. Duplications

if 'duplications' in args.phase[0] or args.phase[0] == []:

    def tocsv(s, k):                                               # saving a .tsv file containing a list of duplications filtered for specie e duplication kind

        dups = dup.duplist(s, k)                                   # retrieving the accessions list for specie and duplication kind
        for l in dups:
            path = cwd + '/duplications/'
            file = path + s + '_' + k + '.tsv'
            print(*l, sep='\t', file=open(file, 'a'))              # writing duplicated pair in a tsv file

    Pool(args.threads).starmap(tocsv, product(args.rs, args.dups)) # activating the tocsv function in multiprocessing

# 5. Orthology (COMPARA)

if 'orthology' in args.phase[0] or args.phase[0] == []:

    def forquery(ID, d):
        
        file = d + '_no_filter.csv'
        path = cwd + '/orthologues/' + file
        df   = ort.df(ID[0], True)                                 # per ogni ID fornito scarica la lista degli ortologhi in formato dataframe, se viene fornita una lista di specie viene filtrato per queste
        df['orthogroup'] = ID[1]                                   # assegna un numero corrispondende all'ortogruppo in base all'ordine di apparizione nel genoma dell'ID in esame
        df   = df.values.tolist()

        for l in df:
            print(*l, sep=';', file=open(path, 'a'))               # dopo la conversione in lista delle righe del dataframe, vengono stampati al momento i risultati in un file con tutte le informazioni non filtrate
            
    def forspecie(rs, d):
            
        dups_df = ort.mergedups(args.rs, d)[[8, 10]]               # unisce i dataframe contenenti i geni duplicati di tutte le specie di riferimento per ciascun tipo di duplicazione
        IDS = [l for l in dups_df.values.tolist()]                 # recupera tutti gli ID e i numeri degli ortogruppi assegnati per ogni tipo di duplicazione
        Pool(5).starmap(forquery, product(IDS, [d]))               # attiva la funzione (in multiprocessing) forquery per ogni ID presente nel dataframe nato dall'unione precedente

    for d in args.dups:
        file = d + '_no_filter.csv'
        path = cwd + '/orthologues/' + file
        forspecie(args.rs, d)                                      # attiva la funzione for specie per ogni specie di riferimento e tipo di duplicazione
        ort.brh(path)                                              # brh interno nel dataframe di ortologia ottenuto, altri commenti presenti in orthology.py


# 6. Database

if 'database' in args.phase[0] or args.phase[0] == []:

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
            
    Pool(args.threads).map(db_func, 
        [s for s in slist])
            
    json.dump(data.copy(), open(cwd + '/database.json', 'w'))      # dumping json database

# 7. FASTAs for alignments

if 'fa_for_alignments' in args.phase[0] or args.phase[0] == []:

    for d in args.dups:                                                 # for each kind of duplication

        df    = db.orthodf(d)                                      # open and keep in memory the orthologues dataframe
        pairs = ffal.pairslist(d)                                  # write a list containing the orthologues dataframe column indexes corresponding to orthogroup pairs 
        suff  = db.suffixes()                                      # write a dictionary containing the species references based on accessions Ensembl coding {ENSP0: 'Homo_sapiens'}
        js    = json.load(open(cwd + '/database.json'))            # open the local database wrote in step 6
        
        for p in pairs:                                            # iterating over dataframe column indexes
            ffal.printfa(df, p, suff, js, d)                       # writing FASTA file (1.fa will contain the FASTA corresponding to columns 1A and 1B)

# 8. Alignments

if 'alignments' in args.phase[0] or args.phase[0] == []:

    for d in args.dups:
        
        folder = 'alignments/' + d
        ogroups = al.fanum(folder, '.fa')
        
        for o in ogroups:                                          # for each orthogroup
            
            fa   = 'alignments/' + d + '/' + str(o) + '.fa'        # FASTA file path
            faln = 'alignments/' + d + '/' + str(o) + '.fasta'     # Aligned FASTA file path
            al.clustifcov(fa, faln, args.th_coverage)              # perform a clustalo if sequence lenghts coverage is higher than config threshold
            
        manager = Manager()
        logs = manager.dict()
        
        def log(ogroup, srefs, matrix):
            
            try:
                faln     = folder + '/' + str(ogroup) + '.fasta'   # Aligned FASTA file path
                falnfile = al.alignmentfile(faln)                  # aligned FASTA file
                if not ogroup in logs.keys():
                    log  = al.log(falnfile, args.srefs, args.matrix) # return a log containing alignments informations
                    logs.update({ogroup: log})                     # storing those information inside a Manager dictionary for multiprocessing
            except:
                pass
        
        Pool(args.threads).starmap(                                # activating log function in multiprocessing
            log, product(
                ogroups, [args.srefs], [args.matrix]))
            
        json.dump(logs.copy(), open(
            cwd + '/' + folder + '/' + d + '.json', 'w'))          # dumping those informations inside a json file

# 9. Features

if 'features' in args.phase[0] or args.phase[0] == []:

    for d in args.dups:

        ortholist = dbinfo.ids_to(                                 # IDs conversion from protein to gene 
            args.features_ref_specie, d, 'gene')
        converted = feat.convert_id(                               # IDs conversion from ENSEMBL to Uniprot
            'ENSEMBL_ID', 'ACC', ortholist)
        converted_IDS = [v for k,v in converted.items()            # storing in a list all uniprot IDS 
                         if not 'ENS' in v]

        manager = Manager()
        appended_data = manager.list()
        def getfeaturesparallel(k):
            data = feat.getfeatures(k)                             # using features class getfeatures function to retrieve a complete features dataset for each ID
            appended_data.append(data)
        Pool(args.threads).map(getfeaturesparallel, converted_IDS) # activating getfeaturesparallel function in multiprocessing

        allfeatures = pd.concat(appended_data)                     # concatenated dataframe with features infos

        path = cwd + '/alignments/' + d + '/'
        aln = json.load(open(path + d + '.json'))                  # opening the json format alignments containing file
        alns = al.threshold_aln(
            aln, args.alignment_threshold)                         # set a threshold based on the alignment scores and return a dataframe
        alns = feat.add_genes_ids(alns, converted)                 # add more infos in the dataframe
        
        features = feat.intersect_alns_features(
            alns, allfeatures)                                     # return an intersection between alignments dataframe and allfeatures dataframe
        # --> feat.filter(allfeatures)
        features.to_csv(cwd + '/features/' + d + '.csv')

# 10. 3D Clustering

if 'clustering' in args.phase[0] or args.phase[0] == []:
    pass


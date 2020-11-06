# TandemNeo

PACCHETTI DA INSTALLARE

# intanto vedi import (vedi se trovi una pip list)

CONFIG FILE

(Script 1)
- LISTA SPECIE [cwd + 'listaspecie.txt'] / se è vuota va avanti

(Script 2, 3)
- QUERIES:
    queries = ['Homo_sapiens', 'Gallus_gallus', 'Danio_rerio']

(script 3)
- PARAMETRI BLAST
    num_threads=20, evalue='10e-6', max_hsps=1, outfmt=7, max_target_seqs=5

(script 5)
- CLASSI DELLA TASSONOMIA PER FILTRARE DATAFRAME 
    ref_classes = ['Sauropsida', 'Mammalia', 'Actinopteri']
    
(script 8)
- MATRICE
    matrix = 'BLOSUM62'
    min_coverage = 60
    ensp = 'ENSP0'


# script 1 genera file input contenente lista specie "filtrata", genera albero cartelle

    # se uno ha lista di specie parte da script 2

# script 2 fa il download

    # 
    
# script 3 isoforme principali e blastp

    # 
    
# script 4 duplicazioni

    #

# script 5 ortologia compara

    #
    
# script 6 database

    #
    
# script 7 FASTAS allineamento

    #
    
# script 8 Allineamenti

    # 
    
# script 9 Features (UNIPROT)
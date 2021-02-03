from tandemneo.clustering import *
import sys
from tandemneo.config import *

organism = sys.argv[1] #'Homo_sapiens'
orthotab = sys.argv[2] #'https://raw.githubusercontent.com/Percud/TandemNeo/master/orthologues/Tandem_orthologues.csv'
cwd=os.getcwd() # define cwd

try:
    os.chdir(organism) # change directory if exists i.e. Homo_sapiens
except:
    os.mkdir(organism) # make directory if not exists
    os.chdir(organism)
    
ortho=pd.read_csv(cwd + '/' + orthotab,';') # read ortologue table
ens_acc = ortho.loc[ortho['Species'] == organism].iloc[:,3:].dropna(axis=1).values.tolist()[0] # get all the ensembl accessions
fromto=(convert_ac(ens_acc,'ENSEMBL_ID','ACC') # convert gene_id ensembl in uniprot ac
        .drop_duplicates('From') # get only the best result for uniprot ac
        .rename(columns = {'From':'ensembl_ac','To':'uniprot_ac'})) #rename from/to in ensembl/uniprot
accession=fromto.uniprot_ac.tolist() # uniprot ac into list

def get_fa_pdb_match(i):
    try:
        urllib.request.urlretrieve(fr'https://swissmodel.expasy.org/repository/uniprot/{i}.pdb', fr'{i}.pdb') # download pdb from SM repository
        urllib.request.urlretrieve(fr'https://www.uniprot.org/uniprot/{i}.fasta', fr'{i}.fasta') #download fasta from uniprot
        pdb2fasta(i+'.pdb',i+'.fa') # convert pdb in fasta format
        return (match_fasta_position(i+'.fasta',i+'.fa')) # match positions between uniprot fasta and pdbfasta
    except:
        pass
    
dfs=Pool(threads).imap(get_fa_pdb_match,set(accession)) # multiprocessing for each uniprot ac

try:
  df=pd.concat(dfs).reset_index() # concatenate dataframes of multiprocessing
  df[['uniprot_ac','fasta']] = pd.DataFrame(df['query'].str.split('.').values.tolist(), index=df.index) # define new columns
  df.drop(['fasta','query'], axis=1, inplace=True) # drop columns 
  s=df.merge(fromto, on='uniprot_ac').rename(columns = {'sequence': 'chain_id', # rename all columns 
                                                        'query_res': 'ensembl_aa',
                                                        'query_num':'ensembl_num',
                                                        'hit_res':'pdb_aa',
                                                        'hit_num':'pdb_num',
                                                       })
  s.to_csv(organism+'.csv', index=False) # save matching position to csv file
except:
    pass
os.chdir(cwd) # back to cwd

import pandas as pd, tarfile, urllib, wget, os, glob, requests, json, random, numpy as np, re
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast import NCBIXML 
from Bio import SeqIO, SearchIO
from Bio.PDB import PDBParser as PDB, PDBIO
from sklearn.cluster import DBSCAN
from multiprocessing import Pool
from biopandas.pdb import PandasPdb

def match_fasta_position(query,subject,num=None):
    """return dataframe of position matching in subject and query fasta, or a dictionary given a list of number"""
    df=[]
    blastp(query=query, subject=subject, out=query+'_'+subject+'.xml', outfmt=5, max_hsps=1)()
    xml = SearchIO.read(query+'_'+subject+'.xml', "blast-xml")
    for n in range(len(xml)):
        x=xml[n][0]
        hit_gap = np.array([index for index, value in enumerate(x.hit) if value == '-'])
        q_gap = np.array([index for index, value in enumerate(x.query) if value == '-'])
        hit_num=list(np.arange(x.hit_start+1,x.hit_end+1))
        q_num=list(np.arange(x.query_start+1,x.query_end+1))
        for i in hit_gap:
            hit_num.insert(i,np.nan)
        for i in q_gap:
            q_num.insert(i,np.nan)
        df.append(list(zip(len(x.query)*[query],len(x.query)*xml[n].description,x.query,q_num,x.hit,hit_num)))
    df = pd.DataFrame([j for i in df for j in i],columns=['query','sequence','query_res','query_num','hit_res','hit_num'])
    df.set_index(['query'],inplace=True)
    os.remove(query+'_'+subject+'.xml')
    if num==None:
        return df
    else:
        return df[df.query_num.isin(list(map(int,num)))].to_dict('records')

def convert_ac(list,From,to):
    """convert accession list"""
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
    'from': From,
    'to': to,
    'format': 'tab',
    'query': ' '.join(list)
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    return pd.read_csv(urllib.request.urlopen(req), sep='\t')

def pdb2fasta(pdb,fasta):
    """convert pdb file to fasta format for each chain"""
    with open(fasta,'w') as outfa:
        for record in SeqIO.parse(pdb, 'pdb-atom'):
            if record.annotations['start'] < 1:
                print('>' +record.annotations['chain']+'\n'+record.seq[1:],file = outfa)
            else:
                print('>' +record.annotations['chain']+'\n'+'X'*(record.annotations['start']-1)+record.seq,file = outfa)

def getfeatures(uniprot_id):
    """return the Uniprot complete dataset for a given Uniprot ID"""
    try:
        r = requests.get("https://www.ebi.ac.uk/proteins/api/proteins?size=-1&accession=" + uniprot_id, 
                         headers={ "Accept" : "application/json"})
        data = pd.json_normalize(r.json())
        return data
    except:
        return str(uniprot_id) + ' not_found'

def threshold_aln(aln,threshold):
    '''Return dataframe with position > threshold, given a alignment.json'''
    dfs=[]
    for k in aln:
        scores = np.array([n['Score difference'] for n in aln[k]['Positions'].values()])
        positions = np.array([n[0] for n in aln[k]['Positions'].items()])
        tags = [i[0] for i in aln[k]['Positions']['1'].items() if type(i[1])==dict]
        for tag in zip(tags, list(reversed(tags))):
            aca = ''.join([k[0] for k in aln[k]['Positions']['1'][tag[0]].items() if type(k[1])==dict])
            posa = np.array([n[1][tag[0]][aca]['Positions'] 
                             if aca in n[1][tag[0]].keys() and aca!=''
                             else aca
                             for n in aln[k]['Positions'].items() 
                            ])
            resa = np.array([n[1][tag[0]][aca]['Residue'] 
                             if aca in n[1][tag[0]].keys() and aca!=''
                             else aca
                             for n in aln[k]['Positions'].items() 
                            ])
            cola = np.array([n[1][tag[0]]['Column'] for n in aln[k]['Positions'].items()])
            acb = ''.join([k[0] for k in aln[k]['Positions']['1'][tag[1]].items() if type(k[1])==dict])
            posb = np.array([n[1][tag[1]][acb]['Positions'] 
                             if acb in n[1][tag[1]].keys() and acb!=''
                             else acb 
                             for n in aln[k]['Positions'].items()
                            ])
            resb = np.array([n[1][tag[1]][acb]['Residue'] 
                             if acb in n[1][tag[1]].keys() and acb!=''
                             else acb
                             for n in aln[k]['Positions'].items()
                            ])
            colb = np.array([n[1][tag[1]]['Column'] for n in aln[k]['Positions'].items()])
            cons = np.array([n['Consensus'] for n in aln[k]['Positions'].values()])
            consa = np.array([n[tag[0]]['Consensus'] for n in aln[k]['Positions'].values()])
            consb = np.array([n[tag[1]]['Consensus'] for n in aln[k]['Positions'].values()])
            cut_off = scores>=int(threshold)
            dfs.append(pd.DataFrame([[k]*len(posa[cut_off]),[tag[0]]*len(posa[cut_off]),cons[cut_off],positions[cut_off],
                                     [aca]*len(posa[cut_off]),posa[cut_off],cola[cut_off],resa[cut_off],consa[cut_off],
                                     [acb]*len(posa[cut_off]),posb[cut_off],colb[cut_off],resb[cut_off],consb[cut_off],
                                     scores[cut_off],[sum(scores[cut_off])]*len(posa[cut_off]),[sum(scores)]*len(posa[cut_off])],
                                    index=['pair', 'ab', 'consensus', 'position',
                                           'ensembl_ac', 'ensembl_num', 'col', 'res', 'consensus_ac',
                                           'ensembl_ac_partner', 'ensembl_num_partner', 'col_partner', 'res_partner', 'consensus_ac_partner',
                                           'scores','total_scores>1', 'total_scores']).T.fillna(''))
    return pd.concat(dfs, ignore_index=True)

def get_coord_ca(directory):
    '''return dataframe with coordinates of only carbon alpha, given a pdb directory'''
    repo=[]
    for ac in glob.glob(directory+'/*pdb'):
        d = PandasPdb().read_pdb(ac).df['ATOM']
        d = d[d['atom_name']=='CA'][['residue_number', 'chain_id', 'x_coord', 'y_coord', 'z_coord']]
        d['uniprot_ac']=[os.path.basename(ac.split('.')[0])]*len(d)
        repo.append(d)
    return pd.concat(repo)

def selection_pymol(pdb_dir,tabella,entry,cluster=None):
    '''Show on pymol -R the selection given Uniprot ID and table of clustering'''
    try:
        tabella=pd.read_csv(tabella,sep=';')
        name = (tabella[tabella['Entry']==entry]['Gene names'].tolist()[0].split(' ')[0])
        import xmlrpc.client as xmlrpclib
        try:
            cmd = xmlrpclib.ServerProxy('http://localhost:9123')
            cmd.reinitialize()
            cmd.load(pdb_dir+'/'+entry+'.pdb')
        except:
            print('Pymol not connected')
        if cluster == None:
            cluster=set(tabella[tabella['Entry']==entry]['Cluster'].values)

        for c in cluster:
            selection=' or '.join(['(resi '+str(int(x[0]))+' in chain '+str(x[1])+')' for x in tabella[(tabella['Cluster']==int(c))&(tabella['Entry']==entry)][['pdb_num','chain_id']].values.tolist()])+' in '+entry
            try:
                cmd.create(entry+'_'+name+'_'+str(c), selection)
                cmd.center(selection)
                cmd.set_color('color'+str(c), random.sample(range(0, 255), 3))
                cmd.color('color'+str(c), '(name C*) and '+entry+'_'+name+'_'+str(c))
                cmd.show_as('spheres', entry+'_'+name+'_'+str(c))
                cmd.center()
                cmd.zoom()
#                 cmd.png(pdb_dir+'/'+entry+'_'+name+'_'+str(cluster)+'.png', 1000, 1000, 300, 1)
        #         print('load '+entry+'.pdb\nselect '+entry+'_'+name+'_'+str(cluster)+', '+selection+'\ncolor yellow, (name C*) in '+entry+'_'+name+'_'+str(cluster)+'\nshow spheres, '+entry+'_'+name+'_'+str(cluster))
            except:
                print('Invalid Selection')
    except:
        print(None)

def distance(x,y,z):
    '''Return max distance diven arrays in xyz of points'''
    return (((max(x)-min(x))**2+(max(y)-min(y))**2+(max(z)-min(z))**2)**0.5)

def clustering(ca, matching, aln_json, threshold, min_samples, eps):
    '''Return clusters in structure given, threshold of DBSCAN and dataframes for coordinates and alignment scores'''
    clusters={}
    pdbca,matching = pd.read_csv(ca), pd.read_csv(matching)
    sc = pd.merge(matching,threshold_aln(aln_json,threshold).astype({'ensembl_num':'float64'}), on=['ensembl_ac','ensembl_num']).dropna()
    sc_pdb = pdbca.merge(sc,left_on = ['uniprot_ac','chain_id','residue_number'], 
                       right_on = ['uniprot_ac','chain_id','pdb_num'],
                       how = 'inner')
    for ac in set(sc_pdb['uniprot_ac'].to_list()):
        points = sc_pdb[sc_pdb['uniprot_ac']==ac][['pdb_num','chain_id','scores']].to_numpy()
        coord = sc_pdb[sc_pdb['uniprot_ac']==ac][['x_coord','y_coord','z_coord']].to_numpy()
        x = sc_pdb[sc_pdb['uniprot_ac']==ac][['x_coord']].to_numpy()
        y = sc_pdb[sc_pdb['uniprot_ac']==ac][['y_coord']].to_numpy()       
        z = sc_pdb[sc_pdb['uniprot_ac']==ac][['z_coord']].to_numpy()

        ca=pdbca[pdbca['uniprot_ac']==ac]
        max_dist = distance(ca.x_coord.tolist(),ca.y_coord.tolist(),ca.z_coord.tolist())

        if len(coord) >0:
            pred = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(coord)
            for clu in list(dict.fromkeys(pred)):
                if clu > -1:
                    resn=points[pred==clu]
                    x_clu, y_clu, z_clu = x[pred==clu], y[pred==clu], z[pred==clu] 
                    volume_clu = (max(x_clu)-min(x_clu))*(max(y_clu)-min(y_clu))*(max(z_clu)-min(z_clu))
                    max_dist_clu = distance(x_clu,y_clu,z_clu)[0]
                    value = (list(pred).count(clu)/len(list(pred))*(list(pred).count(clu)/max_dist_clu)/max_dist_clu)*1000
                    clusters.setdefault(ac,{}).update({clu:{
                                                      'Cluster Residues':list(pred).count(clu),
                                                      'Residues':resn,
                                                      'xyz_dimension':[max(x_clu)-min(x_clu), max(y_clu)-min(y_clu), max(z_clu)-min(z_clu)],
                                                      'Cluster Space':volume_clu,
                                                      'Max Distance Percentage (%)':round(max_dist_clu/max_dist*100, 3),
                                                      'Cluster Density':list(pred).count(clu)/max_dist_clu,
                                                      'Cluster Fraction (%)':round(list(pred).count(clu)/len(list(pred))*100,3),
                                                      'Max Distance (A)':max_dist_clu,
                                                      'Value':round(list(pred).count(clu)/(max_dist_clu/max_dist*100),3),  
                                                        }})  
    ddf=pd.DataFrame.from_dict({(i,j): clusters[i][j] 
                           for i in clusters.keys() 
                           for j in clusters[i].keys()},
                       orient='index')
    human_proteome = pd.read_csv('https://www.uniprot.org/uniprot/?query=organism:9606&columns=id,entry%20name,reviewed,protein%20names,genes,organism,length,ec&format=tab', sep = '\t')
    ddf.reset_index(inplace=True)
    ddf.rename(columns={'level_0': 'uniprot_ac','level_1':'Cluster'}, inplace=True)
    ddf=ddf.explode('Residues')
    ddf[['pdb_num','chain_id','scores']] = pd.DataFrame(ddf['Residues'].tolist(),index=ddf.index)
    ddf=ddf.drop('Residues', axis=1)
    finale=pd.merge(ddf,
                    sc_pdb,on=['uniprot_ac','chain_id','pdb_num','scores'])
    clustering=finale.rename(columns={'uniprot_ac': 'Entry'})
    tabella = pd.merge(human_proteome[['Entry','Gene names','Protein names','EC number']],
                       clustering, 
                       on='Entry')
    output = 'clustering_1_'+str(min_samples)+'_'+str(eps) + os.path.basename(aln_json).split('.')[0] + '.csv'
    tabella.to_csv(output, sep=';',index=False) ###SAVE CSV
    return tabella
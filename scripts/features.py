import requests
import pandas as pd
import urllib.request
import os, sys, json

cwd = os.path.dirname(os.getcwd())
sys.path.append(cwd)
from config import *

class features:
    
    def convert_id(from_, _to, ID_list):
      """ENSEMBL to UNIPROT ID conversion 
      for a given IDs list"""

        # retrieve/ID Uniprot mapping requires 
        # a space separated accessions list
        ID_list_joined = ' '.join(map(str, ID_list))

        # setting up accessions for conversion 
        # phase (ensembl --> uniprot)
        params = {'from': from_, 'to': _to, 
                  'format': 'tab', 'query': ID_list_joined}
        url = 'https://www.uniprot.org/uploadlists/'
        req = urllib.request.Request(url, 
                                     urllib.parse.urlencode(params)
                                     .encode('utf-8'))
        reqopen = urllib.request.urlopen(req)

        df = pd.read_table(reqopen)
        df = df.drop_duplicates('From')
        df = df.values.tolist()

        converted = {
            **{l[0]: l[1] for l in df}, 
            **{l[1]: l[0] for l in df}
            }

        return converted

        # Ensembl GENE_ID = ENSEMBL_ID | Uniprot PROT_ID = ACC
        # converted = features.convert_id('ENSEMBL_ID', 'ACC', homo)

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

    def getuniprotinfo(uniprot_id):
      """return the Uniprot complete dataset 
      (json format) for a given Uniprot ID"""
        
        try:
            r = requests.get("https://www.ebi.ac.uk/proteins/api/proteins?size=-1&accession=" 
                             + uniprot_id, headers={ "Accept" : "application/json"})
            data = pd.json_normalize(r.json())
            return data
        except:
            return str(uniprot_id) + ' not_found'

        # features.getuniprotinfo('Q86YH2')
        
    def getfeatures(uniprot_id):
      """return the normalized uniprot 
      features column as dataframe"""
        
        data = features.getuniprotinfo(uniprot_id)
            
        cols = data.columns.tolist()

        feat = pd.DataFrame(
            data['features'][0]) if 'features' in cols else pd.DataFrame([])

        if not feat.empty:

            gene = data[
                'gene'][0][0
                          ]['name']['value'] if 'gene' in cols else None

            name = data[
                'protein.recommendedName.fullName.value'
            ][0] if 'protein.recommendedName.fullName.value' in cols else None

            ECN = data[
                'protein.recommendedName.ecNumber'
            ] if 'protein.recommendedName.ecNumber' in cols else None

            feat['protein id'] = uniprot_id
            feat['gene id'] = gene
            feat['product'] = name
            feat['EC number'] = ECN

            return feat

        else:
            return None

        # features.getfeatures('Q86YH2')
        
    def filter(df, to_remove):
      """filter the dataframe based 
      on a given features list"""
        
        feats = ['SIGNAL',
                 'PROPEP',
                 'CHAIN',
                 'DOMAIN',
                 'SITE',
                 'CARBOHYD',
                 'DISULFID',
                 'VARIANT',
                 'CONFLICT',
                 'TURN',
                 'STRAND',
                 'HELIX',
                 'LIPID',
                 'VAR_SEQ',
                 'REGION',
                 'ACT_SITE',
                 'METAL',
                 'MOD_RES',
                 'TOPO_DOM',
                 'TRANSMEM',
                 'CA_BIND',
                 'MUTAGEN',
                 'COMPBIAS',
                 'BINDING',
                 'MOTIF',
                 'ZN_FING',
                 'PEPTIDE',
                 'INIT_MET',
                 'REPEAT',
                 'NP_BIND',
                 'CROSSLNK',
                 'TRANSIT']

        #for l in to_remove:
        #    feats.remove(l)

        df = df[df['type'].isin(feats)]

        return df
        # features.filter(allfeatures)

    def add_genes_ids(alns, converted):
      """using phase 6 database to 
      convert protein ID to gene ID, 
      and add partner informations"""
	    
	    js = json.load(open(cwd + '/database.json'))

	    alns['ensembl_ac'] = alns[
	        'ensembl_ac_protein'].apply(
	        lambda x: js.get(x)['gene id'] 
	        if str(x) != '' else x)

	    alns['ensembl_ac_partner'] = alns[
	        'ensembl_ac_protein_partner'].apply(
	        lambda x: js.get(x)['gene id'] 
	        if str(x) != '' else x)

	    alns['uniprot_id'] = alns[
	        'ensembl_ac'].apply(
	        lambda x: converted.get(x))

	    alns['uniprot_id_partner'] = alns[
	        'ensembl_ac_partner'].apply(
	        lambda x: converted.get(x))

	    return alns

    def rearrange(df):
      """simply a dataframe rearrangment"""
        
        df['pair'] = df['pair'].astype(int)
        df = df.sort_values('pair')
        df = df[[
                'pair',
                'ab',
                'ensembl_ac_protein',
                'ensembl_ac_protein_partner',
                'ensembl_ac',
                'ensembl_ac_partner',
                'uniprot_id',
                'uniprot_id_partner',
                'gene id',
                'product',
                'begin',
                'end',
                'position',
                'scores',
                'res',
                'res_partner',
                'ensembl_num',
                'ensembl_num_partner',
                'type',
                'category',
                'description',
                'alternativeSequence',
                'EC number',
                'total_scores>1',
                'total_scores',
                'col',
                'col_partner']]
        
        return df

    def intersect_alns_features(alns, allfeatures):

    	found_indexes  = []
    	found_position = []
    	allfeatures = allfeatures.reset_index(drop=True)
    	idxs = allfeatures.index.tolist()

    	groups = list(alns[['uniprot_id', 'ensembl_num']].groupby('uniprot_id'))
    	for x in range(len(groups)):

    		al_ac  = groups[x][1]['uniprot_id'].values[0]
    		al_pos = groups[x][1]['ensembl_num'].values

    		feats = allfeatures[allfeatures['protein id'] == al_ac]

    		for i in idxs:
    			if i in feats.index.tolist():
    				for p in al_pos:
    					try:
    						ifeatures = feats.loc[i]
    						begin = int(ifeatures['begin'])
    						end = int(ifeatures['end'])
    						if not p == '':
    							if begin-pos_range <= int(p) <= end+pos_range:
    								found_indexes.append(i)
    								found_position.append(p)
    					except:
    						pass

    	df = allfeatures[allfeatures.index.isin(found_indexes)].reset_index()
    	posdict = pd.DataFrame(list(zip(found_indexes, found_position))).groupby(0)[1].apply(list).to_dict()
    	df['found'] = df['index'].apply(lambda x: posdict.get(x))
    	df = df.explode('found')
    	df = pd.merge(alns, df, left_on=['uniprot_id', 'ensembl_num'], right_on=['protein id', 'found'])
    	df = features.rearrange(df)
    	return df


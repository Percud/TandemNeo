import os, sys
from Bio import AlignIO, Align, SeqIO
from Bio.Align import substitution_matrices, MultipleSeqAlignment as MSA, AlignInfo
from random import sample
import statistics as st
from itertools import combinations
from Bio.Align.Applications import ClustalOmegaCommandline as clustalo
import numpy as np
import pandas as pd

cwd = os.path.dirname(os.getcwd())
sys.path.append(cwd)

class alignments:
        
    def fastafile(fa):

        if fa[0] == '/':
            fa = fa
        else:
            fa = cwd + '/' + fa

        file = list(SeqIO.parse(fa, 'fasta'))

        return file
    
        # al.fastafile('2.fa')
        
    def alignmentfile(faln):
            
        if faln[0] == '/':
            faln = faln
        else:
            faln = cwd + '/' + faln

        alignment = AlignIO.read(faln, 'fasta')

        return alignment
    
    def fanum(fold, pattern):
        
        f = os.listdir(cwd + '/' + fold)
        f = [l for l in f if pattern in str(l)]
        n = sorted([int(l.split('.')[0]) for l in f])
        return n

        # alignments.fanum('alignments/convergent', '.fa')
    
    def conservation(msa, matrix): 

        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load(matrix)

        numbers = []
        columns = []
        scores = []

        for n in range(msa.get_alignment_length()):

            columns.append((msa[:,n]))
            numbers.append(n+1)

        for c in columns:

            c = list(c)
            if 'X' in c:
                c.remove('X')

            pairs = list(combinations(c,2))
            score = []
            try: 
                for p in pairs:  
                    if not '-' in p:
                        score.append(
                            aligner.score(p[0],p[1]))
                scores.append(sum(score)/len(pairs))
            except:
                pass

        return (list(zip(numbers,columns,scores)))
    
    def lenght(fafile):

        lena = [len(f.seq) for f in fafile 
                if 'A' in f.description[-1]]
        lenb = [len(f.seq) for f in fafile 
                if 'B' in f.description[-1]]
        
        return [lena, lenb]
        
    def coverage(fafile): ######### coverage ########
        
        lena = alignments.lenght(fafile)[0]
        lenb = alignments.lenght(fafile)[1]
        
        if not lena == [] and not lenb == []:
        
            stma = st.mean(lena)
            stmb = st.mean(lenb)

            coverage = round(min(stma,stmb)/
                             max(stma,stmb)
                             *100,2)
        
        else:
            
            coverage = 0

        return coverage
    
    def clustifcov(fa, faln, th_coverage):

        fafile = alignments.fastafile(fa)
        
        if fa[0] == '/':
            fa = fa
        else:
            fa = cwd + '/' + fa
        
        if faln[0] == '/':
            faln = faln
        else:
            faln = cwd + '/' + faln

        if alignments.coverage(fafile) > th_coverage and not os.path.exists(faln):
            clustalo(infile = fa, outfile = faln)()
            
    def log(falnfile, srefs, matrix):

        tags = {record.description.split('\t')[-1] 
                for record in falnfile}

        alns=[[s for s in falnfile 
               if tag in s.description.split('\t')[-1]] 
              for tag in tags]

        samples = [MSA(sample(s, min([(len(s)) 
                                      for s in alns]))) 
                   for s in alns]

        total = MSA(j for i in samples 
                    for j in i)

        diff = list(zip(*[
            alignments.conservation(total, matrix)]+[
            alignments.conservation(aln, matrix) 
            for aln in  samples]))

        scores = [min([x[2] for x in d[1:]])-d[0][2] 
                  for d in diff]

        consensus = [
            AlignInfo.SummaryInfo(falnfile).dumb_consensus()]+[
            AlignInfo.SummaryInfo(aln).dumb_consensus() 
            for aln in samples]

        positions = {p+1:
                     {tag:
                      {record.name:
                       {'Positions':p-list(record.seq[:p]).count('-')+1 
                        if record.seq[p] != '-' else None,
                        'Residue':record.seq[p]
                       }
                      for record in falnfile 
                       for ref in srefs 
                       if ref in record.name 
                       if tag in record.description.split('\t')[-1]}
                     for tag in tags}
                    for p in range(len(falnfile[0]))}

        for pos in positions.keys():
            for tag in zip(range(len(tags)),tags):

                p=pos-1
                
                positions[pos].update({
                    'Score difference': scores[p], 
                    'Column': diff[p][0][1],
                    'Consensus':consensus[0][p]})
                
                positions[pos][tag[1]].update({
                    'Score': diff[p][tag[0]+1][2],
                    'Column':diff[p][tag[0]+1][1],
                    'Consensus':consensus[tag[0]+1][p]})

                log={
                    'Coverage':0,
                    'Sequences':len(falnfile),
                    'Length':len(falnfile[0]),
                    'Positions':positions}

        return log

        #json.dump(log,open('prova.json','w'))
        
    def threshold_aln(aln,threshold):
        '''Return dataframe with position > threshold, given a alignment.json'''

        dfs=[]
        for k in aln:

            scores = np.array([
                n['Score difference'] 
                for n in aln[k]['Positions'].values()
                ])

            positions = np.array([
                n[0] for n in aln[k]['Positions'].items()
                ])

            tags = [
                i[0] for i in aln[k]['Positions']['1'].items() 
                if type(i[1])==dict
                ]

            for tag in zip(tags, list(reversed(tags))):

                aca = ''.join([
                    k[0] for k in aln[k]['Positions']['1'][tag[0]].items() 
                    if type(k[1])==dict
                    ])

                posa = np.array([
                    n[1][tag[0]][aca]['Positions'] 
                    if aca in n[1][tag[0]].keys() and aca!='' else aca
                    for n in aln[k]['Positions'].items()
                    ])

                resa = np.array([
                    n[1][tag[0]][aca]['Residue'] 
                    if aca in n[1][tag[0]].keys() and aca!='' else aca
                    for n in aln[k]['Positions'].items() 
                    ])

                cola = np.array([
                    n[1][tag[0]]['Column'] for n in aln[k]['Positions'].items()
                    ])

                acb = ''.join([
                    k[0] for k in aln[k]['Positions']['1'][tag[1]].items() 
                    if type(k[1])==dict
                    ])

                posb = np.array([
                    n[1][tag[1]][acb]['Positions'] 
                    if acb in n[1][tag[1]].keys() and acb!='' else acb 
                    for n in aln[k]['Positions'].items()
                    ])

                resb = np.array([
                    n[1][tag[1]][acb]['Residue'] 
                    if acb in n[1][tag[1]].keys() and acb!='' else acb
                    for n in aln[k]['Positions'].items()
                    ])

                colb = np.array([
                    n[1][tag[1]]['Column'] for n in aln[k]['Positions'].items()
                    ])

                cons = np.array([
                    n['Consensus'] for n in aln[k]['Positions'].values()
                    ])

                consa = np.array([
                    n[tag[0]]['Consensus'] for n in aln[k]['Positions'].values()
                    ])

                consb = np.array([
                    n[tag[1]]['Consensus'] for n in aln[k]['Positions'].values()
                    ])

                cut_off = scores>=int(threshold)

                dfs.append(
                    pd.DataFrame([

                        [k]*len(posa[cut_off]),
                        [tag[0]]*len(posa[cut_off]),
                        cons[cut_off],
                        positions[cut_off],
                        [aca]*len(posa[cut_off]),
                        posa[cut_off],
                        cola[cut_off],
                        resa[cut_off],
                        consa[cut_off],
                        [acb]*len(posa[cut_off]),
                        posb[cut_off],
                        colb[cut_off],
                        resb[cut_off],
                        consb[cut_off],
                        scores[cut_off],
                        [sum(scores[cut_off])]*len(posa[cut_off]),
                        [sum(scores)]*len(posa[cut_off])

                    ],                 

                        index=['pair', 'ab', 'consensus', 'position',
                               'ensembl_ac', 'ensembl_num', 'col', 'res', 'consensus_ac',
                               'ensembl_ac_partner', 'ensembl_num_partner', 'col_partner', 
                               'res_partner', 'consensus_ac_partner',
                               'scores','total_scores>1', 'total_scores']

                    ).T.fillna(''))

        return pd.concat(dfs, ignore_index=True)

        
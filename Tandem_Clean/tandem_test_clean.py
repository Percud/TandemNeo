# /usr/bin/env python

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from ftplib import FTP
import re, itertools, pandas as pd, os, gzip

cwd = os.getcwd()
dirs, dirt, dirtc, dirtm, dirtmf, dirg = os.path.join(cwd + '/main'), os.path.join(cwd + '/tandem'), os.path.join(cwd + '/tandem_clust'), os.path.join(cwd + '/tandem_vs_main'), os.path.join(cwd + '/tandem_vs_main/fasta'), os.path.join(cwd + '/gtf_faa')
if not os.path.exists(dirs) and not os.path.exists(dirt) and not os.path.exists(dirtc) and not os.path.exists(dirtm) and not os.path.exists(dirtmf) and not os.path.exists(dirg):
    os.mkdir(dirs), os.mkdir(dirt), os.mkdir(dirtc), os.mkdir(dirtm), os.mkdir(dirtmf), os.mkdir(dirg)

ftp = FTP('ftp.ncbi.nih.gov')
ftp.login()

print('########### DOWNLOADING #############')
for line in open('species_list.txt', 'r'):
    if not os.path.exists(cwd + '/assembly_summary_refseq.txt'):
        ftp.cwd('/genomes/refseq/')
        ftp.retrbinary("RETR " + 'assembly_summary_refseq.txt', open('assembly_summary_refseq.txt', 'wb').write)
    if not '#' in line:
        taxid, species = re.split('\t', line)[0], re.split('\t', line)[1].strip()
        specie = (re.split('\s', species)[0] + '_' + re.split('\s', species)[1])
        for line2 in open('assembly_summary_refseq.txt'):
            if (taxid + '\t' + species) in line2 and not os.path.exists(cwd + '/gtf_faa/' + specie + '.gtf.gz'):
                print(specie)
                source = str((re.compile('/genomes/all/\w{3}/\d{3}/\d{3}/\d{3}/GCF_\d+.\d+_.+\t')).findall(line2))[2:-6]
                version = str((re.compile('/genomes/all/\w{3}/\d{3}/\d{3}/\d{3}/(GCF_\d+.\d+_.+\t)')).findall(line2))[2:-6]
                ftp.cwd(source)
                ftp.retrbinary("RETR " + version + '_genomic.gtf.gz', open(cwd + '/gtf_faa/' + specie + '.gtf.gz', 'wb').write)
                ftp.retrbinary("RETR " + version + '_protein.faa.gz', open(cwd + '/gtf_faa/' + specie + '.faa.gz', 'wb').write)
ftp.quit()

print('\n########### PROCESSING ##############')
for line in open('species_list.txt', 'r'):
    if not '#' in line:
        specie = (re.split('\s', re.split('\t', line)[1].strip())[0] + '_' + re.split('\s', re.split('\t', line)[1].strip())[1])
        string = {}
        with gzip.open(cwd + '/gtf_faa/' + specie + '.gtf.gz', 'rt') as gtf:
            for lineg in gtf:
                if 'NC_' in lineg:
                    if 'protein_coding' in lineg:
                        gene, chromosome, strand = (re.split('[\t"]', lineg)[9]).replace(' ', ''), re.split('[\t]', lineg)[0], re.split('[\t]', lineg)[6]
                        names = (chromosome + '\t' + gene + '\t' + strand)
                    elif 'start_codon' in lineg:
                        if 'protein_id' in lineg:
                            ac, start, product = re.split('"', (re.search('; protein_id "(.*)', lineg)).group(1))[0], (re.split('\t', lineg))[3], (re.search('; product "(.*)"; protein', lineg)).group(1)
                    elif 'stop_codon' in lineg:
                        stop = (re.split('\t', lineg))[4]
                        string.setdefault(names, {}).update({ac: [start, stop, product]})

        if not os.path.exists(cwd + '/main/' + specie + '.tsv'):
            print('Chromosome\tGene\tStrand\tCDS\tStart\tStop\tProduct', file=open(cwd + '/main/' + specie + '.tsv', "w"))
            for e, i in string.items():
                first, out = list(i.items())[0], dict(itertools.islice(i.items(), 1))
                print(e + '\t' + '\t'.join("{}\t{}".format(k, '\t'.join(v)) for k, v in out.items()), file=open(cwd + '/main/' + specie + '.tsv', "a"))
            tab = pd.read_csv(cwd + '/main/' + specie + '.tsv', sep='\t', header=0)
            tab = tab.sort_values(['Chromosome', 'Start'], ascending=(True, True))
            tab.to_csv(cwd + '/main/' + specie + '.tsv', sep='\t', index=False)

        header, sequence = [], []
        fasta_sequences = SeqIO.parse(gzip.open(cwd + '/gtf_faa/' + specie + '.faa.gz', 'rt'), 'fasta')
        for fasta in fasta_sequences:
            header.append(fasta.id), sequence.append(str(fasta.seq))
        fa = dict(zip(header, sequence))
        if not os.path.exists(cwd + '/main/' + specie + '_main.faa'):
            print(specie + ': ')
            processed = 0
            for e, i in string.items():
                protid, gene = (list(i.items())[0])[0], re.split('\t', e)[1]
                if protid in header:
                    processed += 1
                    print('>' + protid + ' ' + gene + ' [' + specie + ']' + '\n' + fa.get(protid), file=open(cwd + '/main/' + specie + '_main.faa', 'a'))
            print('- Parsing: ' + str(processed) + ' main isoforms found')

        if 'query' in line and not os.path.exists(cwd + '/tandem/' + specie + '_tandem' + '.faa'):
            ac_seen, ac2_seen, faa_seen, processed = [], [], {}, 0
            print('Gene 1\tAc 1\tGene 2\tAc 2\tid %\tsom %\tE-value', file=open(cwd + '/tandem_clust/' + specie + '_tandemclust' + '.tsv', 'w'))
            with open(cwd + '/main/' + specie + '.tsv') as tab:
                lines = [line.rstrip() for line in tab]
                res = [[lines[i], lines[i + 1]] for i in range(len(lines) - 1)]
                for i in res:
                    line1, line2 = i[0], i[1],
                    [chromosome, name, strand, ac, start, stop, product] = (re.split('\t', line1))
                    [chromosome2, name2, strand2, ac2, start2, stop2, product2] = (re.split('\t', line2))
                    if strand == strand2:
                        if (strand == '+' and stop < start2) or (strand == '-' and start < stop2):
                            pair = (ac, ac2)
                            print('>' + pair[0] + '\n' + (fa.get(pair[0])), file=open(specie + "1.fa", "w")), print('>' + pair[1] + '\n' + (fa.get(pair[1])), file=open(specie + "2.fa", "w"))
                            blastp_cline = blastp(query=specie + '1.fa', subject=specie + '2.fa', evalue='10e-6', max_hsps=1, out=specie + 'pair.txt')
                            blastp_cline()
                            with open(specie + 'pair.txt') as out:
                                for line in out:
                                    if 'Expect' in line:
                                        E = (re.search('Expect = (.*), ', line)).group(1)
                                    elif 'Identities' in line:
                                        Identities = (re.search('Identities = .* \((.*)%\), P', line)).group(1)
                                        Positives = (re.search('Positives = .* \((.*)%\), G', line)).group(1)
                                        print(name + '\t' + ac + '\t' + name2 + '\t' + ac2 + '\t' + Identities + '\t' + Positives + '\t' + E, file=open(cwd + '/tandem_clust/' + specie + '_tandemclust' + '.tsv', 'a'))
        #### rendere uniche le sequenze in tutti i faa
                                        print('>' + pair[0] + ' ' + name + ' [' + specie + ']\n' + (fa.get(pair[0])) + '\n' + '>' + pair[1] + ' ' + name2 + ' [' + specie + ']\n' + (fa.get(pair[1])), file=open(cwd + '/tandem_clust/' + specie + '_tandemclust' + '.faa', "a"))
                                        ac_seen.append(ac), ac2_seen.append(ac2), faa_seen.update({('>' + pair[0] + ' ' + name + ' [' + specie + ']' + '\n' + (fa.get(pair[0]))): ('>' + pair[1] + ' ' + name2 + ' [' + specie + ']' + '\n' + (fa.get(pair[1])))})
                                        processed += 1
            print('- Blast: ' + str(processed) + ' tandem duplications found', end= ', ')
            os.remove(specie + '1' + '.fa'), os.remove(specie + '2' + '.fa'), os.remove(specie + 'pair' + '.txt')

            xprocessed = 0
            print('Gene 1\tAc 1\tGene 2\tAc 2\tid %\tsom %\tE-value', file=open(cwd + '/tandem/' + specie + '_tandem.tsv', 'w'))
            for line in open(cwd + '/tandem_clust/' + specie + '_tandemclust' + '.tsv', 'r'):
                [name, ac, name2, ac2, id, som, Eval] = (re.split('\t', line))
                if ac not in ac2_seen and ac2 not in ac_seen and not 'Ac 1' in line:
                    open(cwd + '/tandem/' + specie + '_tandem.tsv', 'a').write(line)
                    open(cwd + '/tandem/' + specie + '_tandem.faa', 'a').write('>' + ac + ' ' + name + ' [' + specie + ']\n' + fa.get(ac) + '\n' + '>' + ac2 + ' ' + name2 + ' [' + specie + ']\n' + fa.get(ac2) + '\n'), open(cwd + '/tandem_vs_main/query.fa', 'a').write('>' + ac + ' ' + name + ' [' + specie + ']\n' + fa.get(ac) + '\n' + '>' + ac2 + ' ' + name2 + ' [' + specie + ']\n' + fa.get(ac2) + '\n')
                    xprocessed +=1
            print('(' + str(xprocessed) + ' isolated)')

print('\n############# BLAST #################')
for file in os.listdir(cwd + '/main/'):
    if '.faa' in file and not 'phr' in file and not 'pin' in file and not 'psq' in file:
        specie = os.path.splitext(re.split('_main', file)[0])[0]
        if not os.path.exists(cwd + '/main/' + specie + '_main.faa.phr'):
            makeblastdb_cline = NcbimakeblastdbCommandline(input_file=(cwd + '/main/' + specie + '_main.faa'), dbtype='prot')
            makeblastdb_cline()
        if not os.path.exists(cwd + '/tandem_vs_main/' + specie + '_qvs.txt') and not 'phr' in file and not 'pin' in file and not 'psq' in file:
            print('Query vs ' + specie + ': ', end='')
            blastp_cline = blastp(query=(cwd + '/tandem_vs_main/query.fa'), db=(cwd + '/main/' + specie + '_main.faa'), num_threads = 10, evalue='10e-6', max_hsps=1, outfmt=7, out=(cwd + '/tandem_vs_main/' + specie + '_qvs.txt'))
            blastp_cline()
            print('done')

        seenheaderb, fastas = set(), {}
        if not os.path.exists(cwd + '/tandem_vs_main/fasta/' + specie + '_qvs.faa'):
            for line in open(cwd + '/tandem_vs_main/' + specie + '_qvs.txt', 'r').readlines():
                if '#' not in line:
                    seenheaderb.add(re.split('\t', line)[1])
            fasta_sequences = SeqIO.parse(open(cwd + '/main/' + file), 'fasta')
            for fast in fasta_sequences:
                fastas.update({fast.id: str(fast.seq)})
            for header in seenheaderb:
                for i, s in fastas.items():
                    if header in i:
                        print('>' + i + '\n' + s, file=open(cwd + '/tandem_vs_main/fasta/' + specie + '_qvs.faa', 'a'))
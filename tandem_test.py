# /usr/bin/env python

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
import re
import itertools
import sys
import csv
import os
import pandas as pd
import glob
import gzip
import shutil
from ftplib import FTP

def download():

    print('############## Downloading ################')
    cwd = os.getcwd()

    ftp = FTP('ftp.ncbi.nih.gov')
    ftp.login()

    if not os.path.exists(cwd + '/assembly_summary_refseq.txt'):
        print('Assembly summary refseq: downloading...', end=' ')
        source = ('/genomes/refseq/')
        filename = ('assembly_summary_refseq.txt')
        ftp.cwd(source)
        localfile = open(filename, 'wb')
        ftp.retrbinary("RETR " + filename, localfile.write)
        print('done')
    else:
        print('Assembly summary refseq: already exists')

    with open('species_list.txt', 'r') as list:
        with open('assembly_summary_refseq.txt', 'r') as summary:

            list_lines = list.readlines()
            version_specie = {}
            for line in summary:
                for line2 in list_lines:
                    if (re.split('\s', line2)[0] + ' ' + re.split('\s', line2)[1]) in line and re.search(
                            ('representative genome|reference genome'), line):
                        pattern = re.compile('/genomes/all/\w{3}/\d{3}/\d{3}/\d{3}/GCF_\d+.\d+_.+\t')
                        x = pattern.findall(line)
                        vers = str(x)[2:-6]
                        version_specie.update({vers: ((re.split('[\s]', line2)[0] + '_' + re.split('[\s]', line2)[1]).strip())})

            for v, s in version_specie.items():
                print(re.split('_', s)[0] + '_' + re.split('_', s)[1] + ':', end = ' ')
                version = re.split('[/]', v)[7]

                if not os.path.exists(cwd + '/query/' + s + '.faa') and not os.path.exists(cwd + '/sub/' + s + '.faa') and not os.path.exists(cwd + '/' + s + '.faa'):
                    print('downloading...', end = ' ')
                    ftp.cwd(v)
                    gtf = (version + '_genomic.gtf.gz')
                    faa = (version + '_protein.faa.gz')
                    localfile = open(gtf, 'wb')
                    ftp.retrbinary("RETR " + gtf, localfile.write)
                    localfile = open(faa, 'wb')
                    ftp.retrbinary("RETR " + faa, localfile.write)

                    with gzip.open(version + '_genomic.gtf.gz', 'r') as gtf, open((s + '.gtf'), 'wb') as localfile:
                        shutil.copyfileobj(gtf, localfile)
                        with gzip.open(version + '_protein.faa.gz', 'r') as faa, open((s + '.faa'),'wb') as localfile:
                            shutil.copyfileobj(faa, localfile)
                    os.remove(version + '_genomic.gtf.gz')
                    os.remove(version + '_protein.faa.gz')
                    print('done')
                else:
                    print('already exists')
            else:
                print('Done')
    ftp.quit()
download()

def parsing():

    print('\n################ Parsing ##################')
    cwd = os.getcwd()

    for g in os.listdir():
        if '.gtf' in g:
            string = {}
            specie = os.path.splitext(g)[0]
            if not os.path.exists(cwd + '/' + specie + '.tsv') and not os.path.exists(cwd + '/sub/' + specie + '.tsv') and not os.path.exists(cwd + '/sub/main/' + specie + '_tandem.faa'):
                with open(g) as gtf:
                    for line in gtf:
                        if 'NC_' in line:
                            if 'protein_coding' in line:
                                gene = (re.split('[\t"]', line)[9])
                                chromosome = re.split('[\t]', line)[0]
                                strand = re.split('[\t]', line)[6]
                                name = (chromosome + '\t' + gene + '\t' + strand)
                            elif 'start_codon' in line:
                                if 'protein_id' in line:
                                    ac = (re.search('; protein_id "(.*)', line)).group(1)
                                    ac = re.split('"', ac)[0]
                                    start = (re.split('\t', line))[3]
                                    product = (re.search('; product "(.*)"; protein', line)).group(1)
                            elif 'stop_codon' in line:
                                stop = (re.split('\t', line))[4]
                                string.setdefault((name), {}).update({ac: [start, stop, product]})

                print('Chromosome\tGene\tStrand\tCDS\tStart\tStop\tProduct', file=open(specie + '.tsv', "w"))
                for e, i in string.items():
                    first = list(i.items())[0]
                    out = dict(itertools.islice(i.items(), 1))
                    print(e + '\t' + '\t'.join("{}\t{}".format(k, '\t'.join(v)) for k, v in out.items()), file=open(specie + '.tsv', "a"))

                tab = pd.read_csv(specie + '.tsv', sep='\t', header=0)
                tab = tab.sort_values(['Chromosome', 'Start'], ascending=(True, True))
                tab.to_csv(specie + '.tsv', sep='\t', index=False)

                for f in os.listdir():
                    if '.faa' in f:
                        if os.path.splitext(g)[0] == os.path.splitext(f)[0]:
                            zprocessed = 0
                            header, sequence = [], []
                            fasta_sequences = SeqIO.parse(open(f), 'fasta')
                            for fasta in fasta_sequences:
                                header.append(fasta.id)
                                sequence.append(str(fasta.seq))
                            fa = dict(zip(header, sequence))

                            for e, i in string.items():
                                protid = (list(i.items())[0])[0]
                                gene = re.split('\t', e)[1]
                                if protid in header:
                                    print('>' + gene + '_' + specie + '_' + protid + '\n' + fa.get(protid), file=open(specie + '_main.faa', 'a'))
                                    zprocessed += 1
                                    sys.stdout.write('\r')
                                    sys.stdout.write((specie + ": {} main isoforms found").format(zprocessed))
                                    sys.stdout.flush()
                sys.stdout.write('\n')
        else:
            if '.faa' in g:
                specie = os.path.splitext(g)[0]
                print(specie + ': already exists')
    else:
        print('Done')
parsing()

def movefiles():

    cwd = os.getcwd()
    dirs, dirm = os.path.join(cwd + '/sub'), os.path.join(cwd + '/sub/main')
    if not os.path.exists(dirs) and not os.path.exists(dirm):
        os.mkdir(dirs), os.mkdir(dirm)

    subbies, queries = [], []
    with open('species_list.txt', 'r') as list:
        list_lines = list.readlines()
        for line in list_lines:
            if 'query' not in line:
                subbies.append((re.split(' ', line)[0] + '_' + re.split(' ', line)[1]).strip())
            else:
                queries.append((re.split(' ', line)[0] + '_' + re.split(' ', line)[1]).strip())

    for line in subbies:
        if not os.path.exists(cwd + '/sub/main/' + line + '_main.faa'):
            try:
                shutil.move(cwd + '/' + line + '_main.faa', cwd + '/sub/main/' + line + '_main.faa')
                shutil.move(cwd + '/' + line + '.tsv', cwd + '/sub/' + line + '.tsv')
                shutil.move(cwd + '/' + line + '.faa', cwd + '/sub/' + line + '.faa')
            except FileNotFoundError:
                pass
    for line in queries:
        if not os.path.exists(cwd + '/sub/main/' + line + '_main.faa'):
            try:
                shutil.move(cwd + '/' + line + '_main.faa', cwd + '/sub/main/' + line + '_main.faa')
            except FileNotFoundError:
                pass

    for g in os.listdir():
        if '.gtf' in g:
            os.remove(cwd + '/' + g)
movefiles()

def blast_query():

    cwd = os.getcwd()
    dirtc, dirt = os.path.join(cwd + '/tandem_clust'), os.path.join(cwd + '/tandem')
    if not os.path.exists(dirtc) and not os.path.exists(dirt):
        os.mkdir(dirtc), os.mkdir(dirt)

    print('\n############# Blast queries ###############')

    for f in os.listdir():
        if '.faa' in f:
            name_seen, name2_seen, faa_seen = [], [], {}
            specie = os.path.splitext(f)[0]
            if not os.path.exists(cwd + '/tandem/' + specie + '_tandem.faa') and not os.path.exists(cwd + '/tandem/' + specie + '_tandem.tsv'):
                header, sequence = [], []
                fasta_sequences = SeqIO.parse(open(f), 'fasta')
                for fasta in fasta_sequences:
                    header.append(fasta.id)
                    sequence.append(str(fasta.seq))
                fa = dict(zip(header, sequence))

                print('Gene 1\tGene 2\tid %\tsom %\tE-value', file=open(specie + '_tandemclust' + '.tsv', 'w')) ## togliere percentuale
                with open(specie + '.tsv') as tab:
                    lines = [line.rstrip() for line in tab]
                    res = [[lines[i], lines[i + 1]] for i in range(len(lines) - 1)]
                    xprocessed, vprocessed = 0, 0
                    for i in res:
                        rt = 100 * vprocessed / (len(res))
                        vprocessed += 1
                        line1, line2 = i[0], i[1]
                        [chromosome, name, strand, ac, start, stop, product] = (re.split('\t', line1))
                        [chromosome2, name2, strand2, ac2, start2, stop2, product2] = (re.split('\t', line2))
                        if strand == strand2:
                            if (strand == '+' and stop < start2) or (strand == '-' and start < stop2):
                                pair = (ac, ac2)
                                print('>' + pair[0] + '\n' + (fa.get(pair[0])), file=open(specie + "1.fa", "w"))
                                print('>' + pair[1] + '\n' + (fa.get(pair[1])), file=open(specie + "2.fa", "w"))
                                blastp_cline = blastp(query=specie + '1.fa', subject=specie + '2.fa', evalue='10e-6', max_hsps=1, out=specie + 'pair.txt')
                                stdout, stderr = blastp_cline()
                                with open(specie + 'pair.txt') as out:
                                    for line in out:
                                        if 'Expect' in line:
                                            E = (re.search('Expect = (.*), ', line)).group(1)
                                        elif 'Identities' in line:
                                            Identities = (re.search('Identities = .* \((.*)\), P', line)).group(1)
                                            Positives = (re.search('Positives = .* \((.*)\), G', line)).group(1)
                                            ## mettere ac senza parentesi ma con i tab
                                            print(name + ' (' + ac + ')' + '\t' + name2 + ' (' + ac2 + ')' + '\t' + Identities + '\t' + Positives + '\t' + E, file=open(specie + '_tandemclust' + '.tsv', 'a')) ## ac con tab senza parentesi
                                            print('>' + name + '_' + specie + '_' + pair[0] + '\n' + (fa.get(pair[0])) + '\n' + '>' + name2 + '_' + specie + '_' + pair[1] + '\n' + (fa.get(pair[1])), file=open(specie + '_tandemclust' + '.faa', "a"))
                                            name_seen.append(name), name2_seen.append(name2)
                                            faa_seen.update({('>' + name + '_' + specie + '_' + pair[0] + '\n' + (fa.get(pair[0]))): ('>' + name2 + '_' + specie + '_' + pair[1] + '\n' + (fa.get(pair[1])))})
                                            xprocessed += 1
                                            sys.stdout.write('\r')
                                            sys.stdout.write((specie + ": {} pairs found [progress: {}%]").format(xprocessed, round(rt, 1)))
                                            sys.stdout.flush()
                sys.stdout.write('\n')

                with open(specie + '_tandemclust' + '.tsv', 'r') as t, open(specie + '_tandem' + '.tsv', 'w') as tn, open(specie + '_tandem' + '.faa', 'w') as fn:
                    yprocessed = 0 - 1
                    for linet in t:
                        [name, name2, id, som, Eval] = (re.split('\t', linet))
                        name_only, name2_only = re.split(' ', name)[0], re.split(' ', name2)[0]
                        if name_only not in name2_seen and name2_only not in name_seen:
                            tn.write(linet)
                            yprocessed += 1
                    for k, v in faa_seen.items():
                        if (re.split('_', k)[0][1:]) in name_seen and (re.split('_', k)[0][1:]) in name2_seen:
                            continue
                        elif (re.split('_', v)[0][1:]) in name_seen and (re.split('_', v)[0][1:]) in name2_seen:
                            continue
                        else:
                            fn.write(k + '\n' + v + '\n')
                    print("- " + str(yprocessed) + " isolated duplications found")

                os.remove(specie + '1' + '.fa'), os.remove(specie + '2' + '.fa'), os.remove(specie + 'pair' + '.txt')

                if not os.path.exists(cwd + '/tandem/' + specie + '_tandem.faa') and not os.path.exists(cwd + '/tandem/' + specie + '_tandem.tsv') and not os.path.exists(cwd + '/tandem_clust/' + specie + '_tandemclust.faa') and not os.path.exists(cwd + '/tandem_clust/' + specie + '_tandemclust.tsv'):
                    try:
                        shutil.move(cwd + '/' + specie + '_tandem.faa', cwd + '/tandem/' + specie + '_tandem.faa')
                        shutil.move(cwd + '/' + specie + '_tandem.tsv', cwd + '/tandem/' + specie + '_tandem.tsv')
                        shutil.move(cwd + '/' + specie + '_tandemclust.faa', cwd + '/tandem_clust/' + specie + '_tandemclust.faa')
                        shutil.move(cwd + '/' + specie + '_tandemclust.tsv', cwd + '/tandem_clust/' + specie + '_tandemclust.tsv')
                    except FileNotFoundError:
                        pass
            else:
                if (specie + '.faa') in f:
                    print(specie + ': already exists')
    else:
        print('Done')
blast_query()

def blast_all():

    cwd = os.getcwd()
    dir = os.path.join(cwd + '/tandem_vs_main')
    if not os.path.exists(dir):
        os.mkdir(dir)
    tandem, tandem_main, sub, sub_main = (cwd + '/tandem/'), (cwd + '/tandem_vs_main/'), (cwd + '/sub/'), (cwd + '/sub/main/')

    print('\n########## Blast tandem vs main ###########')

    if not os.path.exists(tandem_main + 'query.txt'):
        os.chdir(tandem)
        print('Merging queries...', end=' ')
        for q in os.listdir():
            if '_tandem.faa' in q:
                    with open(q, 'r') as quer, open(tandem_main + 'query.txt', 'a') as output:
                        queries = quer.readlines()
                        for line in queries:
                            output.write(line)
        print('Done')
    else:
        print('Merging queries... already exists')

    os.chdir(sub_main)
    for m in os.listdir():
        if '_main.faa' in m and not 'phr' in m and not 'pin' in m and not 'psq' in m:
            specie = os.path.splitext(re.split('_main', m)[0])[0]

            print(specie + '_main: creating database...', end=' ')
            if not os.path.exists(specie + '_main.faa.phr'): # cerca di fare ogni volta il database per ogni iterazione, va messo fuori
                makeblastdb_cline = NcbimakeblastdbCommandline(input_file=(specie + '_main.faa'), dbtype= 'prot')
                makeblastdb_cline()
                print('Done')
            else:
                print('already exists')

            print('query_vs_' + specie + '_main: processing...', end=' ')
            if not os.path.exists(tandem_main + 'query_vs_' + specie + '.txt'):
                blastp_cline = blastp(query=(tandem_main + 'query.txt'), db=(specie + '_main.faa'), evalue='10e-6', max_hsps=1, outfmt=7, out=(tandem_main + 'query_vs_' + specie + '.txt'))
                blastp_cline()
                print('Done')
            else:
                print('already exists')
blast_all()

def parsingblastout():

    cwd = os.getcwd()
    dirs = os.path.join(cwd + '/tandem_vs_main/fasta')
    if not os.path.exists(dirs):
        os.mkdir(dirs)
    main = (cwd + '/sub/main/')
    tandem_vs_main = (cwd + '/tandem_vs_main/')
    fastad = (tandem_vs_main + '/fasta/')

    for t in os.listdir(tandem_vs_main):
        seenheaderb = set()
        if 'query' in t and not 'query.txt' in t:
            s = os.path.splitext(t)[0]
            specieb = re.split('_', s)[2] + '_' + re.split('_', s)[3]
            with open(tandem_vs_main + t) as b:
                blastout = b.readlines()
                for line in blastout:
                    if '#' not in line:
                        seenheaderb.add(re.split('\t', line)[1])
            for f in os.listdir(main):
                if specieb in f and not 'phr' in f and not 'pin' in f and not 'psq' in f:
                    sf = os.path.splitext(f)[0]
                    specief = re.split('_', sf)[0] + '_' + re.split('_', sf)[1]
                    with open(main + f) as famain:
                        fasta_sequences = SeqIO.parse(famain, 'fasta')
                        for fasta in fasta_sequences:
                            for headers in seenheaderb:
                                if fasta.id in headers:
                                    print('>' + fasta.id + '\n' + fasta.seq, file=open(fastad + 'query_vs_' + specief + '.faa', 'a'))
parsingblastout()
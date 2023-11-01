#!/usr/bin/env python
#_*_coding:utf-8_*_

import time
import argparse
import re, sys, os, platform
import itertools
from collections import Counter
import pandas as pd
from pickle import load
import matplotlib.pyplot as plt

print('\n')
print('''******************************
*      NIGAB ProbioPred      *
*  Machine Learning Tool for *
* Yeast Probiotic Prediction *
******************************''')

time.sleep(3)

def merge(file):
    con = ''
    con += '>' + file.split('.')[0] + '\n'
    with open(file, 'r') as f:
        r = f.read()
        split = r.split('>')
        for x in split:
            con += ''.join(x.split('\n')[1:])
    with open(file, 'w') as w:
        w.write(con)

def read_nucleotide_sequences(file):
    if os.path.exists(file) == False:
        print('Error: file %s does not exist.' % file)
        sys.exit(1)
    with open(file) as f:
        records = f.read()
    if re.search('>', records) == None:
        print('Error: the input file %s seems not in FASTA format!' % file)
        sys.exit(1)
    records = records.split('>')[1:]
    fasta_sequences = []
    for fasta in records:
        array = fasta.split('\n')
        header, sequence = array[0].split()[0], re.sub('[^ACGTU-]', '-', ''.join(array[1:]).upper())
        header_array = header.split('|')
        name = header_array[0]
        label = header_array[1] if len(header_array) >= 2 else '0'
        label_train = header_array[2] if len(header_array) >= 3 else 'training'
        sequence = re.sub('U', 'T', sequence)
        fasta_sequences.append([name, sequence, label, label_train])
    return fasta_sequences

def kmerArray(sequence, k):
    kmer = []
    for i in range(len(sequence) - k + 1):
        kmer.append(sequence[i:i + k])
    return kmer

def Kmer(fastas, k=3, type="DNA", upto=False, normalize=True, **kw):
    encoding = []
    header = ['#', 'label']
    NA = 'ACGT'
    if type in ("DNA", 'RNA'):
        NA = 'ACGT'
    else:
        NA = 'ACDEFGHIKLMNPQRSTVWY'

    if k < 1:
        print('Error: the k-mer value should larger than 0.')
        return 0

    if upto == True:
        for tmpK in range(1, k + 1):
            for kmer in itertools.product(NA, repeat=tmpK):
                header.append(''.join(kmer))
        encoding.append(header)
        for i in fastas:
            name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
            count = Counter()
            for tmpK in range(1, k + 1):
                kmers = kmerArray(sequence, tmpK)
                count.update(kmers)
                if normalize == True:
                    for key in count:
                        if len(key) == tmpK:
                            count[key] = count[key] / len(kmers)
            code = [name, label]
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    else:
        for kmer in itertools.product(NA, repeat=k):
            header.append(''.join(kmer))
        encoding.append(header)
        for i in fastas:
            name, sequence, label = i[0], re.sub('-', '', i[1]), i[2]
            kmers = kmerArray(sequence, k)
            count = Counter()
            count.update(kmers)
            if normalize == True:
                for key in count:
                    count[key] = count[key] / len(kmers)
            code = [name, label]
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    return encoding

def save_file(encodings, file):
    with open(file, 'w') as f:
        for line in encodings[1:]:
            line = line[1:]
            f.write('%s' % line[0])
            for i in range(1, len(line)):
                f.write(',%s' % line[i])
            f.write('\n')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="Generating Kmer feature vector for nucleotide sequences")
    parser.add_argument("--file", required=True, help="input fasta file")
    args = parser.parse_args()
    kw = {}
    print('\nPreprocessing of Sequence...\n----------------------------')
    concat = merge(args.file)
    fastas = read_nucleotide_sequences(args.file)
    print('\nCalulating Features...\n----------------------')
    encodings = Kmer(fastas, **kw)
    save_file(encodings, 'encodings.csv')
    read_feature = pd.read_csv('encodings.csv', header=None)
    del read_feature[0]
    print('\nPredicting...\n-------------')
    model = load(open('NIGAB_ProbioPred', 'rb'))
    pred = model.predict(read_feature)
    proba = model.predict_proba(read_feature)
    encode = ''
    if pred == 1:
        encode = 'Probiotic'
    else:
        encode = 'Non Probiotic'

    print(encode)
    print(proba[:, 1])
    print('\nProgram Finished\n----------------')

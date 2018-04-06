import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os
import gzip
import numpy as np

from Bio import SeqIO
from collections import defaultdict
from scipy.signal import argrelextrema


class KmerSpec:
    '''
    The class created for cutting read's-sequences into k-mers of a certain length
    and their spectrum visualisating.
    '''

    def __init__(self):
        self.fastq = ''
        self.k = 1
        self.quality = 1
        self.kmer_dict = defaultdict(int)
        self.kmer_spectr = {}

    def fastq_parse(self, fastq, k, q):
        '''
        Reading input fastq file.
        :param fastq: input .fastq file, or .fastq.gz
        :param k: k-mer length
        :param q: quality threshold
        :return nothing
        '''

        self.fastq, self.k, self.quality = fastq, k, q

        # may repack .gz
        if '.gz' in fastq:
            with gzip.open(fastq,'rt') as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    self.search_and_filter(str(record.seq), record.letter_annotations['phred_quality'])

        else:
            with open(fastq) as fq:
                for record in SeqIO.parse(fq, "fastq"):
                    self.search_and_filter(str(record.seq), record.letter_annotations['phred_quality'])

    def search_and_filter(self, seq, qual):
        '''
        Slice sequence into k-mers
        :param seq: string of 1 read
        :param qual: string of quality
        :return updated kmer_dict: {'kmer_name': [quantity]}
        '''

        length_seq = len(seq)
        for i in range(length_seq-self.k+1):
            current_kmer = seq[i:(i + self.k)]
            current_kmers_qual = qual[i:(i + self.k)]
            for q in current_kmers_qual:
                condition = True
                if q < self.quality:
                    condition = False
                    break
            if condition == True:
                self.kmer_dict[current_kmer] += 1

    def spectum_building(self):

        # k-mers frequencies
        val = list(self.kmer_dict.values())
        # frequencies occurrence
        self.kmer_spectr = {i: val.count(i) for i in list(set(val))}

    def vis_spectr(self, xmin='', xmax='', ymin='', ymax='', w=False, ax='auto'):
        '''
        Visualize k-mer spectrum. (before vizualize you should use spectrum_building method.
        :param xmin, xmax, ymin, ymax: plot borders chords
        :param w: use it if you want save plot
        :return: nothing
        '''


        sns.set(style="white", palette="BuGn_r")

        plt.title(f'Kmer spectr k={self.k} q={self.quality},\nAppr. genome size is {self.genome_size()}')
        # x = frequency, y = frequncies occurrence
        df = pd.DataFrame(list(self.kmer_spectr.items()), columns=['A', 'B'])
        df['B'].plot(logy=True)
        plt.ylabel('Frequence')
        plt.xlabel('Number of kmers')

        # noise cutoff
        plt.axvline(x=self.cut_of_search(), color='r')

        # more fine-size tuning
        if xmin and xmax and ymin and ymax:
            plt.axis([xmin, xmax, ymin, ymax])
        else:
            plt.axis(ax)
        if w:
            plt.savefig(self.fastq.split('.')[0] + f'_k{self.k}-q{self.quality}.png')

        plt.show()

    def cut_of_search(self):
        '''
        :return x[i]: local function minimum
        '''

        x = np.array(list(self.kmer_spectr.keys()))
        y = np.array(list(self.kmer_spectr.values()))

        val = y[argrelextrema(y, np.greater)[0]]

        ind = argrelextrema(y, np.less)[0]

        di = np.diff(val)

        maxval = max(di)
        j = di[:di.tolist().index(maxval)].tolist()
        i, n = 0, 0
        while i >= 0:
            i = j.pop()
        n = len(j) + 1
        n = ind[n]

        return x[n]

    def genome_size(self):
        '''
        :return: approximated genome size
        '''
        x, y, z = 0, 0, 0

        for k, v in self.kmer_spectr.items():
            if k > int():
                x += k *v
                y += k
                z += 1
        y = y/z
        return round(x/y)

if __name__ == '__main__':

    def str2bool(v):
        '''
        string variable converter
        :param v : string value, that should be boolean
        :return  : True/False
        '''

        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = argparse.ArgumentParser(description='''Kmer spectrum tool''')
    parser.add_argument('-i', '--input', help='Paste your path to input file here.', type=str)
    parser.add_argument('-k', '--kmer', help='k-mer length', default=23, type=int)
    parser.add_argument('-q', '--quality', help='quality threshold', default=25, type=int)
    parser.add_argument('-w', '--write', help='write spectrum image to .png', type=str2bool, default=True)
    parser.add_argument('-xmin', '--borderxmin', help='set graphic borders (xmin)', default='', type=str)
    parser.add_argument('-xmax', '--borderxmax', help='set graphic borders (xmax)', default='', type=str)
    parser.add_argument('-ymin', '--borderymin', help='set graphic borders (ymin)', default='', type=str)
    parser.add_argument('-ymax', '--borderymax', help='set graphic borders (ymax)', default='', type=str)
    parser.add_argument('-ax', '--axis', help='auto axis setting', default='auto', type=str)
    parser.add_argument('-gs', '--get_spectr', help='write .tsv file with spectrum for visualize somewhere else',
                        default=False, type=str2bool)

    args = parser.parse_args()
    input, kmer, quality = args.input, args.kmer, args.quality
    borderxmin, borderxmax, borderymin, borderymax = args.borderymin, args.borderxmax, args.borderymin, args.borderymax
    wr, axis, spec = args.write, args.axis, args.get_spectr

    #### for debagging ######
    # input, kmer, quality = 'SRR292678sub_S1_L001_R1_001.fastq', 20, 24
    # borderxmin, borderxmax, borderymin, borderymax = '', '', '', ''
    # wr, axis, sp = False, 'auto', True

    os.getcwd()
    A = KmerSpec()
    A.fastq_parse(input, kmer, quality)
    A.spectum_building()
    A.kmer_spectr
    if wr == True:
        A.vis_spectr(xmin=borderxmin,
                 xmax=borderxmax,
                 ymin=borderymin,
                 ymax=borderymax,
                 w=wr, ax=axis)

    if spec == True:
        for k, v in A.kmer_spectr.items():
            with open('spectr.txt', 'a') as w:
                w.write(f'{k}\t{v}\n')

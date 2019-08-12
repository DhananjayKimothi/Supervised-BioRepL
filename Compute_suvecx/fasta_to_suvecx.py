'''
Objective:
1. To covert fasta file as fastext file with labels
'''

from collections import namedtuple
import re
import numpy as np
from Bio import SeqIO
import warnings
import time
from numpy import random


############### PARAMETERS IN THE SCRIPT ###################
"""
word_size = 3  # size of each word/kmer
window_type = 'non_overlap'  # can choose overlap/non overlap
input_path # path for input fasta file
output_path # path for output fasata file
"""

class seq_slicing(object):
    def __init__(self, w_s, w_type, n_w):
        """
        w_s is word size, w_type is for selecting how you want to
        slice the sequence: 'overlap' or 'non_overlap' ; n_w is special parameter
        required only when non overlap is used, ex: for w_s = 3, n_w = '...?'
        """
        self.w_s = w_s
        self.w_type = w_type
        self.n_w = n_w

    def slices(self, seq):
        if self.w_type == 'overlap':
            words = []
            for i in range(0, len(seq) - self.w_s + 1):
                words.append(seq[i:i + self.w_s])
            return words

        if self.w_type == 'non_overlap':
            words = re.findall(self.n_w, seq)
            seq_len = len(seq)
            p = seq_len // self.w_s  # floored quotient of seq_len
            words = words[0:p]
            seq1 = seq

            xx = np.zeros(self.w_s - 1)  # to delete 1st index, del seq1[i]
            xx = xx.astype(np.int)

            words_list = []
            words2 = []
            for j in xx:
                seq1 = list(seq1)
                del seq1[j]
                seq1 = "".join(seq1)
                seq_len = len(seq1)
                words1 = re.findall(self.n_w, seq1)
                p = seq_len // self.w_s
                words1 = words1[0:p]
                words2.extend(words1)
            words.extend(words2)
            return words


class TaggedDocument(namedtuple('TaggedDocument', 'label words tags')):
    """
    A single document, made up of `words` (a list of unicode string tokens)
    and `tags` (a list of tokens). Tags may be one or more unicode string
    tokens, but typical practice (which will also be most memory-efficient) is
    for the tags list to include a unique integer id as the only tag.

    Replaces "sentence as a list of words" from Word2Vec.

    """

    def __str__(self):
        return '%s(%s,%s, %s)' % (self.__class__.__name__,  self.label, self.words,  self.tags)

class LabeledSentence(TaggedDocument):
    def __init__(self, *args, **kwargs):
        warnings.warn('LabeledSentence has been replaced by TaggedDocument', DeprecationWarning)


class LabeledLineSentence(seq_slicing):
    def __init__(self, filename,ids,w_s, w_type, n_w):
        super(LabeledLineSentence, self).__init__(w_s=w_s, w_type=w_type, n_w=n_w)
        self.filename = filename
        self.SeqDict = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))  # dictonary: keys as fasta ids
        self.ids = ids

    def __iter__(self):
        i = 0
        for key in self.ids:
            cls_lbl = key.split('_')
            key = key.split('\n')
            key = key[0]
            seq_record = self.SeqDict[key]
            seq = seq_record.seq
            seq = str(seq)
            kmer_list = self.slices(seq)  # word size
            tag = seq_record.id
            yield TaggedDocument(cls_lbl[2], kmer_list,tags=[tag])
            i = i+1


def sequence_ids(file_):
    '''
    file_: the file to be processed
    IdsL : list of sequence ids
    '''
    IdsL = []
    file__ = open(file_,'r')
    for line in file__:
        if line[0] == '>':
            line = line.replace('>', '')
            line = line.replace('\r\n','')
            IdsL.append(line)
    file__.close()

    return IdsL




def main(file_path,file_n):
####### PARAMETERS #####################
    word_size = 3  # size of each word/kmer
    window_type = 'non_overlap'  # can choose overlap/non overlap
    """
    Example; sequence --- ABCDEF
    ABCDEF --- overlap split ---- ABC BCD CDE DEF
    ABCDEF --- non overlap split --- ABC DEF
    """
    n_w = '...?'  # used when non overlaping window is selected
    IdsL = sequence_ids(filepath)
    seq_corpus = LabeledLineSentence(filepath, IdsL, w_s=word_size, w_type=window_type, n_w=n_w)
    return seq_corpus

    # Seq2Vec(seq_corpus,file_n)

if __name__ == '__main__':
    import sys
    import os
    input_path = r"C:\Users\Dhananjay\Desktop\Review_PlosOne\Clf_experiments\Sub-CellularLocation\suvecx\fasta_files"
    output_path = r"C:\Users\Dhananjay\Desktop\Review_PlosOne\Clf_experiments\Sub-CellularLocation\suvecx\fasta_files" # save folder
    for i in range(4,5):
        filepath_out = os.path.join(output_path,"Train"+str(i)+"_suvecx.fasta")
        out_file = open(filepath_out,'w')
        file = "Train"+str(i)+".fasta"
        filepath = os.path.join(input_path,file)
        seq_corpus = main(filepath,i)
        # for seq in seq_corpus:
        #     # print(seq.words)
        #     print(len(seq.words))
        for seq in seq_corpus:
            # print(seq.words)
            j = 0
            tag = seq.tags[0]
            # print(tag)
            tag = tag.split('_')
            tag = tag[2]
            tag = "__label__"+tag
            last_word = seq.words[-1]
            out_file.write(tag+' ')
            while(j<len(seq.words)):
                for kmer in seq.words:
                    out_file.write(kmer+' ')
                    j += 1
            out_file.write(last_word+"\n")

    out_file.close()

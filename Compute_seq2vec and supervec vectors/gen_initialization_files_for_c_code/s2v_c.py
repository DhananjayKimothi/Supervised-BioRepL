'''
This code is same as s2v except the rerencing to sequence ids-- so as to match the c implimentation

Objective:
1. For labelling the sequence corpus train, test file are passed to the labled sequence module

'''


from collections import namedtuple
import re
import numpy as np
from Bio import SeqIO
import warnings
import time
import pickle as pkl
from MakeVocab1 import BuildVocab
from DocVecsArray import DocvecsArray
from numpy import random

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
            # yield LabeledSentence(1, kmer_list, tags=['SEQ_%s' % i])
            i = i+1




def sequence_ids(file_train_p):
    IdsL = []
    file_train = open(file_train_p,'r')
    for line in file_train:
        if line[0] == '>':
            line = line.replace('>', '')
            line = line.replace('\r\n','')
            IdsL.append(line)
    file_train.close()

    return IdsL


class Seq2Vec(BuildVocab):
    def __init__(self,sequences,path_no,file_no):
        self.min_count = 0
        self.sample = 0
        self.max_vocab_size = 100000
        self.docvecs = DocvecsArray()
        self.hs = 1
        self.vector_size = 100  # check this --its a vector size
        self.layer1_size = self.vector_size
        self.negative = 0
        self.sorted_vocab = 0  # want to sort vocab on the basis of frequency ?
        self.null_word = 0  # check this
        self.dm = 1
        self.dm_concat = 0
        self.seed = 0
        self.random = random.RandomState(self.seed)
        self.hashfxn = hash
        self.total_words = 0
        self.build_vocab(sequences)  # the function is defined in BuildVocab
        self.path_no = path_no
        self.file_no = file_no
        """
        Saving initializations

        """

        # path = r"E:\SuperVec_Codes_25072018\Data\50_classes\DifferentDimension"


        x = self.file_no
        index2word_py = open(output_path + '\index2word' + str(x) + '.pkl', 'w')
        pkl.dump(self.index2word, index2word_py)

        vocab_py = open(output_path + '\\vocab' + str(x) + '.pkl', 'w')  # kmers, paths and codes
        pkl.dump(self.vocab, vocab_py)

        doctag_py = open(output_path + '\doctag' + str(x) + '.pkl',
                         'w')  # doctag initialization
        pkl.dump(self.docvecs.doctag_syn0, doctag_py)

        kmer_py = open(output_path + '\kmer' + str(x) + '.pkl', 'w')  # kmer initialization
        pkl.dump(self.syn0, kmer_py)

        vocab_py.close()
        doctag_py.close()
        kmer_py.close()
        index2word_py.close()



def main(path_no,file_no):
    # file = "T" + str(file_no)+".fasta"
    # file = "T"+str(file_no)+".fasta"
    file = "T_8.fasta"
    # file = "pos.fasta"

    word_size = 3  # size of each word
    window_type = 'non_overlap'  # can choose overlap or non overlap
    n_w = '...?'  # used when non overlaping window is selected
    # filepath = r"E:\inferringPPI\DATA\Guo Data\E.coli\interacting\T2.fasta"
    filepath = os.path.join(input_path,file)
    # filepath = r"E:\SuperVec_Codes_25072018\Data\50_classes\DifferentDimension\T0.fasta"

    IdsL = sequence_ids(filepath)
    print(IdsL)
    seq_corpus = LabeledLineSentence(filepath, IdsL, w_s=word_size, w_type=window_type, n_w=n_w)
    Seq2Vec(seq_corpus,path_no,file_no)


input_path = r"C:\Users\Dhananjay\Desktop\HSuperVec_BLAST\Largest_8class"
output_path = r"C:\Users\Dhananjay\Desktop\HSuperVec_BLAST\Largest_8class" # save folder



if __name__ == '__main__':
    import sys
    import os
    # param()
    #file_no = 0
    # path_no = sys.argv[1]
    # file_no = sys.argv[2]
    path_no = 0
    for i in range(0,1):
        file_no = i
        #print(file_no)
        main(path_no,file_no)
        # param()
        #filepath_train = "E:\Implementations\DATA\Top25\Train1.fasta"
        #filepath_test = "E:\Implementations\DATA\Top25\Test1.fasta"

        # for i in range(5,6):
        #     main(i)

### Run vocab_text.py
""""

system arguments ---- path and file number

"""

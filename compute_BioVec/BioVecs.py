"""
--- generates the sequence vectors from database.fasta and queries.fasta
--- Algo:
    iterate through the sequence
    Split the sequences
    Search the kmer rep from the db and add

--- input fasta file

--- Output :
Database['num_families'] = num_families
Database['family_names'] = family_name
Database['family_sizes'] = family_size
Database['sequence_vectors'] = sequences_vecs
Database['info'] =

"""

import re
import numpy as np
vec_size = 100


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

def sequence_ids(file_path):
    IdsL = [] # L is for denoting list
    file = open(file_path,'r')
    for line in file:
        if line[0] == '>':
            line = line.replace('>', '')
            line = line.replace('\r\n', '')
            IdsL.append(line)
    file.close()
    return IdsL

def other_stats(file_path):
    """
    Database['num_families'] = num_families
Database['family_names'] = family_name
Database['family_sizes'] = family_size
Database['sequence_vectors'] = sequences_vecs
Database['info'] =


    :param file_path:
    :return:
    """

    class_of_seqs = []
    family_sizes = []
    IdsL = sequence_ids(file_path)

    for id in IdsL:
        id = id.split("\n")
        id = id[0]
        id = id.split('_')
        class_of_seqs.append(id[2])
    class_of_seqs = [int(i) for i in class_of_seqs]
    family_names = set(class_of_seqs)
    family_names = list(family_names)
    family_names = sorted(family_names)
    for cls in family_names:
        family_sizes.append(class_of_seqs.count(cls))
    num_families = len(family_names)

    return {'num_families':num_families, 'family_names': family_names, 'family_sizes':family_sizes}

def add_kmers(kmer_list, kmer_vecs, vec_size):
    """
        -- search the kmers from the list and add the
        corresponding vector available.
        :return:
        """
    column_indxs = kmer_vecs.columns.get_values() # ['Unnamed: 0' '<unk>' 'AAA' ..., 'ZPG' 'ZSY' 'ZZK']
    seq_vec = np.zeros(vec_size) # list of lists
    for kmer in kmer_list:
        if kmer in column_indxs:
            temp = kmer_vecs[kmer]
        else:
            temp = kmer_vecs['<unk>']
        seq_vec = seq_vec + np.asarray(temp)
    return seq_vec

def split_sequence(filename,vec_size):
    """
    Split the sequence
    :return: list of kmers
    """
    # read sequences -- use of Bio library
    from Bio import SeqIO
    SeqDict = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))

    import pandas
    # E:\SuperVec_Codes_25072018\Python_scripts\compute_BioVec
    path = r'E:\SuperVec_Codes_25072018\Python_scripts\compute_BioVec\protVec_100d_3gramsDataframe.csv'
    X = pandas.read_csv(path)

    w_s = 3
    w_type = 'overlap'
    n_w = '...?'
    slicing = seq_slicing(w_s, w_type, n_w)  # make an object
    ids = sequence_ids(filename)
    seq_vecs = []
    print("total_sequences {}".format(len(ids)))
    i = 1
    for key in ids:
        # print(len(ids), i)
        key = key.split('\n')
        key = key[0]
        seq_record = SeqDict[key]
        seq = seq_record.seq
        seq = str(seq)
        kmer_list = slicing.slices(seq)
        temp = add_kmers(kmer_list,X,vec_size)
        temp = list(temp)
        seq_vecs.append(temp)
        i = i+1
    return seq_vecs




def main(input_fasta_file, output_BioVecs):
    import os
    import time
    # in_Dir_path = r"E:\SuperVec_Codes_25072018\Data\pairs_from1-200\pairs_fasta_files"
    # out_Dir_path = r"E:\SuperVec_Codes_25072018\Data\pairs_from1-200"
    # input_fasta_file = os.path.join(in_Dir_path, 'T_0.fasta')
    # output_BioVecs = os.path.join(out_Dir_path, 'T_0_BioVecs.pkl')
    strt = time.time()
    BioVecs = other_stats(input_fasta_file)
    sequence_vectors = split_sequence(input_fasta_file,vec_size)
    BioVecs['sequence_vectors'] = sequence_vectors
    end = time.time()
    print("time taken for generating BioVecs", round((end - strt),3))

    import pickle as pkl
    DB = open(output_BioVecs, 'wb')
    pkl.dump(BioVecs, DB)
    DB.close()

    #print(x)



if __name__ == "__main__":
    import time
    strt = time.time()
    main()
    end = time.time()
    print("time taken for generating BioVecs", round((end - strt),3))

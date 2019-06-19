
import pandas as pd
import pickle as pkl
from Bio import SeqIO
import os
import random

# read dataframes

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
            line = line.replace('\n','')
            IdsL.append(line)
    file__.close()

    return IdsL


def id_seq_dict(file):
    """
    Tuple (Dictionary -- keys : sequenceIDs and value:corresponding sequences, list of sequence ids)
    """
    SeqDict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
    IdsL = sequence_ids(file)
    return SeqDict, IdsL

def gen_txt_file(IdsL,sample_splitsL,SeqDict,path):
    # len(sample_splitsL)
    for i in range(0,5):
        train = []
        test = sample_splitsL[i]
        print(len(test))
        for j in range(0,len(sample_splitsL)):
            print(j)
            if j != i:
                print(sample_splitsL[j][0:5])
                train.extend(sample_splitsL[j])
        print(len(train))
        file_tst = os.path.join(path,"Tst"+str(i)+".fasta")
        file_train = os.path.join(path,"Train"+str(i)+".fasta")
        file_ = open(file_tst,'w')
        file__ = open(file_train,'w')
        for indx in test:
            id = IdsL[indx]
            print(id)
            file_.write(">"+id+"\n")
            file_.write(str(SeqDict[id].seq)+"\n")
        for indx in train:
            id = IdsL[indx]
            print(id)
            file__.write(">"+id+"\n")
            file__.write(str(SeqDict[id].seq)+"\n")
        file_.close()
        file__.close()

def cross_valid(SeqDict,IdsL):
    sample_list = range(0,len(IdsL))
    # sample_list = range(0,10)
    random.shuffle(sample_list)
    n_test_samples_one_split = len(sample_list)//5

    strt = 0
    end = strt+n_test_samples_one_split

    sample_splitsL = []
    for i in range(0,4):
        sample_splits = sample_list[strt:end]
        sample_splitsL.append(sample_splits)
        strt = end
        end = end+n_test_samples_one_split
    sample_splitsL.append(sample_list[strt:])
    return sample_splitsL

    # print(len(sample_splitsL))
        # sample_list[y:y+x]

        # print(sample_list[0:20])
random.seed(0)
if __name__ == "__main__":
    # dataframe -- protein1
    path = r"C:\Users\Dhananjay\Desktop\Review_PlosOne\enzyme-nonenzyme\5-foldCrossValidation\non-enzyme"
    pathA = os.path.join(path, "nonenzyme_preprocessed.fasta")
    SeqDict, IdsL = id_seq_dict(pathA)
    # print(SeqDict.keys()[0])
    # print(IdsL[0:10])
    # print()
    sample_splitsL = cross_valid(SeqDict, IdsL)
    # print(sample_splitsL[0][0:5])
    # print(sample_splitsL[1][0:5])
    # print(sample_splitsL[2][0:5])
    # print(sample_splitsL)
    gen_txt_file(IdsL,sample_splitsL,SeqDict,path)

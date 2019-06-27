import os
import pickle as pkl
import numpy as np
from collections import defaultdict

def sequence_ids(file_path,file_name):
    IdsL = [] # L is for denoting list
    file_path = os.path.join(file_path,file_name+".fasta")
    file = open(file_path,'r')
    for line in file:
        if line[0] == '>':
            IdsL.append(line)
    file.close()
    return IdsL

def fasta_to_dict(file_path,file_name):
    D = defaultdict()
    file_path = os.path.join(file_path,file_name+".fasta")
    file = open(file_path, 'r')
    for line in file:
        if line[0] == '>':
            temp_id = line
            D[temp_id] = ""
            continue
        else:
            temp = line
            #temp = line.split('\n')
            D[temp_id] = D[temp_id]+temp

    file.close()
    return D

def gen_file(main_ids_dict, db_ids_query, filename):
    """
    :main_ids_dict: ids of main file ; DB
    :ids_dict: To be selected
    :filename : filename to be used to store the sequence from selected classes
    :clses
    """
    file = open(filename,'w')
    for id in db_ids_query:
        file.write(id)
        seq = main_ids_dict[id]
        file.write(seq)
    file.close()


def load_data(Dir_path, file_name):
    """
    load database sequence vectors
    """

    temp = open(os.path.join(Dir_path, file_name + ".pkl"), 'rb')
    temp_ = pkl.load(temp)
    vecs = temp_['sequence_vectors']
    vecs = np.asarray(vecs)
    family_sizes = temp_['family_sizes']
    temp.close()

    return (vecs,family_sizes)


def extract_labels(Dir_path, file_name):
    """
    script for extracting labels from a fasta file
    """

    file_path = os.path.join(Dir_path,file_name+".fasta")
    file_in = open(file_path,'r')
    # test_ = file_in.readlines()
    file_iter = iter(file_in) # iterator
    labels = []
    while True:
        try:
            line = file_iter.next()
            if line[0] == ">":
                line = line.replace('\n','')
                line = line.split('_')
                label = int(line[2])
                labels.append(label)
        except StopIteration:
            break
    return(labels)

def unique_labels(labels):
    #To check unique labels
    labels = set(labels)
    labels = list(labels)
    return labels

def inter_pre_values_per_query(qry_sub_distance, total_relevant, db_family_index):
    """
    :param qry_sub_distance: 2D array of subject indices based on their distance with the query.
    :return: interpolated precision value at 100 points
    """

    precisionAt = np.zeros(101)
    found = 0
    retreived = 0
    # qry_sub_distance  # ranked class list to be processed for calculating precision recall values

    while (found < total_relevant):

        if (db_family_index[0] <= qry_sub_distance[0][retreived] < db_family_index[1]):
            found = found + 1
        retreived = retreived + 1
        recall = float(found) / total_relevant
        precision = float(found) / retreived
        recall = int(np.ceil(recall * 100))
        if precisionAt[recall] < precision:
            precisionAt[recall] = precision

    for j in range(99, -1, -1):
        precisionAt[j] = max(precisionAt[j], precisionAt[j + 1])

    return precisionAt




def select_node_weights(score, root_labels, unique_labels_persplit, nbrs = 20):
    ranked_indices = np.argsort(score)
    # print(ranked_indices)
    neighbors = ranked_indices[0][0:nbrs]

    # print(root_labels)
    # print(neighbors)
    out = map(root_labels.__getitem__,neighbors)
    out = list(out)
    right = 0
    left = 0
    for indx in out:
        if indx in unique_labels_persplit[0]:
            left += 1
        else:
            right +=1
    # print(left,right)
    left_weight = left/float(nbrs)
    right_weight = right/float(nbrs)

    return (left_weight, right_weight)

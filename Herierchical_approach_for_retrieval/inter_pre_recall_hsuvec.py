"""
This script compute precision-recall values for H-SuVec
The computation for calculating precision recall values is
done per sequence.

"""

import numpy as np
import time
import pickle as pkl
import os
from sklearn.metrics.pairwise import cosine_distances
from utils_herierchical_approach import *


def dist(q, db):
    score_list = cosine_distances(q, db)
    # X = np.argsort(score_list)
    return (score_list)

def main(Dir_path, DB_file1, Q_file1, DB_file2, Q_file2, DB_file3, Q_file3,root_labels,unique_labels_persplit):
    strt = time.time()
    root_weight = 0.5
    # w = [0.75,0.125,0.125]
    ranked_db_sequecnces_for_queries = []  # indices of database sequences for all queries list of lists
    ranked_db = {}


################# Load the dabase and query sequence vectors computed at different nodes ##############
    db1, db_family_sizes1 = load_data(Dir_path,DB_file1)
    q1, q_family_sizes1 = load_data(Dir_path,Q_file1)

################## right-split ############################################
    db2, db_family_sizes2 = load_data(Dir_path,DB_file2)
    q2, q_family_sizes2 = load_data(Dir_path,Q_file2)

############## left-split ###########################################

    db3, db_family_sizes3 = load_data(Dir_path,DB_file3)
    q3, q_family_sizes3 = load_data(Dir_path,Q_file3)

################

    # len(q)

    x = 0 # for itera
    y = 0
    precision_avg = np.zeros(11)
    families =  len(q_family_sizes1)
    # families = 1
    strt = time.time()
    for i in range(0, families):
        precisionClass = np.zeros(101)
        total_relevant = db_family_sizes1[i]
        db_family_index = [y,y+db_family_sizes1[i]] # starting and the last index for i^th family
        y = y+db_family_sizes1[i]
        # print db_family_index
        #
        queries_in_family = q_family_sizes1[i]
        # queries_in_family = 20
        for j in range(0, queries_in_family):

            ind = x+j
            # print ind
            # computing distance and sorting the subjects
            score1 = dist(q1[[ind], :], db1)
            left_weight, right_weight = select_node_weights(score1,root_labels,unique_labels_persplit)

            w = [root_weight, (1-root_weight)*left_weight,(1-root_weight)* right_weight]


            # classification based
            # if (left_weight > right_weight):
            #     w = [root_weight, 0.5, 0]
            # else:
            #     w = [root_weight, 0, 0.5]
            #

            # print(w)
            score2 = dist(q2[[ind], :], db2)
            score3 = dist(q3[[ind], :], db3)
            ranked_db_sequences_per_query = (w[0]*score1 + w[1]*score2 + w[2]*score3)
            ranked_db_sequences_per_query = np.argsort(ranked_db_sequences_per_query)


            precisionAt = inter_pre_values_per_query(ranked_db_sequences_per_query, total_relevant, db_family_index)

            precisionClass = precisionClass + precisionAt


        x = x+q_family_sizes1[i]
        # print(precisionClass)
        precisionClass = precisionClass / float(q_family_sizes1[i])
        precisionClass_ = [round(val, 4) for val in precisionClass]
        writer.writerow(precisionClass_[0:101:10])
        print("precision_per_class", precisionClass_[0:101:10])
        precision_avg += precisionClass[0:101:10]

    precision_avg = list(precision_avg / float(families))
    precision_avg = [round(val, 4) for val in precision_avg]
    writer.writerow(precision_avg)
    print("precision_avg", precision_avg)


if __name__ == "__main__":
    import csv

    Dir_path = r"C:\Users\Dhananjay\Desktop\Review_PlosOne\largest_eight_class"

    DB_file1 = "T_8_fasttext_vec"
    Q_file1 = "Q_DB_8_fasttext_vec"

    DB_file2 = "T_8_1_fasttext_vec"
    Q_file2 = "Q_DB_8_1_fasttext_vec"

    DB_file3 = "T_8_2_fasttext_vec"
    Q_file3 = "Q_DB_8_2_fasttext_vec"

################# Extract labels from fasta file ########################
    root = "T_8"
    left_split_db = "T_8_1"
    right_split_db = "T_8_2"

    root_labels = extract_labels(Dir_path,root) # for individual sequences
    left_split_unique_labels = unique_labels(extract_labels(Dir_path,left_split_db)) # classes in split1
    right_split_unique_labels = unique_labels(extract_labels(Dir_path,right_split_db)) # classes in split2
    unique_labels_persplit = [left_split_unique_labels,right_split_unique_labels]
##################### To write results in csv file ####################

    result_file = "H_fasttext"
    result_file = open(os.path.join(Dir_path, result_file + ".csv"), 'a')
    writer = csv.writer(result_file, lineterminator='\n')


    main(Dir_path, DB_file1, Q_file1, DB_file2, Q_file2, DB_file3, Q_file3,root_labels,unique_labels_persplit)


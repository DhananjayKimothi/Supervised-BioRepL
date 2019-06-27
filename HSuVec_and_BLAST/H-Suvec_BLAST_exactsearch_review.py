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


def main(path_load_vec, path_save_data,DB_file1, Q_file1, DB_file2, Q_file2, DB_file3, Q_file3, root_labels,unique_labels_persplit,DB_file_dict,DB_IDs):
    strt = time.time()
    root_weight = 0.5
    # w = [0.75,0.125,0.125]
    nn = 2000
    ranked_db_sequecnces_for_queries = []  # indices of database sequences for all queries list of lists
    ranked_db = {}


################# Load the dabase and query sequence vectors computed at different nodes ##############
    db1, db_family_sizes1 = load_data(path_load_vec,DB_file1)
    q1, q_family_sizes1 = load_data(path_load_vec,Q_file1)

################## right-split ############################################
    db2, db_family_sizes2 = load_data(path_load_vec,DB_file2)
    q2, q_family_sizes2 = load_data(path_load_vec,Q_file2)

############## left-split ###########################################

    db3, db_family_sizes3 = load_data(path_load_vec,DB_file3)
    q3, q_family_sizes3 = load_data(path_load_vec,Q_file3)

################

    x = 0 # for itera
    y = 0
    families =  len(q_family_sizes1)
    family_strt = 0
    family_last = 1
    strt = time.time()
    for family in range(family_strt, family_last):
        q_family_size = q_family_sizes1[family]
        print(q_family_size)
        x = sum(q_family_sizes1[0:family])
        for qry_id in range(0, q_family_size):
            db_query_file = os.path.join(path_save_data,"db_query"+str(family+1)+"_"+str(qry_id)+".fasta") # db file generated for each query
            ind = x+qry_id # global id

            # computing distance and sorting the subjects
            score1 = dist(q1[[ind], :], db1)
            left_weight, right_weight = select_node_weights(score1,root_labels,unique_labels_persplit)

            w = [root_weight, (1-root_weight)*left_weight,(1-root_weight)* right_weight]
            score2 = dist(q2[[ind], :], db2)
            score3 = dist(q3[[ind], :], db3)
            ranked_db_sequences_per_query = (w[0]*score1 + w[1]*score2 + w[2]*score3)
            ranked_db_sequences_per_query = np.argsort(ranked_db_sequences_per_query)

            ranked_db_ids_per_query  = [DB_IDs[id] for id in ranked_db_sequences_per_query[0]]
            #print(ranked_db_ids_per_query[0:3])
            ranked_db_ids_per_query = ranked_db_ids_per_query[0:nn]
            gen_file(DB_file_dict, ranked_db_ids_per_query, db_query_file)

        x = x+q_family_sizes1[family]


if __name__ == "__main__":
    import csv

    path_load_vec  = r"/home/iiitd/Dhananjay/200_classes/H_SuVec" # for db, query Vectors
    path_save_data = r"/home/iiitd/Dhananjay/HSuperVec_BLAST/Retrieval_results_200classes_review/db_generated_for_each_query" # path for storing db generated for each query
    path_load_metadata  = r"/home/iiitd/Dhananjay/200_classes/H_SuVec" # for extracting sequence ids

    DB_file_fasta = "T.fasta" # Database file



    # DB_file1 = "T0_fasttext_vec"
    # Q_file1 = "Q_DB_fasttext_vec"
    #
    # DB_file2 = "T0_fasttext_model_100_1_vec"
    # Q_file2 = "Q_DB_fasttext_model_100_1_vec"
    #
    # DB_file3 = "T0_fasttext_model_100_2_vec"
    # Q_file3 = "Q_DB_fasttext_model_100_2_vec"

    DB_file1 = "T_infer_pt5"
    Q_file1 = "Q_DB_pt5"

    DB_file2 = "T_infer_split1"
    Q_file2 = "Q_DB_infer_split1"

    DB_file3 = "T_infer_split2"
    Q_file3 = "Q_DB_infer_split2"

################# Extraxct labels from fasta file ########################

    # root = os.path.join(Dir_path1,"T")
    # left_split_db = os.path.join(Dir_path1,"T_100_1")
    # right_split_db = os.path.join(Dir_path1,"T_100_2")

    root = "T"
    left_split_db = "T_100_1"
    right_split_db = "T_100_2"



    root_labels = extract_labels(path_load_metadata,root) # for individual sequences
    left_split_unique_labels = unique_labels(extract_labels(path_load_metadata,left_split_db)) # classes in split1
    right_split_unique_labels = unique_labels(extract_labels(path_load_metadata,right_split_db)) # classes in split2
    unique_labels_persplit = [left_split_unique_labels,right_split_unique_labels]

##################### for extracting DB sequence ids #######################
    DB_file_dict = fasta_to_dict(path_load_metadata,root) # keys seq_ids , values sequence
    DB_IDs = sequence_ids(path_load_metadata,root)

    main(path_load_vec, path_save_data,DB_file1, Q_file1, DB_file2, Q_file2, DB_file3, Q_file3, root_labels,unique_labels_persplit,DB_file_dict,DB_IDs)

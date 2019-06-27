"""
Querying a database in vector space

calculate distance
sort the distances and corresponding indices


:: for commenting use REM or ::
:: file = DB/Q_SuV, DB/Q_S2V

for /L %%a in (0,1,124) do ( python pre-recall.py "E:\SeqReterivalWork\experiment2\four_classes" "DB_SuV%%a" "Q_SuV%%a" "result_file"
   echo %%a)
"""

import numpy as np
import time
import pickle as pkl
import os
from sklearn.metrics.pairwise import cosine_distances
from Bio import SeqIO

def dist(q, db):
    score_list = cosine_distances(q, db)
    X = np.argsort(score_list)
    return (X)

max_ind_of_a_label_in_db = []
def max_index_list(db_family_sizes):
    first_ind = db_family_sizes[0]-1
    max_ind_of_a_label_in_db.append(first_ind)
    for val in db_family_sizes[1:]:
        max_ind_of_a_label_in_db.append(first_ind+val)
        first_ind = first_ind+val
    return max_ind_of_a_label_in_db

def sequence_ids(file_train_p):
    IdsL = []
    file_train = open(file_train_p,'r')
    # file_train = open(os.path.join(Dir_path, file_train_p),'r')
    for line in file_train:
        if line[0] == '>':
            line = line.replace('>', '')
            # line = line.replace('\n','')
            IdsL.append(line)
    file_train.close()

    return IdsL

# ---------------------------- BLAST ----------------------------------
# Make DB
import os
import time
import subprocess
def make_database(file_in, file_out):
    strt_time = time.time()
    #file_in = root_path+"\DB1.fasta"
    #file_out = root_path+"\DB1_db"
    make_db = "makeblastdb -in" + " " + file_in + " "+"-dbtype prot -out " + file_out
    os.system(make_db)
    #subprocess.Popen(make_db)
    # print make_db
    # os.popen(make_db)


def Query_database(db, qry, rslt):
    # start_time = time.clock()

    # Blast a file with the database

    # db = root_path+"\DB1_db"
    # qry = root_path+"\Q1.fasta"
    # rslt = root_path+"\\rslt_xx.txt"

    # outfmt = "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
    # outfmt = "-outfmt 0"
    outfmt = "-outfmt 10"

    max_target = 100000

    query_db = "blastp -db" + " " + db + " " + "-query" + " " + qry + " " + "-out" + " " + rslt + " " + "-evalue 10" + " " + outfmt + " " + "-max_hsps 1" + " " + "-max_target_seqs" + " " + str(
        max_target)
    #os.popen(query_db)
    os.system(query_db)



if __name__ == "__main__":
    import sys
    import time
    from collections import defaultdict
    make_db = 0
    querying = 1

    Dir_path = r"/home/iiitd/Dhananjay/HSuperVec_BLAST/Retrieval_results_200classes_review/db_generated_for_each_query" # queries
    Dir_path_out = r"/home/iiitd/Dhananjay/HSuperVec_BLAST/Retrieval_results_200classes_review/blast_dbs" # corresponding dbs
    Query_file = r"/home/iiitd/Dhananjay/HSuperVec_BLAST/Retrieval_results_200classes/Q_DB.fasta" # query fasta file
    Dir_path_result = r"/home/iiitd/Dhananjay/HSuperVec_BLAST/Retrieval_results_200classes_review/results_blast"


    Query_Ids = sequence_ids(Query_file)
    cls_wise_queries = defaultdict(list) # key: class; value: sequence id

    for id in Query_Ids:
        cls = id.split('_')
        cls = int(cls[2])
        cls_wise_queries[cls].append(id)

    if make_db:
        logfile = os.path.join(Dir_path_result,"make_db_blast_.txt")
        file = open(logfile,'w')
        # print(cls_wise_queries[78])
        families = len(cls_wise_queries)
        families = 1
        for family in range(0, families):
            for query in range(0,len(cls_wise_queries[family+1])):
                file_in = os.path.join(Dir_path, "db_query"+str(family+1)+"_"+str(query)+".fasta")
                file_out = os.path.join(Dir_path_out, "db_query"+str(family+1)+"_"+str(query))
                # print(str(family+1)+"_"+str(query))
                strt2 = time.time()
                make_database(file_in, file_out)
                end2 = time.time()
                print("db_creation_time {}".format(round(end2-strt2,4)))
                file.write("db_time family{} query{} {}\n".format(family, query, round(end2-strt2,4)))
        file.close()

    if querying:
        Dir_path = r"/home/iiitd/Dhananjay/HSuperVec_BLAST/Retrieval_results_200classes/queries" #queries
        Dir_path_result = r"/home/iiitd/Dhananjay/HSuperVec_BLAST/Retrieval_results_200classes_review/results_blast" #
        logfile = os.path.join(Dir_path_result,"querying_time.txt")
        file = open(logfile,'w')
        total_time = 0
        families = len(cls_wise_queries)
        families = 1
        for family in range(0,families):
            for query in range(0,len(cls_wise_queries[family+1])):
                # print(str(family+1)+"_"+str(query))
                # db_blast = r"/home/iiitd/Dhananjay/HSuperVec_BLAST/Retrieval_results_200class_review/class1"
                db_blast = os.path.join(Dir_path_out, "db_query"+str(family+1)+"_"+str(query))
                rslt =  os.path.join(Dir_path_result, "rslt"+str(family+1)+"_"+str(query))
                q_new_file_path = os.path.join(Dir_path, "q"+str(family+1)+"_"+str(query))
                strt2 = time.time()
                Query_database(db_blast, q_new_file_path, rslt)
                end2 = time.time()
                total_time += end2-strt2
                print(family, query)
                #print("time taken for querying blast on the complete dataset is {}".format(round(end2-strt2,4)))
                file.write("family {}, query {} -- {}\n".format(family+1, query,round(end2-strt2,4)))
        print("time taken for querying blast on the complete dataset is {}".format(round(total_time,4)))
        file.write("total time taken for querying blast on the complete dataset is {}".format(round(total_time,4)))
        #file.write("Average querying time is {}".format(round(total_time/1108,4)))
        file.close()

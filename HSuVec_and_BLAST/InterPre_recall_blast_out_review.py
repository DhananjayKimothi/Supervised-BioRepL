"""
This script is same as InterPre_recall_blast_out except that it employ
iterator for navigating through the text file which is useful in case
the text file is very big.
"""

"""
1. csv file with the query and subject id
2. Relevant sequence in this case are the ids which has same familiy indicated at the end of sequence id
    example - Seq_0_1 ---- > split_label = ['Seq', '0', '1']; split_label[2] ---> class; split_label[1] ---> number
    cls_ids 1,2....25
3. Steps:
--- inputs -- query and the class
            --- noumber of relevant document per class

the query labels are ordered as the class,
---

4 Output:
--- Average interpolated precision- recall values for each family

# blast_out1.txt --- 18412203 lines
File format:
Seq_0_1,Seq_626_1,80.93,367,70,0,1,367,1,367,0.0,610

Seq_0_1,Seq_623_1,66.94,369,118,3,3,367,2,370,6e-178,503

Seq_0_1,Seq_3042_1,64.77,369,128,2,1,367,1,369,3e-174,493

Seq_0_1,Seq_3079_1,59.56,366,144,3,4,367,3,366,2e-151,435

Seq_0_1,Seq_621_1,57.10,366,154,2,4,367,3,367,7e-151,435

"""

import time
import numpy as np
from collections import defaultdict
import csv

def sequence_ids(file_train_p, file_test_p):
    Train_IdsL = []
    Test_IdsL = []
    file_train = open(file_train_p,'r')
    file_test = open(file_test_p,'r')
    for line in file_train:
        if line[0] == '>':
            line = line.replace('>', '')
            line = line.replace('\r\n','')
            Train_IdsL.append(line)
    for line in file_test:
        if line[0] == '>':
            line = line.replace('>', '')
            line = line.replace('\r\n','')
            Test_IdsL.append(line)
    file_train.close()
    file_test.close()

    return Train_IdsL, Test_IdsL


if __name__ == '__main__':
    strt = time.time()

    filepath_train = r"/home/iiitd/Dhananjay/200_classes/T0.fasta" # db i think
    filepath_test = r"/home/iiitd/Dhananjay/200_classes/Q_DB.fasta" # query file i think

    result_file = open(r"/home/iiitd/Dhananjay/HSuperVec_BLAST/Retrieval_results_200classes_review/results_blast/class_wise_blast.csv",'wb')

    wr = csv.writer(result_file, dialect = 'excel')
    Subject_Ids, Query_Ids = sequence_ids(filepath_train, filepath_test)

    cls_wise_queries = defaultdict(list)
    for id in Query_Ids:
        cls = id.split('_')
        cls = int(cls[2])
        cls_wise_queries[cls].append(id)

    relevant_seqs_clswise = defaultdict(list)

    for id in Subject_Ids:
        cls = id.split('_')
        cls = int(cls[2])
        relevant_seqs_clswise[cls].append(id)

    #print ("Relevant sequences class 22", len(relevant_seqs_clswise[22]))

    # Read File
    families = len(cls_wise_queries)
    #print(cls_wise_queries[2])
    families = 1
    for family in range(0, families):
        precisionClass = np.zeros(101)
        precision_avg = np.zeros(11)
        queries = cls_wise_queries[family+1]
        n_queries = len(queries)
        print(n_queries)
        #n_queries = 1
        for query in range(0,n_queries):
            file = open(r"/home/iiitd/Dhananjay/HSuperVec_BLAST/Retrieval_results_200classes_review/results_blast/rslt"+str(family+1)+"_"+str(query), 'r')
            txt_file = iter(file)
            #print(txt_file)
            line_no = 0
            #queries = [queries[0], queries[1], queries[2]]
            qry = queries[query]
            target_idL = []
            target = []  # list of target classes
            precisionAt = np.zeros(101)
            # for line in txt_file[strt_file: end_file]:
            q_retrival = 1 # indicator that there are documents retrieved corresponding to a query
            while (q_retrival == 1):
                try:
                    line = next(txt_file)
                    #print(line)
                except StopIteration:
                    break


                qry_rslt = line.split(',')
                # line_no = line_no+1
                # id
                qry_id = qry_rslt[0]
                target_id = qry_rslt[1]
                #print qry_id
                qry_cls = qry_id.split('_')
                qry_cls = qry_cls[2]
                qry = qry.split('\n')
                qry = qry[0]
                if qry_id == qry:
                    if target_id not in target_idL:
                        target_cls = target_id.split('_')
                        target_cls = target_cls[2]
                        target.append(int(target_cls))
                        target_idL.append(target_id)
                else:
                    q_retrival = 0
                    #print qry_id
                    # strt_file  = line_no
                    # end_file = len(txt_file)
                    # break

            found = 0
            reterived = 0
           #print qry
            #print target
            for target_cls in target:
                if target_cls == family+1:
                    found = found+1
                reterived = reterived+1
                recall = float(found)/len(relevant_seqs_clswise[family+1])
                precision = float(found)/reterived
                recall = int(np.ceil(recall*100))
                if precisionAt[recall] < precision:
                    precisionAt[recall] = precision

            indxs = range(0,100)
            indxs.reverse()
            for i in  indxs:
                if precisionAt[i] < precisionAt[i+1]:
                    precisionAt[i] = precisionAt[i+1]


            precisionClass  = precisionClass+precisionAt


            #print (precisionAt[0:101:10]) # for the query
            #wr.writerow(precisionAt[0:101:10]) # for the query
            file.close()

        precisionClass = precisionClass/n_queries # for the class

        #print (key)
        print (precisionClass[0:101:10])
        precision_avg += precisionClass[0:101:10]
        wr.writerow(precisionClass[0:101:10])
        wr.writerow([])




    precision_avg =precision_avg/len(cls_wise_queries.keys())
    wr.writerow(precision_avg)
    print(precision_avg)
    end = time.time()
    print("Time check", (end - strt))

        #wr.writerow(str(end - strt))

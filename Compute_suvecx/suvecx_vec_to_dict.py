"""
suvecx vectors in other format
"""
import numpy as np
import time
import pickle as pkl
import os
from collections import Counter
from sklearn.metrics.pairwise import cosine_distances



def load_fasttext_vec(file_):
    file = open(file_,'r')
    lines = file.readlines()
    vecs = []
    for i in range(0,len(lines)):
        # print(i)
        line = lines[i]
        # print(line)
        line = line.split(' ')
        line.remove('\n')
        # print(len(line))
        vecs.append([float(i) for i in line])
        # np.asarray(vecs)
    # print(vecs)
    # print(np.shape(vecs))
    return vecs


def extract_family_sizes(file):
    file_in = open(file,'r')
    i = 0
    cls_ids = []
    for line in file_in:
        # print(line)
        if line[0] == '_':
            line_temp = line.split(' ')
            line_temp = line_temp[0]
            # print(line_temp)
            line_temp = line_temp.split('__')
            # print(line_temp)
            line_temp = line_temp[2]
            cls_ids.append(int(line_temp))
            # print(line_temp)
            i += 1
            # if i == 5000:
            #     break
    # print(cls_ids)
    print(Counter(cls_ids).keys())
    family_sizes = Counter(cls_ids).values()
    print(sum(family_sizes))
    file_in.close()
    return family_sizes


# temp = open(os.path.join(Dir_path, DB_file + ".pkl"), 'rb')
# db_temp = pkl.load(temp)
# db = db_temp['sequence_vectors']
# db = np.asarray(db)
# db_family_sizes = db_temp['family_sizes']



if __name__ == "__main__":

    path_ftxt_fasta = r"C:\Users\Dhananjay\Desktop\Review_PlosOne\fasttext\100random_pair_exp\fasttext_files"
    path_ftxt_vec =  r"C:\Users\Dhananjay\Desktop\Review_PlosOne\fasttext\100random_pair_exp\fasttext_vecs"
    path_ftxt_vec_pkl =  r"C:\Users\Dhananjay\Desktop\Review_PlosOne\fasttext\100random_pair_exp\fasttext_vecs_pkl"



    for file_no in range(0,100):
        file_name1 = "T_"+str(file_no)+"_ftxt"
        file_name = "T"+str(file_no)+"_ftxt"
        ftxt_vec = file_name+"_vec"
        ftxt_pkl = ftxt_vec+'.pkl'
        ftxt_fasta = file_name1+".fasta"

        fasta_file = os.path.join(path_ftxt_fasta,ftxt_fasta)
        # print(fasta_file)
        family_sizes = extract_family_sizes(fasta_file)
        print(family_sizes)
        # # print(db[0])
        db = load_fasttext_vec(os.path.join(path_ftxt_vec,ftxt_vec)) # fasttext vectors
        pickle_out = open(os.path.join(path_ftxt_vec_pkl,ftxt_pkl),"wb")
        T0_fasttext_vec = {'sequence_vectors':db, 'family_sizes': family_sizes}
        pkl.dump(T0_fasttext_vec,pickle_out)
        pickle_out.close()



    # for i in range(0,1):
    #     DB_file = "DB0_ftxt_vec"
    #     Q_file = "Q0_ftxt_vec"
    #     main(Dir_path, DB_file, Q_file)

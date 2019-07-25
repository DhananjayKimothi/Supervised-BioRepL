"""
Notes : Using this script we can visuzalize
    (i) Query and DB sequence embeddings in 2D space (using t-SNE plots)


"""

from collections import defaultdict
from numpy import ones
from numpy import array, set_printoptions, linspace, save, load, count_nonzero
from sklearn.manifold import TSNE
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pickle as pkl
import os
import sys

path = r"C:\Users\Dhananjay\Desktop\Review_PlosOne\SUBMITTED\DATA_for_experiments\100random_pair_experiment"

### Database files .pkl
file1_fasta = "T1"
file1_bin = "T1_Bio"

### Query files .pkl
file2_fasta = "Q_DB1"
file2_bin = "Q_DB1_Bio"

def sequence_ids(file_path):
    IdsL = [] # L is for denoting list
    file = open(file_path,'r')
    for line in file:
        if line[0] == '>':
            line = line.replace('>', '')
            line = line.replace('\r\n','')
            IdsL.append(line)
    file.close()
    return IdsL


def LabelVec(lbl_str):
    cls_lbls = []

    # can be optimized here
    for lbl in lbl_str:
        x_str_split = lbl.split("_")
        cls_lbls.append(x_str_split[2])

    # count of samples in each class identified with its
    clslbl_count = defaultdict(int)
    for lbl in cls_lbls:
        if lbl not in clslbl_count:
            clslbl_count[lbl] = 1
        else:
            clslbl_count[lbl] += 1

    _no = 0
    int_lbls = []
    for cls in sorted(clslbl_count):
        _lbls = _no * ones(clslbl_count[str(cls)], dtype=int)
        int_lbls.extend(_lbls)
        _no += 1

    return int_lbls, clslbl_count

def tsne(X,x):

    """
    X: numpy array --- sequence embeddings
    x: labels --- 0: DB1, 1: DB2, 2: Q1, 3: Q2
    """
#########################################################################

    model = TSNE(n_components=2, random_state=0)
    set_printoptions(suppress=True)
    Z = model.fit_transform(X)

    C1 = x.count(0); C2 = x.count(1) # database sequences
    C11 = x.count(2); C12 = x.count(3) # query sequences

    strt = 0;end = C1
    Z1 = Z[strt:end]
    print("# database sequences class1", len(Z1))

    strt = end; end = strt+C2
    Z2 = Z[strt:end]
    print("# database sequences class2", len(Z2))

    strt = end; end = strt+C11
    Z11 = Z[strt:end]
    print("# query sequences class1",  len(Z11))

    strt = end; end = strt+C12
    Z12 = Z[strt:end]
    print("# query sequences class1", len(Z12))

    plt.scatter(Z1[:,0],Z1[:,1], marker = '*', color = 'C0', label = 'DB1')
    plt.scatter(Z2[:, 0], Z2[:, 1], marker='^', color='C1', label = 'DB2')

    plt.scatter(Z11[:, 0], Z11[:, 1], marker='*', color='C2', label='Q1')
    plt.scatter(Z12[:, 0], Z12[:, 1], marker='^', color='C3', label='Q2')

    plt.legend(loc = 2)

    plt.show()


def VecsandIds(path, file_fasta, file_bin):
    file_path = os.path.join(path, file_fasta+'.fasta')
    IdsL = sequence_ids(file_path)

    file_path = os.path.join(path, file_bin+'.pkl')
    in_file = open(file_path, 'rb')
    Database = pkl.load(in_file)
    in_file.close()
    Seqs_vecs = Database['sequence_vectors'] # list of lists

    return Seqs_vecs, IdsL



if __name__ == '__main__':

    # X1 numpy array of embeddings; IdsL1 are the labels (0,1) # for Database sequences, both classes
    X1, IdsL1 = VecsandIds(path,file1_fasta,file1_bin)
    # for query sequences
    X2, IdsL2 = VecsandIds(path,file2_fasta,file2_bin)

    IdsL1, xx = LabelVec(IdsL1);  IdsL2, xx = LabelVec(IdsL2)


    IdsL2 = [x+2 for x in IdsL2]
    IdsL = IdsL1 + IdsL2
    print("Total database and query samples",   len(IdsL))
    X = X1 + X2
    # print(len(X))
    tsne(X, IdsL)

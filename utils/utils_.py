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

def family_Sizes(IdsL):
    family_sizes = []
    for id in IdsL:
#         print(id)
        label = id.split('_')
#         print(label)
        label = int(label[2])
        if label not in family_dict.keys():
            family_dict[label] = 1
        else:
            family_dict[label] += 1
    for key in sorted(family_dict.keys()):
        family_sizes.append(family_dict[key])
    return(family_sizes)
from BioVecs import main
import time
for i in range(0,100):
    import os
    in_Dir_path = r"E:\SuperVec_Codes_25072018\Data\pairs_from1-200\pairs_fasta_files"
    out_Dir_path = r"E:\SuperVec_Codes_25072018\Data\pairs_from1-200\vectors"
    input_fasta_file = os.path.join(in_Dir_path, 'T_'+str(i)+'.fasta')
    output_BioVecs = os.path.join(out_Dir_path, 'T_'+str(i)+'Bio.pkl')
    main(input_fasta_file, output_BioVecs)

    in_Dir_path = r"E:\SuperVec_Codes_25072018\Data\pairs_from1-200\pairs_fasta_files"
    out_Dir_path = r"E:\SuperVec_Codes_25072018\Data\pairs_from1-200\vectors"
    input_fasta_file = os.path.join(in_Dir_path, 'Q_DB_'+str(i)+'.fasta')
    output_BioVecs = os.path.join(out_Dir_path, 'Q_DB_'+str(i) +'Bio.pkl')
    main(input_fasta_file, output_BioVecs)

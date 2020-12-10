'''
Returns selected parts of fasta file as well as SNP file
This allows for us to join our SNP file with a fasta file for easier analysis
'''

import pandas as pd
from Bio import SeqIO


def seq_snip(snp_file, seq_file):
    # read SNP file
    snp = pd.read_csv(snp_file, sep="\t")
    result = pd.DataFrame({"ENSEMBL_ID": [], "SNP": [], "SEQ": [], "SNP_SEQ": [], "SEQ_DB": [], "SNP_SEQ_DB": []})
    new_row = {}
    # read FASTA file for sequences
    fasta_dict = dict()
    #match sequence with info about sequence, seq_file is fasta file SNP file info about sequences
    with open(seq_file, "rt") as handle:
        # the following loops will add all the necessary information
        # to the output file using the .fa and .tsv input files:
        for (key, value) in SeqIO.FastaIO.SimpleFastaParser(handle):
            new_row["ENSEMBL_ID"] = key
            new_row["SEQ"] = value
            fasta_dict[key] = value
            snp_seq = list(fasta_dict[key])
            snps = snp.loc[(snp["ENSEMBL_ID"] == key)]
            for index, row in snps.iterrows():
                temp = snp_seq[row.SNP_REL_POS - 1]
                snp_seq[row.SNP_REL_POS - 1] = row.SNP_TO
                snp_seq_str = "".join(snp_seq)
                snp_seq[row.SNP_REL_POS - 1] = temp
                new_row["SNP"] = int(row["SNP"])
                new_row["SNP_FROM"] = row["SNP_FROM"]
                new_row["SNP_REL_POS"] = int(row["SNP_REL_POS"])
                new_row["SNP_TO"] = row["SNP_TO"]
                new_row["CLINSIG"] = row["CLINSIG"]
                new_row["BIOTYPE"] = row["BIOTYPE"]
                new_row["SNP_SEQ"] = snp_seq_str
                result = result.append(new_row, ignore_index=True)

    return result

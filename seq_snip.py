import pandas as pd
from Bio import SeqIO


def seq_snip(snp_file, seq_file):
    # read SNP file
    snp = pd.read_csv(snp_file, sep="\t")
    result = pd.DataFrame({"ENSEMBL_ID": [], "SNP": [], "SEQ": [], "SNP_SEQ": [], "SEQ_DB": [], "SNP_SEQ_DB": []})
    new_row = {}
    # read FASTA file for sequences
    fasta_dict = dict()
    with open(seq_file, "rt") as handle:
        for (key, value) in SeqIO.FastaIO.SimpleFastaParser(handle):
            new_row["ENSEMBL_ID"] = key
            new_row["SEQ"] = value
            fasta_dict[key] = value
            snp_seq = list(fasta_dict[key])
            snps = snp.loc[(snp["ENSEMBL_ID"] == key)]
            for index, row in snps.iterrows():
                # print(fasta_dict[key])
                temp = snp_seq[row.SNP_REL_POS - 1]
                snp_seq[row.SNP_REL_POS - 1] = row.SNP_TO
                snp_seq_str = "".join(snp_seq)
                # print(snp_seq_str)
                snp_seq[row.SNP_REL_POS - 1] = temp
                new_row["SNP"] = int(row["SNP"])
                new_row["SNP_SEQ"] = snp_seq_str
                result = result.append(new_row, ignore_index=True)

    return result

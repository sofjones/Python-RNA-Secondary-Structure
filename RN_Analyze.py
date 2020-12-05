from energy_min import energy_min
import seq_snip
from subprocess import Popen
import sys
import RNA

def fold(sequence):
    fold = Nussinov.nussinov(sequence)
    dot_bracket = list("".join("%s" % "." for i in range(len(row.SEQ))))
    for pairing in fold:
        # print(pairing)
        dot_bracket[min(pairing) - 1] = "("
        dot_bracket[max(pairing) - 1] = ")"
    return dot_bracket


def main(argv):
    snp_seq_data = seq_snip.seq_snip(argv[0], argv[1])
    for index, row in snp_seq_data.iterrows():
        # fold = energy_min(list(row.SEQ.lower()))
        # snp_fold = energy_min((list(row.SNP_SEQ.lower())))
        dot_bracket = energy_min(list(row.SEQ))
        snp_dot_bracket = energy_min((list(row.SNP_SEQ)))
        dist = RNA.bp_distance(str(dot_bracket),str(snp_dot_bracket))
				# for pairing in fold:
        #     # print(pairing)
        #     dot_bracket[min(pairing)-1] = "("
        #     dot_bracket[max(pairing)-1] = ")"
        # for pairing in snp_fold:
        #     # print(pairing)
        #     snp_dot_bracket[min(pairing)-1] = "("
        #     snp_dot_bracket[max(pairing)-1] = ")"
        snp_seq_data.loc[index, "SEQ_DB"] = "".join(dot_bracket)
        snp_seq_data.loc[index, "SNP_SEQ_DB"] = "".join(snp_dot_bracket)
        snp_seq_data.loc[index, "SEQ_SB"] = "".join("\n {d}".format(d=dist))
    snp_seq_data.to_csv("SEQ_DB.tsv", sep="\t")
    if len(argv) == 3 and argv[2] == "p":
        command = ["python", "RNAvisual.py", "SEQ_DB.tsv"]
        Popen(command)


if __name__ == "__main__":
    main(sys.argv[1:])

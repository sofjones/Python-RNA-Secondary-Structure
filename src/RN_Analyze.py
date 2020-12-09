'''
RN_Analyze allows us to create dot bracket structures for every sequence
Reads through every snip created in seq_snip and prints corresponding structs to .tsv file
'''

# energy minimization algorithm called upon for each sequence
from energy_min import energy_min
#extracts information from data files
import seq_snip
# Popen is used in case the user includes the argument "p"
# in order to automatically run RNAvisual.py
from subprocess import Popen
import sys
#comment out to remove base pair distance
import RNA

def main(argv):
    # seq_snip will organize the .fa and .tsv file data for the program:
    snp_seq_data = seq_snip.seq_snip(argv[0], argv[1])
	
    for index, row in snp_seq_data.iterrows():

        # create dot bracket structures and add to .tsv file:
        dot_bracket = energy_min(list(row.SEQ)) # dot-bracket for the original RNA
        snp_dot_bracket = energy_min(list(row.SNP_SEQ)) # dot-bracket for the SNP
        snp_seq_data.loc[index, "SEQ_DB"] = "".join(dot_bracket)
        snp_seq_data.loc[index, "SNP_SEQ_DB"] = "".join(snp_dot_bracket)			

        # add column for distance (our "variability" score)
        dist = RNA.bp_distance(str(dot_bracket),str(snp_dot_bracket))
        snp_seq_data.loc[index, "SEQ_SB"] = "".join("\n {d}".format(d=dist))
    # we output to the file SEQ_DB.tsv				     
    snp_seq_data.to_csv("SEQ_DB.tsv", sep="\t")
    # if the user includes the argument "p"
    # then visualization is triggered automatically:
    if len(argv) == 3 and argv[2] == "p":
        command = ["python", "RNAvisual.py", "SEQ_DB.tsv"]
        Popen(command)

#start here
if __name__ == "__main__":
    main(sys.argv[1:])

# RNA Secondary Structure Prediction

Using dataset of Non Coding SNPs predict secondary structure with a variation of Nussinovs algorithm, `energy_min.py`
`RN_Analyze.py` is run first and will use `energy_min.py` to return a .tsv file `SEQ_DB.tsv` with sequences and their corresponding dot bracket structures. Base pair distance is calculated using `bp_distance` from the RNA package taken from ViennaRNA.
Visuals are created with `RNAvisual.py`.
## To Run
First run: `python RN_Analyze.py SNP.tsv sequences.fa`
To run visualization: `python RNAvisual.py SEQ_DB.tsv`

## Installations Required

### RN_Analyze
Reading fasta file: `pip install biopython` 
Bp distance: this requires miniconda 
1. `conda create -n viennarna -c bioconda viennarna` 
2. `conda activate viennarna`

Alternatively follow steps provided by viennarna:
			https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/install.html
### RNAvisual
`pip install dash`
`pip install dash_bio`
`pip install pandas`


 ## Dataset

Filtering all clinically significant SNPs

 Biotypes
 - mitochondrial tRNA
 - micro RNA
 - Small nucleolar RNA

Clinical significance

 - Benign, likely benign, uncertain, likely pathogenic, pathogenic

SNP location and substitution

Compare remaining SNPs to all ncRNA to select appropriate genes

## energy\_min.py

Follows nussinovs algorithm with the following additions

 

 1. Minimum loop
	
	This means that pairs must have some distance between them
	Ex. min loop of 2 would reject ..(.).. And accept .(..)..

2. Stacked base pairs

	((.......)) preferred over (.(...)...)

3. Score Pairing

	GC preferred over AU and GU

## RN\_Analyze

run this first with

`python RN\_Analyze SNP.tsv sequences.fa`

outputs `SEQ\_DB.tsv`

uses energy min to create dot structures for RNA sequences.

## RNAVisual

Visual representations using `dash`
to run this you must first run python 	`RN\_Analyze.py SNP.tsv sequences.fa`

after this run
`python RNAvisual.py SEQ\_DB.tsv`


`dash` will run on http://127.0.0.1:8050/


""" head eQTL_orgQTLs_YRI.py -n 14
This finds all possible QTL variants in the AFGR eQTL dataset that have a STARR-seq variant in their open chromatin region

Inputs:
    -QTL filename and path
    -BIRD filename and path
    -Output name for the QTL file
    -maximum distance between a QTL and BIRD SNP
    -maximum p value for the QTL SNPs

Outputs:
    -QTL file with indices for each connection to a BIRD file
    -BIRD file where each variant has a list of QTL indices it is associated with

"""

# %% Import packages
import tabix
import argparse
import csv


# %% Define inputs
parser = argparse.ArgumentParser(
                  prog="eQTL_orgQTLs_YRI",
                  description="This program finds the connections between QTLs and STARR-seq effects",
                  epilog=""
)

parser.add_argument("-Q", "--QTL_filename", help="The path to the input QTL file")
parser.add_argument("-S", "--SS_filename",help="The path to the input BIRD file")
parser.add_argument("-q", "--output_QTL_name", default="QTL.tsv", help="The path to the output QTL file")
parser.add_argument("-p", "--QTL_pval", type=float, default=1, help="The maximum QTL p-value")
parser.add_argument("-d", "--distance", type=int, default=100000, help="The maximum distance between a QTL and a STARR-seq variant")
args = parser.parse_args()

#QTL_filename = sys.argv[1]
#SS_filename = sys.argv[2]
#output_QTL_name = sys.argv[3]
#distance = int(sys.argv[4])
#QTL_pval = float(sys.argv[5])


# %% Open files
SS_file = csv.reader(open(args.SS_filename), delimiter = "\t")
tb = tabix.open(args.QTL_filename)

SS_data = []
for line in SS_file:
    SS_data.append(line)
    SS_data[-1].append([])



total_found = 0

# %% Find the variants that are close together
output_cQTL = open(args.output_QTL_name, 'w')

SS_all = []
QTL_all = []

for line in SS_data:

    # Find the variant
    variant = line[0].split(":")
    chrom = variant[0][3:]
    end_val = int(variant[1])+args.distance
    start_val = int(variant[1])-args.distance

    # Find the corresponding lines
    records = tb.query(chrom, start_val, end_val)
    total_records = 0


    for r in records:
        total_records += 1

        current_variant = []
        current_variant.append("chr" + r[0]) # chromosome
        current_variant.append(r[1]) # position
        current_variant.append(r[3]) # ref
        current_variant.append(r[4]) # alt
        current_variant.append(r[10]) # beta
        current_variant.append(r[9]) # p
        current_variant.append(r[13]) # chrom_start
        current_variant.append(r[14]) # chrom_end
        current_variant.append(r[1]) # var_start
        current_variant.append(r[2]) # var_end
        current_variant.append(r[5]) # effect_af_eqtl
        current_variant.append(r[6]) # variant
        current_variant.append(r[7]) # feature
        current_variant.append(r[8]) # dist_fastql
        current_variant.append(r[11]) # se
        current_variant.append(r[12]) # p_hwe
        current_variant.append(r[15]) # gene_symbol
        current_variant.append(r[16]) # gene_type

        if float(r[9]) < args.QTL_pval:
            output_cQTL.write("\t".join(current_variant)+"\n")





QTL_header_file = open(args.output_QTL_name+".header", 'w')
QTL_header = ["chromosome", "position", "ref", "alt", "beta", "p", "chrom_start", "chrom_end", "start", "end",
              "effect_af_eqtl", "variant", "feature", "log10p", "se", "qstat",
              "df", "p_het", "h_hwe", "gene_symbol", "gene_type"]
QTL_header_file.write("\t".join(QTL_header) + "\n")


output_cQTL.close()


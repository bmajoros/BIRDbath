""" head altius_betas.py -n 14
This finds the connections between the QTL SNPs and the BIRD SNPs with the gzipped altius file using tabix

Inputs:
    -QTL filename and path
    -BIRD filename and path
    -Output name for the BIRD file
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
import statistics as sts


# %% Define inputs
parser = argparse.ArgumentParser(
                  prog="HF_betas",
                  description="This program finds the connections between QTLs and STARR-seq effects",
                  epilog=""
)

parser.add_argument("-Q", "--QTL_filename", help="The path to the input QTL file")
parser.add_argument("-S", "--SS_filename",help="The path to the input BIRD file")
parser.add_argument("-q", "--output_QTL_name", default="QTL.tsv", help="The path to the output QTL file")
parser.add_argument("-s", "--output_bird_name", default = "BIRD.tsv", help="The path to the output BIRD file")
parser.add_argument("-p", "--QTL_pval", type=float, default=1, help="The maximum QTL p-value")
parser.add_argument("-d", "--dist", type=int, default=100000, help="The maximum distance between a QTL and a STARR-seq variant")
args = parser.parse_args()

#DNase_filename = sys.argv[1]
#SS_filename = sys.argv[2]
#output_bird_name = sys.argv[3]
#output_cQTL_name = sys.argv[4]
#distance = int(sys.argv[5])
#pval_thresh = float(sys.argv[6])


# %% Open files
SS_file = csv.reader(open(args.SS_filename), delimiter = "\t")
tb = tabix.open(args.QTL_filename)


output_bird = open(args.output_bird_name, 'w')
output_cQTL = open(args.output_QTL_name, 'w')


header_values_list = ["chr", "SNP", "ref", "alt", "beta", "p"]
header_line = "\t".join(header_values_list)+"\n"
output_cQTL.write(header_line)


total_found = 0

# %% Find the variants that are close together

SS_all = []
QTL_all = []

for line in SS_file:
    
    # Find the variant
    variant = line[0].split(":")
    chrom = variant[0][3:]
    end_val = int(variant[1])+args.dist
    start_val = int(variant[1])-args.dist
    
    # Find the corresponding lines
    records = tb.query("chr"+chrom, start_val, end_val)
    total_records = 0
    
    for r in records:
      if float(r[3]) < args.QTL_pval:
        total_records += 1
        
        current_variant = []
        current_variant.append(r[0]) #chr
        current_variant.append(r[1]) #variant
        current_variant.append(r[4]) #ref
        current_variant.append(r[5]) #alt
        current_variant.append(r[2]) #beta
        current_variant.append(r[3]) #p
        current_variant.append(total_found)
        
        QTL_all.append(current_variant)
    
    if total_records > 0:
        total_found += 1
        SS_line = line.copy()
        SS_line.append([])
        SS_all.append(SS_line)
    



# %% Reorganize for outputting

SS_out = []
QTL_out = []

QTL_used = []

total_printed = 0

for var in QTL_all:
    
    QTL_VAR = var[:-1].copy()
    SS_ind = var[-1]
    
    if QTL_VAR[0]+":"+QTL_VAR[1] not in QTL_used:
        
        total_printed += 1
        QTL_VAR.append(str(total_printed))
        QTL_used.append(QTL_VAR[0]+":"+QTL_VAR[1])
        QTL_out.append(QTL_VAR)
        
        SS_all[SS_ind][-1].append(total_printed)
    
    else:
        
        QTL_ind = QTL_used.index(QTL_VAR[0]+":"+QTL_VAR[1]) + 1
        
        SS_all[SS_ind][-1].append(QTL_ind)
        
    
#%% Save the files

# Save the QTL file
for line in QTL_out:
    output_cQTL.write("\t".join(line)+"\n")

# Save the BIRD file
for line in SS_all:
    if len(line[-1]) != 0:
        variant = line[0].split(":")
        write_line = [variant[0], variant[1]] + line[1:-1] + [variant[2], variant[3]]
        output_bird.write("\t".join(line))
        output_bird.write("\t"+str(line[-1])+"\n")

# Close the output file
output_bird.close()
output_cQTL.close()




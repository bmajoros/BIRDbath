""" head hQTL_betas_notabix.py -n 14
This finds the connections between the QTL SNPs and the BIRD SNPs with the unzipped hQTL file with only one line for each variant

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

#%% Import packages
import argparse
import csv
import gzip


#%% Define inputs
parser = argparse.ArgumentParser(
                  prog="hQTL_betas_notabix",
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

#QTL_filename = sys.argv[1]
#SS_filename = sys.argv[2]
#output_bird_name = sys.argv[3]
#output_cQTL_name = sys.argv[4]
#pval_thresh = float(sys.argv[5])
#dist = int(sys.argv[6])


#%% Load BIRD file

# Open the file
SS_file = csv.reader(open(args.SS_filename), delimiter = "\t")
# Initialize permenant variables
SS_vals = []
SS_pos = []
# Set up to divide by chromosome
for i in range(23):
    SS_vals.append([])
    SS_pos.append([])

# Run through the file and save each line
for line in SS_file:
    # Find the chromosome and position
    SNP_id = line[0]
    variant = line[0].split(":")
    if "X" in SNP_id: chrom=23
    else: chrom = int(variant[0][3:])
    pos = int(variant[1])
    # Save the information    
    save_line = list(line)
    save_line.append([])
    if len(SS_vals[chrom-1]) == 0:
        SS_pos[chrom-1].append(pos)
        SS_vals[chrom-1].append(save_line)
    elif SS_pos[chrom-1][-1] < pos:
        SS_pos[chrom-1].append(pos)
        SS_vals[chrom-1].append(save_line)
    else:
        insert_index = 0
        while SS_pos[chrom-1][insert_index]<pos:
            insert_index += 1
        SS_pos[chrom-1].insert(insert_index, pos)
        SS_vals[chrom-1].insert(insert_index, save_line)


#%% Search QTL file
total_saved = 0
QTL_file = gzip.open(args.QTL_filename)
old_chrom = None
QTL_saved = []
for raw_line in QTL_file:
  line = str(raw_line)[2:-3].split("\\t")
  if "variant_id" not in line:
    # Find the chromosome
    if "X" in line[0]: chrom=23
    else: chrom = int(line[0])
    pos = int(line[1])
    # Check if it changed
    if chrom != old_chrom:
        start_ind = 0
    # Find the region of valid open chromatin
    valid_region = []
    valid_region.append(pos-args.dist)
    valid_region.append(pos+args.dist)
    # Find the corresponding BIRD SNPs
    pval = float(line[4])
    if pval < args.QTL_pval:
        corresponding_BIRD = []
        while start_ind > 0 and SS_pos[chrom-1][start_ind] > valid_region[0]:
            start_ind -= 1
        while start_ind < len(SS_pos[chrom-1]) and SS_pos[chrom-1][start_ind] < valid_region[0]:
            start_ind += 1
        end_ind = start_ind
        while end_ind < len(SS_pos[chrom-1]) and SS_pos[chrom-1][end_ind] <= valid_region[1]:
            end_ind += 1
        
        if end_ind != start_ind:
            total_saved += 1
            for i in range(start_ind, end_ind):
                SS_vals[chrom-1][i][-1].append(total_saved)
            line.append(str(total_saved))
            QTL_saved.append(line)


#%% Save the files

# Save the QTL file
QTL_save_file = open(args.output_QTL_name,'w')
for line in QTL_saved:
    save_line = []
    save_line.append("chr" + line[0])  # chromosome
    save_line.append(line[1])  # position
    save_line.append(line[8])  # ref
    save_line.append(line[9])  # alt
    save_line.append(line[7]) # beta
    save_line.append(line[4])  # p
    save_line.append(line[2])  # epiQTL_rsID
    save_line.append(line[3])  # Mapped_Gene
    save_line.append(line[5])  # H3K27ac_pval
    save_line.append(line[6])  # H2K4mel_pval
    save_line.append(line[10])  # Total_Reads
    save_line.append(line[11]) # Reference_Reads
    save_line.append(line[12]) # Alternate_Reads
    save_line.append(line[13]) # H3K27ac_Replication_pval
    save_line.append(line[-1]) # save number
    QTL_save_file.write("\t".join(save_line)+"\n")

QTL_header_file = open(args.output_QTL_name + ".header", 'w')
QTL_header_file.write("chromosome\tposition\tref\talt\tbeta\tp\epiQTL_rsID\tMapped_Gene\tH3K27ac_pval\tH2K4mel_pval\tTotal_Reads\tReference_Reads\tAlternate_Reads\tH3K27ac_Replication_pval\n")

# Save the BIRD file
BIRD_save_file = open(args.output_bird_name,'w')
for chrom in SS_vals:
    for line in chrom:
        if len(line[-1]) != 0:
            variant = line[0].split(":")
            write_line = [variant[0], variant[1]] + line[1:-1] + [variant[2], variant[3]]
            BIRD_save_file.write("\t".join(write_line))
            BIRD_save_file.write("\t"+str(line[-1])+"\n")

QTL_save_file.close()
QTL_header_file.close()
BIRD_save_file.close()

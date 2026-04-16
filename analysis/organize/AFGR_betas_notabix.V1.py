""" head AFGR_betas_notabix.py -n 14
This finds the connections between the QTL SNPs and the BIRD SNPs with the unzipped AFGR file with only one line for each variant

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


#%% Define inputs
parser = argparse.ArgumentParser(
                  prog="AFGR_betas_notabix",
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
    #at_index = SNP_id.index("@")
    variant = SNP_id.split(":")
    if "X" in SNP_id: chrom=23
    else: chrom = int(variant[0][3:])
    pos = int(variant[1])
    # Save the information    
    #save_line = list(line)
    save_line = variant[:2] + list(line[1:]) + variant[2:]
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
QTL_file = csv.reader(open(args.QTL_filename), delimiter = "\t")
old_chrom = None
QTL_saved = []
for line in QTL_file:
  if line[0] != "chr":
    if "chr" not in line[0]: line[0] = "chr" + line[0]
    # Find the chromosome
    chrom = int(line[0][3:])
    pos = int(line[1])
    # Check if it changed
    if chrom != old_chrom:
        start_ind = 0
    # Find the region of valid open chromatin
    valid_region = []
    #valid_region.append(pos-int(line[9])-dist)
    #valid_region.append(pos-int(line[10])+dist)
    valid_region.append(pos-args.dist)
    valid_region.append(pos+args.dist)
    # Find the p-value
    pval = float(line[11])
    if pval < args.QTL_pval:
        # Find the corresponding BIRD SNPs
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
    save_line.append(line[0])  # chromosome
    save_line.append(line[1])  # position
    save_line.append(line[3])  # ref
    save_line.append(line[4])  # alt
    save_line.append(line[12]) # beta
    save_line.append(line[11]) # p
    save_line.append(line[9])  # start
    save_line.append(line[10]) # end
    save_line.append(line[2])  # position 2
    save_line.append(line[5])  # variant
    save_line.append(line[6])  # effect_af_eqtl
    save_line.append(line[7])  # p_hwe
    save_line.append(line[8])  # feature
    save_line.append(line[13]) # se
    save_line.append(line[-1]) # save number
    QTL_save_file.write("\t".join(save_line)+"\n")

QTL_header_file = open(args.output_QTL_name + ".header", 'w')
QTL_header_file.write("chromosome\tposition\tref\talt\tbeta\tp\tstart\tend\tpos2\tvariant\teffect_af_eqtl\tp_hwe\tfeature\tse\n")

# Save the BIRD file
BIRD_save_file = open(args.output_bird_name,'w')
for chrom in SS_vals:
    for line in chrom:
        if len(line[-1]) != 0:
            BIRD_save_file.write("\t".join(line[:-1]))
            BIRD_save_file.write("\t"+str(line[-1])+"\n")

QTL_save_file.close()
QTL_header_file.close()
BIRD_save_file.close()

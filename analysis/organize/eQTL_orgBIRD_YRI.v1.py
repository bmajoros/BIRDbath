""" head eQTL_orgBIRD_YRI.py -n 14
This finds the connections between the QTL SNPs and the BIRD SNPs with the rganized QTL file with only one line for each variant

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
import argparse
import csv


# %% Define inputs
parser = argparse.ArgumentParser(
                  prog="eQTL_orgBIRD_YRI",
                  description="This program finds the connections between QTLs and STARR-seq effects",
                  epilog=""
)

parser.add_argument("-Q", "--QTL_filename", help="The path to the input QTL file")
parser.add_argument("-S", "--SS_filename",help="The path to the input BIRD file")
parser.add_argument("-q", "--output_QTL_name", default="QTL.tsv", help="The path to the output QTL file")
parser.add_argument("-s", "--output_bird_name", default = "BIRD.tsv", help="The path to the output BIRD file")
parser.add_argument("-p", "--QTL_pval", type=float, default=1, help="The maximum QTL p-value")
parser.add_argument("-d", "--distance", type=int, default=100000, help="The maximum distance between a QTL and a STARR-seq variant")
args = parser.parse_args()

#QTL_filename = sys.argv[1]
#SS_filename = sys.argv[2]
#output_bird_name = sys.argv[3]
#output_QTL_name = sys.argv[4]
#distance = int(sys.argv[5])
#QTL_pval = float(sys.argv[6])


# %% Open files
SS_file = csv.reader(open(args.SS_filename), delimiter = "\t")


# Sort the file
# Set up the output structure
SS_data = []
SS_pos = []
for i in range(22):
    SS_data.append([])
    SS_pos.append([])

# Loop through the variants in the input file
for line in SS_file:

    SNP_id = line[0]
    variant = line[0].split(":")
    use_line = line.copy()
    use_line.append([])

    # Find the chromosome
    current_chromosome = int(variant[0][3:]) - 1

    # Find the position
    current_position = int(variant[1])

    # Find the index of the position in the output list
    # Put it in the right place
    if len(SS_data[current_chromosome]) == 0:
        SS_pos[current_chromosome].append(current_position)
        SS_data[current_chromosome].append(use_line)
    elif SS_pos[current_chromosome][-1] < current_position:
        SS_pos[current_chromosome].append(current_position)
        SS_data[current_chromosome].append(use_line)
    else:
        insert_index = 0
        while SS_pos[current_chromosome][insert_index]<current_position:
            insert_index += 1
        SS_pos[current_chromosome].insert(insert_index, current_position)
        SS_data[current_chromosome].insert(insert_index, use_line)

output_bird = open(args.output_bird_name, 'w')
output_cQTL = open(args.output_QTL_name, 'w')




# %% Find the variants that are close together


total_saved = 0
QTL_file = csv.reader(open(args.QTL_filename), delimiter = "\t")
old_chrom = None
for line in QTL_file:
    # Find the chromosome
    chrom = int(line[0][3:])
    pos = int(line[1])
    sig_val = float(line[5])

    if sig_val < args.QTL_pval:

        # Check if it changed
        if chrom != old_chrom:
            start_ind = 0
        # Find the region of valid open chromatin
        valid_region = []
        valid_region.append(int(line[6])-args.distance)
        valid_region.append(int(line[7])+args.distance)
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
                SS_data[chrom-1][i][-1].append(total_saved)
            write_line = line.copy()
            #write_line[0] = line[0] + ":" + line[1]
            #write_line.pop(1)
            write_line.append(str(total_saved))
            output_cQTL.write("\t".join(write_line)+"\n")


QTL_header_file = open(args.output_QTL_name+".header", 'w')
QTL_header = ["chromosome", "position", "ref", "alt", "beta", "p", "start", "end",
              "effect_af_eqtl", "variant", "feature", "log10p", "se", "qstat",
              "df", "p_het", "h_hwe", "gene_symbol", "gene_type"]
QTL_header_file.write("\t".join(QTL_header) + "\n")




#%% Save the files

# Save the BIRD file
for chrom_list in SS_data:
  for line in chrom_list:
    if len(line[-1]) != 0:
        variant = line[0].split(":")
        write_line = [variant[0], variant[1]] + line[1:-1] + [variant[2], variant[3]] 
        output_bird.write("\t".join(write_line))
        output_bird.write("\t"+str(line[-1])+"\n")

# Close the output file
output_bird.close()
output_cQTL.close()


""" head process_LD.py -n 14
This script finds the counts of all pairs of alleles which can then be calculated into LD statistics

Inputs:
    -The name and location of the file with the pairs of SNPs
    -The 3 letter code for each population to check
        -groups are separated with "^"
        -populations are separated with "," within a group
    -The name of the output file without the file extension
    -The integer number of the chromosome to process

Outputs:
    -remaining_*.txt" : pairs of SNPs where at least one of the two coordinates was not in the Thousand Genomes refrence VCF, as well as the indicator for which SNP was missing
        -0 indicates both were missing
        -1 indicates the first one was missing
        -2 indicates the second one was missing
    -phasing_matrix_*.txt : the file containing the two SNPs, all pairs of alleles that are represented by these SNPs, and the counts in each group of population listed in the order specified in the input

"""

# %% Import packages
import gzip
import csv
import sys
import sample_names
from datetime import datetime

start_time = datetime.now()

# %% Define inputs
inputs_filename = sys.argv[1]
population_all = sys.argv[2]
outname = sys.argv[3]
chrom = int(sys.argv[4])

# %% Find the sample names that will be used for the specified population

# Separate each group of populations
populations_list = population_all.split("^")

# Define the list of sample names
samples_list = []
samples_list_separated = []

# Check each group of populations in the given list of groups
for pop in range(len(populations_list)):

    # Create a new sublist for the associated sample names
    samples_list_separated.append([])
    # Find the string for this grouping of samples
    population = populations_list[pop]


    # Check each population name that is included in the 1000 genomes VCF
    # Save the sample names associated iwth that population to the current sublist

    if "ACB" in population:
        samples_list += sample_names.ACB()
        samples_list_separated[-1] += sample_names.ACB()

    if "ASW" in population:
        samples_list += sample_names.ASW()
        samples_list_separated[-1] += sample_names.ASW()

    if "BEB" in population:
        samples_list += sample_names.BEB()
        samples_list_separated[-1] += sample_names.BEB()

    if "CDX" in population:
        samples_list += sample_names.CDX()
        samples_list_separated[-1] += sample_names.CDX()

    if "CHB" in population:
        samples_list += sample_names.CHB()
        samples_list_separated[-1] += sample_names.CHB()

    if "CHS" in population:
        samples_list += sample_names.CHS()
        samples_list_separated[-1] += sample_names.CHS()

    if "CLM" in population:
        samples_list += sample_names.CLM()
        samples_list_separated[-1] += sample_names.CLM()

    if "ESN" in population:
        samples_list += sample_names.ESN()
        samples_list_separated[-1] += sample_names.ESN()

    if "FIN" in population:
        samples_list += sample_names.FIN()
        samples_list_separated[-1] += sample_names.FIN()

    if "GBR" in population:
        samples_list += sample_names.GBR()
        samples_list_separated[-1] += sample_names.GBR()

    if "GIH" in population:
        samples_list += sample_names.GIH()
        samples_list_separated[-1] += sample_names.GIH()

    if "GWD" in population:
        samples_list += sample_names.GWD()
        samples_list_separated[-1] += sample_names.GWD()

    if "IBS" in population:
        samples_list += sample_names.IBS()
        samples_list_separated[-1] += sample_names.IBS()

    if "ITU" in population:
        samples_list += sample_names.ITU()
        samples_list_separated[-1] += sample_names.ITU()

    if "JPT" in population:
        samples_list += sample_names.JPT()
        samples_list_separated[-1] += sample_names.JPT()

    if "KHV" in population:
        samples_list += sample_names.KHV()
        samples_list_separated[-1] += sample_names.KHV()

    if "LWK" in population:
        samples_list += sample_names.LWK()
        samples_list_separated[-1] += sample_names.LWK()

    if "MSL" in population:
        samples_list += sample_names.MSL()
        samples_list_separated[-1] += sample_names.MSL()

    if "MXL" in population:
        samples_list += sample_names.MXL()
        samples_list_separated[-1] += sample_names.MXL()

    if "PEL" in population:
        samples_list += sample_names.PEL()
        samples_list_separated[-1] += sample_names.PEL()

    if "PJL" in population:
        samples_list += sample_names.PJL()
        samples_list_separated[-1] += sample_names.PJL()

    if "PUR" in population:
        samples_list += sample_names.PUR()
        samples_list_separated[-1] += sample_names.PUR()

    if "STU" in population:
        samples_list += sample_names.STU()
        samples_list_separated[-1] += sample_names.STU()

    if "TSI" in population:
        samples_list += sample_names.TSI()
        samples_list_separated[-1] += sample_names.TSI()

    if "YRI" in population:
        samples_list += sample_names.YRI()
        samples_list_separated[-1] += sample_names.YRI()


# %% Process the VCF file

# This script doesn't loop through multiple chromosomes right now but it can be altered to do that here
# It takes about an hour to process chromosome 1 so I don't have plans to add that in but I could
# All chromosomes have to be checked separately because they are in different VCF files
#for chrom_i in range(1):

loop_start = datetime.now()

#chrom = chrom_i + 1

print("Starting Chromosome "+str(chrom))


# Find the VCF and pull back the relevant info from it

# Open the vcf
vcf_name = "/data/reddylab/kd259/ThousandGenomes/CCDG_14151_B01_GRM_WGS_2020-08-05_chr"+str(chrom)+".filtered.shapeit2-duohmm-phased.vcf.gz"  # input("What is the vcf filename? ")
vcf_file = gzip.open(vcf_name)

# Find the header line
target_chr = "#" # This character indicates a header line
next_line = "" # Initialize the next line
while target_chr == "#": # While the previou line was a header line
    prev_line = next_line # Save the previous line
    next_line = vcf_file.readline() # Read the next line in the file
    target_chr = str(next_line)[2] # Save the target character in this new line
# Delete the target character because it won't be used again
del target_chr


# %% Define the important columns in the file

# Find the locations of the important information in each line
# The first 5 columns have general info about the SNP that are important
# Then you need the column numbers for each sample in the populations that were specified in the input

# Save the header of the file as a list
header_list = str(prev_line).split("\\t") # Split the header line into a list

important_locations = [0,1,2,3,4] # Save the chromosome and position, and the alt and ref sequences
condensed_header_locations = [0,1,2,3,4]

for sample in samples_list: # For each sample in the pool
    try:
        # Put the column number of that sample into the list of columns if it is in the file
        important_locations.append(header_list.index(sample))
        # Save a shorter list with just the samples included in the important locations list
        condensed_header_locations.append(sample)
    except:
        pass

# Create a second list of columns that is sparated by the population group
important_locations_separated = []
for pop in range(len(populations_list)): # For each sample in the pool
    important_locations_separated.append([0,1,2,3,4]) # Start with the first 5 columns
    # Then add the sample columns for this group the same way as before
    for sample in samples_list_separated[pop]:
        try:
            important_locations_separated[-1].append(condensed_header_locations.index(sample)) # Put the column number of that sample into the list of columns
        except:
            pass

# Check that there are samples that were found for the specified populations
if len(important_locations) > 5:

    # Pull the pairs and chromosomal coordinates from the input file for the specified chromosome
    inputs_file = csv.reader(open(inputs_filename,'r'),delimiter='\t') # open file
    pairs = [] # initialize list of pairs
    locations = [] # initialize list of coordinates
    for line in inputs_file: # read the input file
        if "chr" + str(chrom) + ":" in line[0]: # use only pairs in the current chromosome
            pairs.append(line)

            # save each coordinate only once from the first column
            var1 = line[0].split(":")
            loc = int(var1[1])
            if loc in locations:
                pass
            else:
                locations.append(loc)
            # save each coordinate only once from the second column
            var2 = line[1].split(":")
            loc = int(var2[1])
            if loc not in locations:
                locations.append(loc)
    del loc
    locations.sort()

    print("Searching for "+str(len(locations))+" lines in VCF")



    # Pull back the relevant lines


    # Make a variable to hold the phasing data
    phasing_data = []
    phasing_data_locations = []

    # Check the current line
    current_line_list = str(next_line).split("\\t")
    temp_list = [] # Make a temporary list for the curent line
    total_processed = 0
    if int(current_line_list[1]) in locations: # if the current line is one of the bird positions
        total_processed += 1
        temp_list = []
        for loc in important_locations: # For each important column
            temp_list.append(current_line_list[loc]) # put the value from that column into the list for the row
        phasing_data.append(temp_list) # add the current list to the output lst
        phasing_data_locations.append(int(current_line_list[1]))

    # Skim through the file
    for next_line in vcf_file:
        # If it is a variant
        current_line_list = str(next_line).split("\\t")
        if int(current_line_list[1]) in locations: # if the current line is one of the bird or qtl coordiantes
            temp_list = []
            for loc in important_locations: # For each important column
                temp_list.append(current_line_list[loc]) # put the value from that column into the list for the row
            phasing_data.append(temp_list) # add the current list to the output lst
            phasing_data_locations.append(int(current_line_list[1]))



    # Calculate the phasing multipliers
    multipliers = []
    unused_pairs = []
    unused_pairs_reasons = []

    print("\nFinding values\n")
    # Loop through each of the pairs
    not_in_phasing_data_locations = 0
    #reasons = []
    for pair in pairs:
        pair_coord_1 = int(pair[0].split(":")[1])
        pair_coord_2 = int(pair[1].split(":")[1])
        if pair_coord_1 in phasing_data_locations and pair_coord_2 in phasing_data_locations:
            ind1 = phasing_data_locations.index(pair_coord_1)
            ind2 = phasing_data_locations.index(pair_coord_2)

            # Find all the phased lines
            phased_lines = [phasing_data[ind1], phasing_data[ind2]]
            whole_matrix = []

            # Define the ID values
            ID_values = []
            ID_values.append(phasing_data[ind1][3]+"_"+phasing_data[ind2][3])
            ID_values.append(phasing_data[ind1][4]+"_"+phasing_data[ind2][3])
            ID_values.append(phasing_data[ind1][3]+"_"+phasing_data[ind2][4])
            ID_values.append(phasing_data[ind1][4]+"_"+phasing_data[ind2][4])

            for pop in range(len(populations_list)):
                # matrix holds the counts of each combination of alleles in the file
                # matrix[0] = both alleles are 0
                # matrix[1] = the first allele is 1, the second allele is 0
                # matrix[2] = the first allele is 0, the second allele is 1
                # matrix[3] = both alleles are 1
                matrix = [0,0,0,0]
                # Save every double heterozygous line
                for value_int in important_locations_separated[pop]:
                    if "|" in phased_lines[0][value_int]:
                        # Check each combo of alleles and add the count to the matrix

                        if (phased_lines[0][value_int] == "0|0") and (phased_lines[1][value_int] == "0|0"):
                            matrix[0] += 2

                        if (phased_lines[0][value_int] in ["0|1","1|0"]) and (phased_lines[1][value_int] == "0|0"):
                            matrix[0] += 1
                            matrix[1] += 1

                        if (phased_lines[0][value_int] == "0|0") and (phased_lines[1][value_int] in ["0|1", "1|0"]):
                            matrix[0] += 1
                            matrix[2] += 1

                        if (phased_lines[0][value_int] == "0|1" and phased_lines[1][value_int] == "0|1") or (phased_lines[0][value_int] == "1|0" and phased_lines[1][value_int] == "1|0"):
                            matrix[0] += 1
                            matrix[3] += 1

                        if phased_lines[0][value_int] == "1|1" and phased_lines[1][value_int] == "0|0":
                            matrix[1] += 2

                        if (phased_lines[0][value_int] == "0|1" and phased_lines[1][value_int] == "1|0") or (phased_lines[0][value_int] == "1|0" and phased_lines[1][value_int] == "0|1"):
                            matrix[1] += 1
                            matrix[2] += 1

                        if (phased_lines[0][value_int] == "1|1") and (phased_lines[1][value_int] in ["0|1","1|0"]):
                            matrix[1] += 1
                            matrix[3] += 1

                        if (phased_lines[0][value_int] == "0|0") and (phased_lines[1][value_int] == "1|1"):
                            matrix[2] += 1

                        if (phased_lines[0][value_int] in ["0|1","1|0"]) and (phased_lines[1][value_int] == "1|1"):
                            matrix[2] += 1
                            matrix[3] += 1

                        if (phased_lines[0][value_int] == "1|1") and (phased_lines[1][value_int] == "1|1"):
                            matrix[3] += 1

                # Save this matrix and all others in a growing list
                whole_matrix += matrix


            # If there were variants found, save to the output file for locations where data exists
            if sum(whole_matrix) > 0 and min(whole_matrix) >= 0: # old version: sum(whole_matrix) > 0 and min(whole_matrix) > 0
                for mat_ind in range(len(whole_matrix)):
                    whole_matrix[mat_ind] = str(whole_matrix[mat_ind])
                # Count the number of things in either direction
                multipliers.append([str(pair[0]), str(pair[1])]+ID_values+whole_matrix)
            else:
                unused_pairs.append([str(pair[0]), str(pair[1]),"3"])

        # If no variants were found, save to the file to indicate that this pair wasn't used
        else:
            not_in_phasing_data_locations += 1
            if pair_coord_1 in phasing_data_locations:
                unused_pairs.append([str(pair[0]), str(pair[1]), "2"])
            elif pair_coord_2 in phasing_data_locations:
                unused_pairs.append([str(pair[0]), str(pair[1]), "1"])
            else:
                unused_pairs.append([str(pair[0]), str(pair[1]), "0"])


    print("\nOutputting now \n")

    output_file_bad_name = "remaining_"+outname+".txt"
    output_file_good_name = "phasing_matrix_"+outname+".txt"


    # Make the output files
    output_file_good = open(output_file_good_name,'a') # Open a file for things that worked
    output_file_bad = open(output_file_bad_name,'a') # Open a file for things that didn't work

    # Write the file of things that didn't work
    for line in unused_pairs:
        write_line = "\t".join(line) + "\n"
        output_file_bad.write(write_line)
    for line in multipliers:
        write_line = "\t".join(line)+"\n"
        output_file_good.write("ch"+write_line[2:])

loop_end = datetime.now()
print("Chromosome "+str(chrom)+" time: "+str(loop_end-loop_start))

print("\n")

print("Header Line: \n")
output_header_line_text = "SNP1, SNP2, ID1, ID2, ID3, ID4, "
for pop in range(len(populations_list)):
    output_header_line_text += "pop"+str(pop)+".1, "
    output_header_line_text += "pop"+str(pop)+".2, "
    output_header_line_text += "pop"+str(pop)+".3, "
    output_header_line_text += "pop"+str(pop)+".4, "
print(output_header_line_text)

print("\n")

for pop in range(len(populations_list)):
    print("Pop " + str(pop) + " = " + str(populations_list[pop]))

print("\n")

print("Output File Names: ")
print("  "+output_file_good_name)
print("  "+output_file_bad_name)

print("\n\n\n")


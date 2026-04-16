""" head calculate_LD.py -n 10
This script takes the counts of pairs of SNPs and calculates LD statistics for each pair
This can be edited to include more stastics as they become applicable

Inputs:
    -the input file name of counts of pairs of alleles
    -the output file name for the LD statistics

Outputs:
    -A file with each pair of SNPs listed and then the D' and R statistics

"""

# %% Import packages
import sys
import math
from datetime import datetime

start_time = datetime.now()

# %% Define the inputs
input_filename = sys.argv[1]
output_filename = sys.argv[2]

if "*" in input_filename:
    chrs = sys.argv[3]
    start_c = int(chrs[:chrs.index(",")])
    end_c = int(chrs[chrs.index(",")+1:])

    input_filename_list = []
    output_filename_list = []
    in_filename_start = input_filename[:input_filename.index("*")]
    in_filename_end = input_filename[input_filename.index("*")+1:]
    out_filename_start = output_filename[:output_filename.index("*")]
    out_filename_end = output_filename[output_filename.index("*")+1:]

    for ind in range(start_c, end_c+1):
        temp = in_filename_start + str(ind) + in_filename_end
        input_filename_list.append(temp)
        temp = out_filename_start + str(ind) + out_filename_end
        output_filename_list.append(temp)

else:
    input_filename_list = [input_filename]
    output_filename_list = [output_filename]


# %% Define functions
def find_probs(counts_list_input, alleles, pseudocounts = True):
    """
    Parameters
    ----------
    AB : TYPE
        DESCRIPTION.
    Ab : TYPE
        DESCRIPTION.
    aB : TYPE
        DESCRIPTION.
    ab : TYPE
        DESCRIPTION.

    Returns
    -------
    individual_probs : TYPE
        PA, Pa, PB, Pb.
    group_probs : TYPE
        PAB, PAb, PaB, Pab.

    """

    if pseudocounts:
        counts_list = []
        for item in counts_list_input:
            counts_list.append(float(item+1))
    else:
        counts_list = counts_list_input.copy()

    total_counts = sum(counts_list)

    if alleles[0][0] == alleles[1][0]:
        # The listed order is AB, Ab, aB, ab
        io = [0, 1, 2, 3]
    else:
        # The listed order is AB, aB, Ab, ab
        io = [0, 2, 1, 3]

    group_probs = [(counts_list[io[0]])/total_counts,
                   (counts_list[io[1]])/total_counts,
                   (counts_list[io[2]])/total_counts,
                   (counts_list[io[3]])/total_counts]

    individual_counts = [counts_list[io[0]]+counts_list[io[1]],
                         counts_list[io[2]]+counts_list[io[3]],
                         counts_list[io[0]]+counts_list[io[2]],
                         counts_list[io[1]]+counts_list[io[3]]]

    individual_probs = [individual_counts[0]/total_counts,
                        individual_counts[1]/total_counts,
                        individual_counts[2]/total_counts,
                        individual_counts[3]/total_counts]

    return [individual_probs, group_probs]


def LD_scores(counts_list, alleles):

    probs = find_probs(counts_list, alleles)

    # Calculate the D value
    # D = PAB - (PA * PB)
    D = probs[1][0] - (probs[0][0] * probs[0][2])

    if D > 0:
        D_max = min((probs[0][0])*(probs[0][3]),
                    (probs[0][1])*(probs[0][2]))
    else:
        D_max = min((probs[0][0])*(probs[0][2]),
                    (probs[0][1])*(probs[0][3]))

    D_prime = D / D_max

    R = D / math.sqrt(probs[0][0] * probs[0][1] * probs[0][2] * probs[0][3])


    return D_prime, R


# %% Set up for the looping
in_file = open(input_filename_list[0])
first_line = in_file.readline()
first_line = first_line.split("\t")
num_pops = int((len(first_line)-6)/4)

indices = []
start = 6
for i in range(num_pops):
    end = start + 4
    indices.append([start,end])
    start = end



# %% Loop through the chromosomes
for i in range(len(input_filename_list)):

    loop_start = datetime.now()

    in_name = input_filename_list[i]
    out_name = output_filename_list[i]

    in_file = open(in_name,'r')
    out_file = open(out_name,'w')

    # Loop through the input file
    for line in in_file:

        # Read the line
        read_line = line.split("\t")
        alleles = read_line[2:6]

        # Pull out the two variants from the list
        var1 = read_line[0].split(":")
        var2 = read_line[1].split(":")

        # Define the expected double reference and double alternate alleles
        ref_allele = var1[2] + "_" + var2[2]
        alt_allele = var1[3] + "_" + var2[3]

        # Check the expected double reference and double alternate alleles
        if True: #(read_line[2] == ref_allele and read_line[5] == alt_allele) or ("N" in ref_allele or "N" in alt_allele):

            write_line = read_line[:6]

            # Loop through each
            for j in range(num_pops):

                # Define the list of counts
                counts_list = read_line[indices[j][0]:indices[j][1]]
                if "\n" in counts_list[-1]:
                    counts_list[-1] = counts_list[-1][:-1]
                counts_list_int = []
                for item in counts_list:
                    counts_list_int.append(int(item))

                # Calculate the LD scores
                Scores = LD_scores(counts_list_int,alleles)

                write_line.append(str(Scores[0]))
                write_line.append(str(Scores[1]))

            write_line_str = "\t".join(write_line) + "\n"

            out_file.write(write_line_str)


    in_file.close()
    out_file.close()

    loop_end = datetime.now()
    print("Chromosome "+str(i+1)+" time: "+str(loop_end-loop_start))


end_time = datetime.now()
print("Total Time: "+str(end_time-start_time))


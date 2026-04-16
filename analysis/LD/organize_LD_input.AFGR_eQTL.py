""" head make_LD_inuputs.py -n 10
This creates a shell script that will find the LD scores for all of the pairs listed in the input file

Inputs:
    -A previously generated file to append to
    -The new BIRD file with indicators for the correlated QTL SNPs
    -The new QTL file with indicators for the correlated BIRD SNPs

Outputs:
    -LD_inputs.tsv

"""

# %% Import packages
import argparse
import numpy as np
import csv

# %% Define inputs
parser = argparse.ArgumentParser(
                  prog="make_LD_inputs",
                  description="This program finds all pairs of variants to calculate LD for",
                  epilog=""
)
parser.add_argument("-b", "--BIRD", help="The path to the file with BIRD variants")
parser.add_argument("-q", "--QTL", help="The path to the file with QTLs")
parser.add_argument("-o", "--output", help="The path to the output file")
parser.add_argument("-t", "--output_type", default='w', choices=['a', 'w'], help="Write or append an output file")

args = parser.parse_args()

# %% Define the script inputs
#new_bird_filename = sys.argv[2]
#new_qtl_filename = sys.argv[3]
new_bird_filename = args.BIRD
new_qtl_filename = args.QTL

# %% Define funcitons
def find_insert_index(value, list1, restrictions):
    """
    This finds the location in a list where the next value would chronologically fit

    Parameters
    ----------
    value : integer
        The value that needs to be placed into the list.
    list1 : list
        The list of values to be appended.
    restrictions : list
        The start and end coordinates of valid placements of this.

    Returns
    -------
    output_index : int
        The index where the current value should be appended in to.

    """

    # Apply the restrictions
    if len(restrictions) > 0:
        usable_list = []
        for val in list1[restrictions[0]:restrictions[1]]:
            usable_list.append(int(val[1]))
    else:
        usable_list = []
        for val in list1:
            usable_list.append(int(val[1]))

    # Search through the remaining list until you get to a value that is greater than the current searching value
    output_index = 0
    while output_index < len(usable_list) and value > usable_list[output_index]:
        output_index += 1

    if len(restrictions) > 0:
        output_index += restrictions[0]

    return output_index


def read_bird_line(file):
    new_line = file.readline()
    new_ls = new_line.replace("\n", "").split("\t")
    new_data = new_ls[:-1].copy()
    new_indices = new_ls[-1][1:-1].split(", ")
    new_data.append([])
    #print(new_indices)
    if new_indices != [""]:
        for i in new_indices:
            new_data[-1].append(int(i))
    else:
        new_data[-1].append(0)
    #print(new_data)
    return new_data


def print_variant(BIRD, QTL, file):
    pos1 = int(BIRD[1])
    pos2 = int(QTL[1])
    var1 = BIRD[0] + ":" + BIRD[1] + ":" + BIRD[6] + ":" + BIRD[7]
    var2 = QTL[0] + ":" + QTL[1] + ":" + QTL[2] + ":" + QTL[3]
    if pos1 < pos2:
        file.write(var1 + "\t" + var2 + "\n")
    elif pos1 > pos2:
        file.write(var2 + "\t" + var1 + "\n")
    elif var1 < var2:
        file.write(var1 + "\t" + var2 + "\n")
    elif var1 > var2:
        file.write(var2 + "\t" + var1 + "\n")


# %% Load all BIRD variants
BIRD_variants = []
prev_chrom = "NA"
new_line = "A"
BIRD_file = open(new_bird_filename)
while new_line != [[0]]:
    new_line = read_bird_line(BIRD_file)
    if new_line[0] != prev_chrom:
        BIRD_variants.append([])
        prev_chrom = new_line[0]
    BIRD_variants[-1].append(new_line)

BIRD_variants.pop(-1)


# %% Load all QTL variants
new_qtl_reader = csv.reader(open(new_qtl_filename),delimiter="\t")
output_file = open(args.output, args.output_type)
for row in new_qtl_reader:
    chrom_ind = int(row[0][3:]) - 1
    for var in BIRD_variants[chrom_ind]:
        if int(row[-1]) in var[-1]:
            print_variant(var, row, output_file)

"""
# %% Load all QTL variants
current_var_ID = 1
BIRD_file = open(new_bird_filename)
BIRD_lines = [read_bird_line(BIRD_file)]
while current_var_ID in BIRD_lines[-1][-1]:
    BIRD_lines.append(read_bird_line(BIRD_file))

# Load the new QTL file
new_qtl_reader = csv.reader(open(new_qtl_filename),delimiter="\t")
output_file = open(args.output, args.output_type)
for row in new_qtl_reader:
    while len(BIRD_lines) == 0 or current_var_ID in BIRD_lines[-1][-1]:
        BIRD_lines.append(read_bird_line(BIRD_file))
    for var in BIRD_lines:
        if current_var_ID in var[-1]:
            print_variant(var, row, output_file)
    current_var_ID += 1
    while len(BIRD_lines) > 0 and current_var_ID > max(BIRD_lines[0][-1]):
        a = BIRD_lines.pop(0)
"""

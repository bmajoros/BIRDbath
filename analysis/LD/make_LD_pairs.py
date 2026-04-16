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
parser.add_argument("-o", "--output", help="")

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


# %% Load and process the new files

# Load the new QTL file
new_qtl = [[],[],[],[]]
new_qtl_reader = csv.reader(open(new_qtl_filename),delimiter="\t")
for row in new_qtl_reader:
    # Save just the chromosome and position because nothing else will be used in the LD calculations
    new_qtl[0].append(int(row[0][3:]))
    new_qtl[1].append(int(row[1]))
    new_qtl[2].append(row[2])
    new_qtl[3].append(row[3])


old_file = [[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]],[[],[]]]
# Load the new bird file
new_bird = [[],[],[]]
new_bird_reader = csv.reader(open(new_bird_filename),delimiter="\t")
matching_location = False
for row in new_bird_reader:

    # Split the coordinate value into chromosome and position
    coordinate = [row[0], row[1], row[6], row[7]]
    new_bird[0].append(int(row[0][3:]))
    new_bird[1].append(int(row[1]))

    # Find the connections to the QTL coordinates
    connections = row[-1][1:-1].split(", ")

    # Loop through each of the pairs associated with this BIRD position
    for pair_s in connections:

        # Convert the index to a integer
        pair = int(pair_s)

        # Find the QTL coordinate
        QTL_coordinate = [new_qtl[0][pair-1], new_qtl[1][pair-1], new_qtl[2][pair-1], new_qtl[3][pair-1]]
        chrom = QTL_coordinate[0] - 1

        # Check that the chromosomes are the same
        if QTL_coordinate[0] != int(coordinate[0][3:]):

            # Print this pair and the index value to a problems file
            concerns_file = open("Problems_in_input.tsv",'a')

            concerns_file.write("\t".join(row) + "\t" + pair_s + "\n")

            concerns_file.close()


        # If the chromosomes line up
        else:

            # Check the locations of both coordinate values in the old LD pairs file variable

            BIRD_indices = np.where(np.array(old_file[chrom][0]) == int(coordinate[1]))
            QTL_indices = np.where(np.array(old_file[chrom][1]) == int(QTL_coordinate[1]))

            #matching_location = find_same(BIRD_indices[0], QTL_indices[0])

            # If the bird value doesn't exist in the file
            if len(BIRD_indices[0]) == 0:

                # Set the global changes variable to True
                global_changes = True

                # Find the correct bird location
                insert_index = find_insert_index(int(coordinate[1]),old_file[chrom][0],[])

                # put this pair into the output file
                old_file[chrom][0].insert(insert_index, coordinate)
                old_file[chrom][1].insert(insert_index, QTL_coordinate)

            # If the QTL value doesn't exist in the file
            # or
            # If there is no overlap in their indices
            elif len(QTL_indices[0]) == 0 or matching_location == False:

                # Set the global changes variable to True
                global_changes = True

                # Find the correct QTL location
                insert_index = find_insert_index(int(QTL_coordinate[1]),old_file[chrom][1],[BIRD_indices[0][0], BIRD_indices[0][-1]+1])

                # Put this pair into the output file
                old_file[chrom][0].insert(insert_index, coordinate)
                old_file[chrom][1].insert(insert_index, QTL_coordinate)

            # If there is an overlap in their indices
            elif matching_location:

                # This pair is in the input file
                # pass
                pass


# %% Check for changes to add to the new file
print(len(old_file))
print(len(old_file[0]))
print(len(old_file[0][0]))
print(type(old_file[0][0][0]))
print(len(old_file[0][0][0]))
print(old_file[0][0][0])
if global_changes:

    # Overwrite the old generated file to update it
    new_file = open(args.output,'w')

    # Loop through each chromosome
    for chrom_ind in range(len(old_file)):

        chrom = str(chrom_ind + 1)

        # Loop through each pair in the chromosome
        for pair_index in range(len(old_file[chrom_ind][0])):

            # Write this pair into the file
            #new_file.write("chr" + chrom + ":" + str(old_file[chrom_ind][0][pair_index]) + "\t")
            #new_file.write("chr" + chrom + ":" + str(old_file[chrom_ind][1][pair_index]) + "\n")
            new_file.write(":".join(old_file[chrom_ind][0][pair_index]) + "\t")
            #new_file.write(":".join(old_file[chrom_ind][1][pair_index]) + "\n")
            new_file.write("chr" + str(old_file[chrom_ind][1][pair_index][0]) + ":")
            new_file.write(str(old_file[chrom_ind][1][pair_index][1]) + ":")
            new_file.write(old_file[chrom_ind][1][pair_index][2] + ":")
            new_file.write(old_file[chrom_ind][1][pair_index][3] + "\n")

    # Close and save the overwritten file
    new_file.close()


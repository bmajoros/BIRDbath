""" head plotting_groups.py -n 22
This script takes the combinations of BIRD and QTL SNPs and plots them with the specified paramters from the inputs

Inputs:
    -The filename for the BIRD file with indicators for the matching QTL variants
    -The filename for the QTL file with indicators for the matching BIRD variants
    -The threshold for the minimum accepted BIRD p-value
    -The threshold for the maximum accepted QTL p-value
    -The threshold for the minimum accepted D' value
    -The threshold for the minimum accepted R2 value
    -The maximum distance between points (only used as an indicator in the csv file)
    -The nickname of the QTL dataset
        -must be a known dataset because it determines the indices of the beta and p values
        -current options: "AFGR", "altius"
    -The method of phasing to be used:
        -D = multiply by D'
        -R = multiply by R2
        -S = use only the sign
    -Averaging (T/F)
    -Print the resulting graph (T/F)

Outputs:
    -One plot with one point for each
    -One new line in the in the file 'correlations.csv' with the correlation coefficient for this set of parameters

"""

# %% Import packages
from datetime import datetime
import argparse
import csv
import math
import functions
import matplotlib.pyplot as plt
import os
#from scipy.stats import gaussian_kde

start_time = datetime.now()

# %% Define the inputs

parser = argparse.ArgumentParser(
                  prog="write_bedfile",
                  description="This program writes a bedfile from a BIRD file",
                  epilog=""
)
parser.add_argument("-B", "--BIRD", help="The path to the input BIRD file", required=True)
parser.add_argument("-Q", "--QTL", help="The path to the input QTL file", required=True)
parser.add_argument("-b", "--bird_sig", help="The minimum BIRD Preg value", type=float, default=0)
parser.add_argument("-q", "--qtl_sig", help="The maximum QTL p-value", type=float, default=1)
parser.add_argument("-D", "--D_threshold", help="The minimum D value for the LD value", type=float, default=0)
parser.add_argument("-R", "--R_threshold", help="The minimum R value for the LD value", type=float, default=0)
parser.add_argument("-d", "--distance", help="The distance extended outside of the open chromatin region", default="0")
parser.add_argument("-t", "--data_type", help="The label for this dataset", required=True)
parser.add_argument("-L", "--LD_type", help="The method of assessing LD between variants", choices=["S", "D", "R"], default = "S")
parser.add_argument("-A", "--averaging", help="Flag to average across an open chroamtin region", action='store_true')
parser.add_argument("-g", "--print_graphs", help="Flag to print graphs", action='store_true')
parser.add_argument("-c", "--print_csvs", help="Flag to print csv files of the graphs", action='store_true')
parser.add_argument("--LD_path", help="The path to the LD files", required=True)


args = parser.parse_args()

bird_filename = args.BIRD
gtex_filename = args.QTL
bird_threshold = args.bird_sig
qtl_threshold = args.qtl_sig
D_threshold = args.D_threshold
R_threshold = args.R_threshold
dist_val = args.distance
dat_type = args.data_type
phase_type = args.LD_type
if args.averaging:
    averaging = "T"
else:
    averaging = "F"
print_graphs = args.print_graphs
save_csv = args.print_csvs
LD_filename_input = args.LD_path


# Define the percentages of each population
# These are the populations for GTEx, I haven't separated the populations since then
# ORDER : Asian, White, Black
# Blood Populations: 12, 811, 122
# Colon Populations: 8, 73, 482
# Liver Populations: 5, 20, 216
pop_counts = [1]
pop_percentages = []
for item in pop_counts:
    pop_percentages.append(item/sum(pop_counts))


# %% Functions

# Define functions

def find_phase(ld, bird, gtex):
    """
    This searches through the LD phasing file to find the D' and R value for this pair

    Parameters
    ----------
    ld : list
        The input LD files.
    bird : str
        The BIRD SNP in format chr#:#.
    gtex : str
        The GTEx SNP in format chr#:#.

    Returns
    -------
    phase_vals : list
        The D' and R values for this pair.

    """

    # Find the chromosome
    chrom = int(bird[3:bird.index(":")])

    # Find the index within the chromosome
    ID = bird+","+gtex
    ID_ind = ld[chrom-1][0].index(ID)

    # Find the D and R values
    phase_vals = [ld[chrom-1][1][ID_ind], ld[chrom-1][2][ID_ind]]

    return phase_vals

# %% Loading

print("Loading...")

# Load the BIRD data
# The variable has two lists
# The first list is all the indices of the QTL variants with open chromatin regions that overlap this variant
# The second list is all of the data from the BIRD file, where the numerical values are
bird_file = csv.reader(open(bird_filename),delimiter="\t")
bird_inputs = [[],[]]
for line in bird_file:
    if float(line[5]) > bird_threshold:
        bird_inputs[0].append([])
        for item in line[-1][1:-1].split(","):
            bird_inputs[0][-1].append(int(item))
        bird_inputs[1].append([line[0]])
        for item in line[1:-1]:
          try:
            bird_inputs[1][-1].append(float(item))
          except:
            bird_inputs[1][-1].append(item)

# Load the GTEx data and set up the plotting variable (variable described under organizing section)
# This variable has two lists
# The first list is all of the indices of the variant
# The second list is all of the data included that relates to the variant
gtex_file = csv.reader(open(gtex_filename),delimiter="\t")
gtex_inputs = [[],[]]
graphing_variables = [[],[],[],[],[]]
qtl_indices = []
for line in gtex_file:
    if line[0] != "variant":
        if float(line[5]) < qtl_threshold:
            gtex_inputs[0].append(int(line[-1]))
            gtex_inputs[1].append(line[:-1])
            graphing_variables[0].append(line[0] + ":" + line[1] + ":" + line[2] + ":" + line[3]) # QTL variant ID
            graphing_variables[1].append([]) # BIRD variant IDs
            graphing_variables[2].append(float(line[4])) # QTL aFC
            graphing_variables[3].append(0.0) # BIRD aFC
            graphing_variables[4].append(0.0) # aFC difference
            qtl_indices.append(int(line[-1])) # GTEx index


# Load the LD data
# This variable has three separate lists for each chromosome
# The first list has the ID values for the two variants with the given LD values
# The second list has the D' value (weighted average based on populations)
# The third list has the R value (weighted average based on populations)
ld_filename_base = LD_filename_input + "/LD_values_chr"
ld_files = []
for c in range(22):

    chrom = c+1
    # Make a new list to put this LD file into
    ld_files.append([[],[],[]])

    # Load the file
    ld_file = csv.reader(open(ld_filename_base+str(chrom)+".txt"),delimiter="\t")

    for line in ld_file:

        if len(line) > 0:
            # Put the pair ID into the file object
            ld_files[c][0].append(line[0]+","+line[1])
            ld_files[c][1].append(float(0))
            ld_files[c][2].append(float(0))

            # Find the R and D' for each population type
            allR = [float(line[7])]
            allD = [float(line[6])]

            # Put the R and D' value into the file object
            for i in range(len(allR)):
                ld_files[c][1][-1] += allD[i] * pop_percentages[i]
                ld_files[c][2][-1] += allR[i] * pop_percentages[i]


# %% Organizing

print("Organizing...")
# The usual graphing variables variable
# 0. GTEX Location
# 1. BIRD Locations
# 2. GTEx aFC
# 3. BIRD aFC
# 4. Difference in aFC
#graphing_variables = [[],[],[],[],[]]

# Loop through each group as defined by the GTEx SNPs
for b_i in range(len(bird_inputs[0])):
    total_bird = 0

    # Check each correlated QTL variant
    for conn_i in bird_inputs[0][b_i]:
        try:
            q_i = qtl_indices.index(conn_i)
            bird_loc_temp = bird_inputs[1][b_i][0] + ":" + str(int(bird_inputs[1][b_i][1])) + ":" + bird_inputs[1][b_i][6] + ":" + bird_inputs[1][b_i][7]
            gtex_loc_temp = gtex_inputs[1][q_i][0] + ":" + gtex_inputs[1][q_i][1] + ":" + gtex_inputs[1][q_i][2] + ":" + gtex_inputs[1][q_i][3]
            phase = find_phase(ld_files, bird_loc_temp, gtex_loc_temp)

            # Check the phasing threshold
            if abs(phase[0]) > D_threshold and abs(phase[1]) > R_threshold:
                if phase_type == "D":
                    phase_use = phase[0]
                elif phase_type == "R":
                    phase_use = phase[1]
                else:
                    if phase[0] < 0:
                        phase_use = -1
                    elif phase[0] > 0:
                        phase_use = 1
                    else:
                        phase_use = 0
                total_bird += 1
                bird_afc_temp = math.log(bird_inputs[1][b_i][2],2) * phase_use


                graphing_variables[1][q_i].append(bird_loc_temp) # BIRD name
                graphing_variables[3][q_i] += bird_afc_temp # BIRD aFC

        except:
            pass

delete_inds = []
for q_i in range(len(graphing_variables[0])):
    if len(graphing_variables[1][q_i]) > 0:
        if args.averaging:
            bird_afc = graphing_variables[3][q_i]
            graphing_variables[3][q_i] = bird_afc/float(len(graphing_variables[1][q_i]))
        graphing_variables[4][q_i] = graphing_variables[2][q_i] - graphing_variables[3][q_i]

    else:
        delete_inds.append(q_i)

delete_inds.reverse()
for d_i in delete_inds:
    graphing_variables[0].pop(d_i)
    graphing_variables[1].pop(d_i)
    graphing_variables[2].pop(d_i)
    graphing_variables[3].pop(d_i)
    graphing_variables[4].pop(d_i)

"""
for ind in range(len(gtex_inputs[0])):

    i = gtex_inputs[0][ind]

    # Find the relevant lines
    gtex_relevant = gtex_inputs[1][ind]

    # Save the GTEx SNP
    gtex_loc_temp = gtex_relevant[0]+":"+gtex_relevant[1]

    # Save the GTEx aFC
    gtex_afc_temp = float(gtex_relevant[2])
    gtex_pval = float(gtex_relevant[3])

    if gtex_pval < qtl_threshold:

        total_bird = 0
        bird_afc_temp = 0
        bird_locs_temp = []
        # Save the bird information
        for b in range(len(bird_inputs[0])):

            try:
                # Check that this bird value is in the current group
                if i in bird_inputs[0][b]:

                    # Save the BIRD SNPs
                    bird_loc_temp = bird_inputs[1][b][0].replace("@",":")
                    bird_locs_temp.append(bird_loc_temp)

                    # Check the bird p threshold
                    if bird_inputs[1][b][4] >= bird_threshold:

                        # Save the BIRD
                        phase = find_phase(ld_files, bird_loc_temp, gtex_loc_temp)

                        # Check the phasing threshold
                        if abs(phase[0]) > D_threshold and abs(phase[1]) > R_threshold:
                            if phase_type == "D":
                                phase_use = phase[0]
                            elif phase_type == "R":
                                phase_use = phase[1]
                            else:
                                if phase[0] < 0:
                                    phase_use = -1
                                elif phase[0] > 0:
                                    phase_use = 1
                                else:
                                    phase_use = 0
                            total_bird += 1
                            bird_afc_temp += math.log(bird_inputs[1][b][1],2) * phase_use


            except:
                pass

        # Average if asked to
        if total_bird > 0 and averaging == "T":
            graphing_variables[0].append(gtex_loc_temp)
            graphing_variables[1].append(bird_locs_temp)
            graphing_variables[2].append(gtex_afc_temp)
            graphing_variables[3].append(bird_afc_temp/total_bird)
            graphing_variables[4].append(graphing_variables[2][-1]-graphing_variables[3][-1])

        elif total_bird > 0:
            graphing_variables[0].append(gtex_loc_temp)
            graphing_variables[1].append(bird_locs_temp)
            graphing_variables[2].append(gtex_afc_temp)
            graphing_variables[3].append(bird_afc_temp)
            graphing_variables[4].append(graphing_variables[2][-1]-graphing_variables[3][-1])
"""

# %% Statistics

# Check if the output file has been started yet
files_in_directory = os.listdir()
if "correlations.csv" not in files_in_directory:
    # Write the header if it hasn't been started
    header_file = open("correlations.csv", 'w')
    header_file.write("QTL_name,distance,bird_p,qtl_p,D_min,R_min,phase_type,avg,B+Q+,B+Q-,B-Q+,B-Q-,score,slope,intercept\n")
    header_file.close()

# Open an output file to append to
output_file = open("correlations.csv",'a')



# Start the output line with the parameters
output_line = "" # initialize the output line
output_line += str(dat_type) + "," # put the nickname of the open chromatin dataset into the otuput line
output_line += str(dist_val) + "," # put the max distance value into the output line
output_line += str(bird_threshold) + "," # put the bird threshold into the output line
output_line += str(qtl_threshold) + "," # put the qtl theshold into he output line
output_line += str(D_threshold) + "," # put the D prime threshold into the output line
output_line += str(R_threshold) + "," # put the R squared threhsold into the output line
output_line += str(phase_type) + "," # put the phasing type into the output line
output_line += averaging + "," # put an indication of the averaging across the OC into the output line

# Assess the sign accuracy

# Make a table for the chi squared test
chi_2_table = [[0,0],[0,0],[0,0]] # initialize a table variable, the first row is positive bird and the first column is positive GTEx

# Place everything that would work into the table

# Loop through all of the haplotypes
for index_num in range(len(graphing_variables[0])):

    # Check the bird sign
    if graphing_variables[3][index_num] > 0:
        index_1 = 0
    elif graphing_variables[3][index_num] < 0:
        index_1 = 1
    else:
        index_1 = 2

    # Check the GTEx sign
    if graphing_variables[2][index_num] > 0:
        index_2 = 0
    else:
        index_2 = 1

    # Change the count in the table
    chi_2_table[index_1][index_2] += 1

# Calculate the chi squared value

chi2_array = [chi_2_table[0][0], chi_2_table[0][1], chi_2_table[1][0],
              chi_2_table[1][1], chi_2_table[2][0], chi_2_table[2][1]]

# Put the chi squared table into the output line
output_line += str(chi2_array[0]) + ","
output_line += str(chi2_array[1]) + ","
output_line += str(chi2_array[2]) + ","
output_line += str(chi2_array[3]) + ","


# Create a regression model to get a little closer to the actual value
try: # this may not work because there may be too few values, which is the reason there is a try/except

    slope, intercept = functions.Linear_Fitting(graphing_variables[2], graphing_variables[3]) # Fit a linear regression model and get the slope and intercept
    predicted_values = [] # Make a list for the predicted values
    for input_value in graphing_variables[2]: # For each input X value for the linear regression
        predicted_values.append( (input_value * slope) + intercept) # Predict the Y value from the linear regression

    # Assess the accuracy of the regression model
    score = functions.R_squared_calculation(graphing_variables[3], predicted_values)
except:
    score = "Can't compute"
    slope = "NA"
    intercept = "NA"


output_line += str(score) + ","
output_line += str(slope) + ","
output_line += str(intercept)
# Finish the summary line and write it to the file
output_line += "\n"
output_file.write(output_line)

# Close the output file
output_file.close()


# %% Plotting

#if print_graphs and score != "Can't compute":
if score != "Can't compute":
    fig_name = dist_val+"_"+str(bird_threshold)+"_"+str(qtl_threshold)+"_"
    fig_name+= str(D_threshold)+"_"+str(R_threshold)+"_"+phase_type+"_"
    fig_name+= "avg"+averaging

    # Plot beta v theta
    x_temp = graphing_variables[2]
    y_temp = graphing_variables[3]
    #xy = np.vstack([x_temp, y_temp])
    #z = gaussian_kde(xy)(xy)
    #plt.scatter(x_temp,y_temp, c=z, s=1)
    #plt.colorbar()
    plt.scatter(x_temp,y_temp, s=3)
    plt.plot(x_temp, predicted_values, c='black')
    plt.plot(x_temp, [0]*len(x_temp), c='black', alpha=0.5)
    plt.plot([0]*len(y_temp), y_temp, c='black', alpha=0.5)
    plt.xlabel('caQTL effect size', fontsize=18)
    plt.ylabel('BIRD effect size', fontsize=18)
    coord1 = 0.05
    if slope > 0:
        coord2 = 0.95
    else:
        coord2 = 0.05
    #plt.text(coord1, coord2, str(score))
    #plt.title(titles[plot_num]+" effect sizes")
    #mini_dataframe.plot(kind='scatter',x=x_temp,y=y_temp,title=titles[plot_num]+" p values",xlabel="bird effect size",ylabel="eQTL effect size",color=z)
    plt.savefig("graphs/"+fig_name+"_efsz.png")
    plt.clf()

if save_csv and score != "Can't compute":
    output_filename = "graphs_csvs/" + dist_val+"_"+str(bird_threshold)+"_"+str(qtl_threshold)+"_"
    output_filename+= str(D_threshold)+"_"+str(R_threshold)+"_"+phase_type+"_"
    output_filename+= "avg"+averaging + ".csv"
    output_file = open(output_filename, 'w')

    write_line = "QTL_SNP,BIRD_SNP,QTL_afc,BIRD_afc,difference,total_BIRD\n"
    output_file.write(write_line)

    for i in range(len(graphing_variables[0])):
        write_line = str(graphing_variables[0][i]) + ","
        write_line+= str(graphing_variables[1][i]).replace(",",";") + ","
        write_line+= str(graphing_variables[2][i]) + ","
        write_line+= str(graphing_variables[3][i]) + ","
        write_line+= str(graphing_variables[4][i]) + "\n"

        output_file.write(write_line)

    output_file.close()

# %% Finishing

print("Done!")
end_time = datetime.now()
print("Time: "+str(end_time-start_time))
print("\n\n")


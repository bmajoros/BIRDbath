# Import packages
import argparse 

# Define inputs
parser = argparse.ArgumentParser(
                  prog="write_bedfile",
                  description="This program writes a bedfile from a BIRD file",
                  epilog=""
)
parser.add_argument("-r", "--remove", help="Flag to indicate that the colons should be removed", action="store_true")
parser.add_argument("-i", "--input", help="The path to the input file", required=True)
parser.add_argument("-o", "--output", help="The path to the output file", required=True)

args = parser.parse_args()

if args.remove:
    def process(in_line):
        outline = in_line.replace(":", "\t")
        return outline

else:
    def process(in_line):
        ls = in_line.split("\t")
        outline = ":".join(ls[:4]) + "\t" + ":".join(ls[4:])
        return outline 

in_file = open(args.input)
out_file = open(args.output, 'w')
for line in in_file:
    out_file.write(process(line))

out_file.close()

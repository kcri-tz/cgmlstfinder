#!/usr/bin/env python3
import sys, os, shutil, subprocess, argparse
import numpy as np
from ete3 import Tree
from ete3.parser.newick import NewickError

def get_distance(mlst1, mlst2):

    mlst1 = np.array(mlst1)
    mlst2 = np.array(mlst2)

    # Mask the uncalled alleles
    mask = (mlst1 !="-") & (mlst2 !="-")

    mlst1 = mlst1[mask]
    mlst2 = mlst2[mask]

    call_length = len(mlst1)

    dist_1_2 = mlst1==mlst2
    dist = call_length - dist_1_2.sum()
    return dist

def distance_matrix(infile, outfile_name):
    allele_profiles = []
    dis_matrix = []
    node_IDs = {}
    real_names = {}
    sample_no = 0
    # Open allele matrix file
    with open(infile) as allel_file:
        header = allel_file.readline()    
        samples_profile = allel_file.readlines()

    loci = header.strip().split("\t")[1:]
    first_sample = True
    first_profile = np.array([])
    for profile in samples_profile:
        profile = profile.strip().split("\t")
        sample = profile[0]
        alleles = profile[1:]

        # Sampel ID for distance matrix must be 10 char, save real sample name in names dict
        sample_ID = "0000000000" + str(sample_no)
        sample_ID = sample_ID[-10:]
        node_IDs[sample_no] = sample_ID
        real_names[sample_ID] = sample
        sample_no += 1
        # Calculate ditance to first sample
        if first_sample == True:
            first_alleles = alleles
            dis = "0"
            first_sample = False
        else:
            dis = str(get_distance(first_alleles, alleles))
        allele_profiles.append(alleles)

        dis_matrix.append([dis])
    # Open outfile to write
    outfile = open(outfile_name, "w")
    outfile.write(str(len(node_IDs)) + "\n")
    outfile.write(node_IDs[0] + "\n")

    # fill in distances in the lower triangle and write the outfile
    for i in range(1, len(dis_matrix)):
        for j in range(1, len(dis_matrix)):
            if j < i:
                i_j_dis = get_distance(allele_profiles[j], allele_profiles[i])
                dis_matrix[i].append(str(i_j_dis))
        outfile.write(node_IDs[i] + "\t" + "\t".join(dis_matrix[i]) + "\n")

    outfile.close()

    return real_names

def create_tree(matrix_path, nj_program):

    # Create "tmp" directory
    try:
       tmp_dir = os.path.join(os.getcwd(), "tmp")
       os.makedirs(tmp_dir)
    except OSError as error: 
        if "outfile" in os.listdir("tmp") or "outtree" in os.listdir("tmp"):
            sys.exit("The directory, {}, already contains the files 'outfile' and 'outtree'. In order to make a new tree please remove these files.".format(tmp_dir) )
    os.chdir("tmp")

    # Use of neighbor program from python:
    matpath = os.path.abspath(matrix_path)
    inputstr = "{0}\nL\nY\n".format(matpath).encode('utf-8')
    proc = subprocess.Popen(nj_program, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    proc.stdin.write(inputstr)

    # Wait for it. it gets stuck in simple wait() if fails
    (stdoutdata, stderrdata) = proc.communicate()
    if proc.returncode:
        sys.exit("Neighbor program failed.")

    os.chdir(outdir)

def decode_dist_matrix(seqid2name, dist_mat_file):
    mat_cont = []
    with open(dist_mat_file, "r") as mat_f:
        for line in mat_f:
            first_word = line.split()[0]
            try:
                mat_cont.append(line.replace(first_word, seqid2name[first_word]))
            except KeyError:
                mat_cont.append(line)

    with open(dist_mat_file, "w") as op:
        print("".join(mat_cont), file=op)

    return

########
# MAIN #
########

parser = argparse.ArgumentParser(description="")

parser.add_argument("-i", "--infile", help="Allele matrix output for cgMLSTFinder", default="results.txt")
parser.add_argument("-o", "--outdir", help="Output directory", default = ".")
parser.add_argument("-n", "--NJ_path", help="Path to executable Neighbor Joining program") #, default = "/home/projects/cge/apps/NDtree/neighbor")

args = parser.parse_args()

# Define input and output files and neighbor joining program path
infile = os.path.abspath(args.infile)
outdir = os.path.abspath(args.outdir)
dis_matrix_file = outdir + "/dis_matrix.txt"
newick_file = outdir + "/allele_tree.newick"

# Check that nj program is executable
nj_program = args.NJ_path
if not shutil.which(nj_program):
   sys.exit("Path to neighbor joining program '{}' must be extecuatbe".format(nj_program))

os.chdir(outdir)

# Generate distance matix
real_names = distance_matrix(infile, dis_matrix_file)

create_tree(dis_matrix_file, nj_program)

fulltree_file = os.path.join(outdir, "tmp/outtree")

try:
    tree = Tree(fulltree_file)
except NewickError as e:
    sys.exit("Couldn't open {0}".format(fulltree_file))

for node in tree.traverse():
    if node.name:
        node.name = real_names[node.name]

# Clean up, remove tmp dir
rm_cmd = "rm -r tmp/outtree tmp/outfile"
os.system(rm_cmd)
if os.listdir("tmp") == []: 
    os.rmdir("tmp")

with open(newick_file, "w") as nf:
    print(tree.write(), file=nf)

decode_dist_matrix(real_names, dis_matrix_file)

#print(tree)

print("Done")

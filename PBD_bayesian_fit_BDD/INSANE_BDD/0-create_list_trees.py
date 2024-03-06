import os
import glob


outdir = "C:/Users/pveron/Output_clusters/PBD_analog/12152"


# Save wd and change it
wd = os.getcwd()
os.chdir(outdir)

# List of file names matching pattenr
list_files = glob.glob("stree_random*.nwk")

# Change again wd and save file 
os.chdir(wd)
with open("list_trees.txt", "w") as f:
    f.writelines(line + '\n' for line in list_files)


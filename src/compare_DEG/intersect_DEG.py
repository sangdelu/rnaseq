import numpy as np
import pandas as pd
import sys

from DEtable import DEtable

# build instances
df_all = DEtable(sys.argv[1])
df_pc = DEtable(sys.argv[2])

# open output
out_handle = open(sys.argv[3],"w")

# Extract information
df_all.load_table()
df_all.extract_DEG()
df_pc.load_table()
df_pc.extract_DEG()


# compare DEG
intersect_up = set(df_all.up).intersection(set(df_pc.up))
intersect_down = set(df_all.down).intersection(set(df_pc.down))
intersect_info = f"{len(intersect_up)},{len(intersect_down)},{len(intersect_up) + len(intersect_down)}"

all_info = f"{len(df_all.up)},{len(df_all.down)},{len(df_all.up) + len(df_all.down)}"
pc_info = f"{len(df_pc.up)},{len(df_pc.down)},{len(df_pc.up) + len(df_pc.down)}"

# Write statistics to output
out_handle.write("intersect_up,intersect_down,intersect_all;all_up,all_down,all_DE;protein_coding_up,protein_coding_down,protein_coding_all\n")
out_handle.write(";".join([intersect_info,all_info,pc_info])+"\n")

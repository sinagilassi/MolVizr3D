# import packages/modules
import os
import json
import numpy as np
import molvizr3d as mv3d
import pprint

# check version
# print(mv3d.__version__)

# sdf file
# Methanol
# Conformer3D_COMPOUND_CID_887
# Benzene
# Conformer3D_COMPOUND_CID_241
# Glycol Dimethacrylate
# Conformer3D_COMPOUND_CID_7355
# 2-(2-Methoxyethoxy)ethanol
# Conformer3D_COMPOUND_CID_8134
# Bromobenzyl cyanide
# Conformer3D_COMPOUND_CID_22044
# Amineptine
# Conformer3D_COMPOUND_CID_34870
# Amidephrine
# Conformer3D_COMPOUND_CID_15010
# N,N-Dimethylformamide
# Conformer3D_COMPOUND_CID_6228
# Acetanilide
# Conformer3D_COMPOUND_CID_904
sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_904.sdf'
sdf_file = os.path.join(os.getcwd(), sdf_file_name_1)

# visualize compound by sdf file
# mv3d.td(sdf_file)

# visualize compound by inchi
# mv3d.td_by_inchi(
#     'InChI=1S/C14H22O6/c1-11(2)13(15)19-9-7-17-5-6-18-8-10-20-14(16)12(3)4/h1,3,5-10H2,2,4H3', display_legend=False)

# network
# mv3d.check_functional_group(
#     'test\Conformer3D_COMPOUND_CID_7979.sdf', functional_groups=['ether'])


res = mv3d.check_functional_group(sdf_file, res_format='dataframe')
print(res)

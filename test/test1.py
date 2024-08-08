# import packages/modules
import os
import json
import numpy as np
import molvizr3d as mv3d
import pprint

# check version
# print(mv3d.__version__)

# sdf file
sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_7979.sdf'
sdf_file = os.path.join(os.getcwd(), sdf_file_name_1)

# visualize compound by sdf file
# mv3d.td(sdf_file)

# visualize compound by inchi
mv3d.td_by_inchi(
    'InChI=1S/C14H22O6/c1-11(2)13(15)19-9-7-17-5-6-18-8-10-20-14(16)12(3)4/h1,3,5-10H2,2,4H3', display_legend=False)

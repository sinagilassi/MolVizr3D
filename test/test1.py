# import packages/modules
import os
import json
import numpy as np
import molvizr3d as mv3d
import pprint

# check version
# print(mv3d.__version__)

# sdf file
# sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_22044.sdf'
sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_7979.sdf'
# sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_887.sdf'
# sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_8134.sdf'
sdf_file = os.path.join(os.getcwd(), sdf_file_name_1)
# print(sdf_file)

# visualize
mol = mv3d.td(sdf_file)
# print(mol['atom_elements'])
# print(mol['xyz_list'])
# print(mol['bond_numbers'])
# atom block
# atom_block = mol['atom_block']
# bond block
# bond_block = mol['bond_block']
# xyz list
# atom_xyz = mol.atom_xyz
# save
# Save to a text file (comma-separated by default)
# np.savetxt(os.path.join(os.getcwd() + '/test/' + 'xyz_list.txt'), atom_xyz)
# np.save(os.path.join(os.getcwd() + '/test/' + 'xyz_list.txt'), atom_xyz)

# save
# Save as JSON
# with open('atom_block.json', 'w') as f:
#     json.dump(atom_block, f)

# with open('bond_block.json', 'w') as f:
#     json.dump(bond_block, f)
# print(type(mol))
plot_summary = mol.view3d()
# pprint.pprint(plot_summary)

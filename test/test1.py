# import packages/modules
import os
import json
import molvizr3d as mv3d

# check version
# print(mv3d.__version__)

# sdf file
sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_241.sdf'
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
# save
# Save as JSON
# with open('atom_block.json', 'w') as f:
#     json.dump(atom_block, f)

# with open('bond_block.json', 'w') as f:
#     json.dump(bond_block, f)
print(type(mol))
mol.view3d()

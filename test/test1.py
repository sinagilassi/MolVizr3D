# import packages/modules
import os
import molvizr3d as mv3d

# check version
# print(mv3d.__version__)

# sdf file
sdf_file_name_1 = 'test\Conformer3D_COMPOUND_CID_241.sdf'
sdf_file = os.path.join(os.getcwd(), sdf_file_name_1)
# print(sdf_file)

# visualize
mol = mv3d.td(sdf_file)
print(mol['atom_elements'])
print(mol['xyz_list'])
print(mol['bond_numbers'])

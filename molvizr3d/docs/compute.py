
# COMPUTATIONAL
# ---------------

# import libs
import numpy as np


class Compute():
    '''
    computational chemistry
    '''

    def __init__(self):
        pass

    @staticmethod
    def cal_atoms_distance(xyzAtom1, xyzAtom2):
        '''
        calculate distance between two atoms (center to center)
        '''
        return np.linalg.norm(xyzAtom1-xyzAtom2)

    @staticmethod
    def atoms_distance_matrix(xyzList, atomName):
        '''
        build a matrix containing a matrix of distance between two different atoms

        args:
            xyzList: xyz list of atoms
            atomName: atom name list such as ['C','H','H','H','H']
        '''
        # atom no
        atomNo = len(atomName)
        atomLength = np.zeros((atomNo, atomNo))
        for i in range(atomNo):
            # atom1
            atom1XYZ = xyzList[i]
            for j in range(atomNo):
                # atom2
                atom2XYZ = xyzList[j]
                if j != i:
                    _length = Compute.cal_atoms_distance(atom1XYZ, atom2XYZ)
                else:
                    _length = 0
                # save
                atomLength[i, j] = _length

        # res
        return atomLength

    @staticmethod
    def atoms_distance(xyzList, atomName, atom_symbols, atom_index):
        '''
        build a matrix containing a matrix of distance between two different atoms

        args:
            xyzList: xyz list of atoms
            atomName: atom name list such as ['C','H','H','H','H']
        '''
        try:
            # size
            atom_symbols_size = len(atom_symbols)
            atom_index_size = len(atom_index)

            # check
            if atom_symbols_size > 2 or atom_index_size > 2:
                raise Exception(
                    'atom symbol/index list must have two elements.')

            if atom_index_size == 0:
                # find the index of the first element
                atom1Index = np.where(atomName == atom_symbols[0])[0][0]
                atom2Index = np.where(atomName == atom_symbols[1])[0][0]
                # set
                atom1XYZ = xyzList[atom1Index]
                atom2XYZ = xyzList[atom2Index]
                _length = Compute.cal_atoms_distance(atom1XYZ, atom2XYZ)
            elif atom_index_size == 2:
                # set
                atom1XYZ = xyzList[atom_index[0]]
                atom2XYZ = xyzList[atom_index[1]]
                _length = Compute.cal_atoms_distance(atom1XYZ, atom2XYZ)
            elif atom_index_size > 2:
                raise Exception(
                    'atom symbol/index list must have two elements.')

            # res
            return _length
        except Exception as e:
            raise


# UTILITY FUNCTION
# ------------------


def CalculateMolecularMass(atom_elements, element_source):
    '''
    calculate molecular mass [g/mol]

    args:
        atom_elements: such as C, H
        element_source: periodic element table
    '''
    try:
        # # check
        # if element_source is None:
        #     # load periodic table of elements
        #     element_source = elements()
        # res
        res = 0
        # atom element size
        atomElementsSize = len(atom_elements)
        # atom property
        _atom_property = ['AtomicMass']
        # loop
        for i in range(atomElementsSize):
            _atom_symbol = atom_elements[i]
            _atom_prop_res = element_source.atom_properties(
                _atom_symbol, _atom_property)
            _atom_prop_val = _atom_prop_res.get(_atom_property[0])
            res += float(_atom_prop_val)

        return res
    except Exception as e:
        raise

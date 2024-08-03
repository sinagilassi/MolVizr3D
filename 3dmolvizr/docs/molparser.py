# MATERIAL PARSER
# ----------------

# import libs
import re
import numpy as np
import copy
# internals
from .config import OBS_POSITIONS
from .structure import Structure
from .element import Element
from .utility import Utility


class MolParser():
    '''
    parse different formats of molecule files such as sdf, json, ...
    '''
    _mat = {}

    def __init__(self, filepath):
        self.filepath = filepath

    @property
    def mat(self):
        return self._mat

    @mat.setter
    def mat(self, value):
        self._mat = copy.deepcopy(value)

    def read_file(self, sourceContent={}):
        '''
        read structure file 
            1. from file
            2. download from api

        args:
            filepath: full file name with directory 
        '''
        try:
            # get file path
            filePath = self.filepath

            # check
            if filePath:
                # open file and get content
                fileContent, fileDir, fileName, fileFormat = Utility.OpenFile(
                    filePath)
            else:
                # read string/json and ... files
                fileContent, fileFormat = Utility.ReadContent(
                    sourceContent)

            # method selection
            parserFun = {
                'sdf': self.sdf_parser,
                'json': self.json_parser
            }

            # parse file
            parserSelection = parserFun.get(fileFormat)
            parserRes = parserSelection(fileContent)
            return parserRes

        except Exception as e:
            raise Exception(e)

    def sdf_parser(self, sdfSource, sdfVersion='V2000'):
        '''
        parse sdf file

        args:
            sdfSource: file read return
            sdfVersion: the compatible version 
        '''
        # decode binary
        # if binary file
        # sdfSource = sdfSource.decode('utf-8')
        # if text file
        sdfSource = sdfSource
        # create list
        sdfSourceList = sdfSource.splitlines()

        # header block
        headerBlock = sdfSourceList[0:3]

        # counts line
        countsLine = sdfSourceList[3]
        if countsLine.find('V2000') == -1:
            raise Exception(
                'SDF file version is not compatible with this method, import 2000 version.')
        # check
        countsLineList = countsLine.split()

        # find 'M END'
        MENDi = -1
        MENDi = sdfSourceList.index('M  END')
        # connection table
        if MENDi != -1:
            connectionTable = sdfSourceList[4:MENDi]

        # *** check
        _countsLineListSize = len(countsLineList)
        if _countsLineListSize == 10:
            # atom no
            atomNo = int(countsLineList[0])
            # bond no
            bondNo = int(countsLineList[1])
        elif _countsLineListSize == 9:
            _connectionTableLines = connectionTable
            _connectionTableLinesSize = len(_connectionTableLines)
            # first record
            _firstRecord = _connectionTableLines[0]
            _firstRecordSize = len(_firstRecord)
            _secondRecord = _connectionTableLines[0]
            _secondRecordSize = len(_secondRecord)
            if _firstRecordSize != _secondRecordSize:
                raise Exception('sdf file is not coded correctly.')
            # loop
            for l in range(_connectionTableLinesSize):
                _loopRecordSize = len(_connectionTableLines[l])
                if _loopRecordSize < _firstRecordSize:
                    # atom no
                    atomNo = int(l)
                    # bond no
                    bondNo = int(str(countsLineList[0]).split(str(atomNo))[1])
                    # break
                    break
        else:
            raise Exception('3rd line of the sdf file is not coded correctly.')

        # element rows
        elementRows = connectionTable[0:atomNo]
        # bond rows
        bondRows = connectionTable[atomNo:]

        # other vars
        compoundPropertiesList = MolParser.__var_finder(sdfSource)
        # dict vars
        compoundProperties = MolParser.__var_analyzer(compoundPropertiesList)
        # extract:
        PUBCHEM_COMPOUND_CID = compoundProperties.get('PUBCHEM_COMPOUND_CID')
        PUBCHEM_IUPAC_NAME = compoundProperties.get('PUBCHEM_IUPAC_NAME')
        PUBCHEM_EXACT_MASS = compoundProperties.get('PUBCHEM_EXACT_MASS')
        PUBCHEM_MOLECULAR_FORMULA = compoundProperties.get(
            'PUBCHEM_MOLECULAR_FORMULA')
        PUBCHEM_MOLECULAR_WEIGHT = compoundProperties.get(
            'PUBCHEM_MOLECULAR_WEIGHT')

        atoms = []
        atomList = []
        xyzList = []

        # atoms position
        for i in range(atomNo):
            _atomRow = elementRows[i].split()
            # position
            _x = float(_atomRow[0])
            _y = float(_atomRow[1])
            _z = float(_atomRow[2])
            # name
            _name = _atomRow[3]

            # atom list
            atomList.append(_name)

            # atom info
            atom = {
                'id': i,
                'symbol': _name,
                'index': i,
                'x': _x,
                'y': _y,
                'z': _z,
                'position': {
                    'x': _x,
                    'y': _y,
                    'z': _z
                },
                'xyz': [_x, _y, _z]
            }
            # save atom info
            atoms.append(atom)
            # xyz
            xyzList.append([_x, _y, _z])

        # set
        xyzList = np.array(xyzList)

        # object base
        objectBaseCoordinate = Structure.CenterPoints(xyzList)
        # print(f"objectBaseCoordinate: {objectBaseCoordinate}")

        # move to the center [0,0,0]
        xyzCenterList, movingCoordinate = Structure.CenterObject(
            xyzList, objectBaseCoordinate)
        # print(f"xyzCenterList: {xyzCenterList} \n movingCoordinate: \n {movingCoordinate}")

        # create mat formula
        matFormula = Structure.create_formula(atomList)
        # calculate molecular mass
        # matMass = 1

        # bond analysis
        bondList = []
        atomBonds = []

        # atomNo: *** atom id in the structure ***
        for i in range(atomNo):
            # name
            _nameAtom1 = str(atomList[i])
            # index [not duplicated element]
            # _nameAtom1Index = atomList.index(_nameAtom1)
            _nameAtom1Index = i+1
            # check bond type
            for j in range(bondNo):
                _bondRow = bondRows[j].split()
                # starts from 1 to ... (not 0)
                _nameAtomBondRow = int(_bondRow[0])
                # check bond exist
                if _nameAtom1Index == _nameAtomBondRow:
                    # atom 2 index
                    _indexAtom = int(_bondRow[1])
                    # atom 2 name
                    _nameAtom2 = str(atomList[_indexAtom-1])
                    # bound type
                    _bondType = int(_bondRow[2])
                    # str bond
                    _bondName = _nameAtom1 + _nameAtom2
                    # atom bond
                    atomBonds.append(
                        (_indexAtom, _nameAtom2, _bondName, _bondType))

            if len(atomBonds) > 0:
                # save
                _bondList = {
                    'id': _nameAtom1Index,
                    'symbol': _nameAtom1,
                    'bonds': atomBonds
                }
                bondList.append(_bondList)

            # reset
            atomBonds = []
            _bondList = {}

        # res
        res = {
            'header_block': headerBlock,
            'counts_line': countsLine,
            'atom_numbers': atomNo,
            'mat_cid': PUBCHEM_COMPOUND_CID,
            'mat_name': PUBCHEM_IUPAC_NAME if PUBCHEM_IUPAC_NAME is not None else '',
            'mat_formula': PUBCHEM_MOLECULAR_FORMULA if PUBCHEM_MOLECULAR_FORMULA is not None else matFormula,
            'mat_mass': PUBCHEM_EXACT_MASS if PUBCHEM_EXACT_MASS is not None else PUBCHEM_MOLECULAR_WEIGHT,
            'atom_names': atomList,
            'atom_elements': atomList,
            'bond_numbers': bondNo,
            'atom_block': atoms,
            'bond_block': bondList,
            'xyz_list': xyzList,
            'xyz_center_list': xyzCenterList,
            'compound_properties': compoundProperties
        }

        # return
        return res

    def json_parser(self, jsonSource, id_sort=True):
        '''
        parse json file

        args:
            jsonSource: dict file 
        '''
        try:
            dictSource = jsonSource['PC_Compounds'][0]

            # analyze json file
            # *** id node
            if 'id' in dictSource:
                _idNode = dictSource['id']
                cid = self.__json_parser_id(_idNode)

            # *** atoms
            if 'atoms' in dictSource:
                _atomsNode = dictSource['atoms']
                # atom ids, element atomic number
                atomIds, _elementAtomicNumber = self.__json_parser_atoms(
                    _atomsNode)

            # atom numbers
            atomNo = len(atomIds)

            # atom list (element list)
            elementListRes, elementList = self.__find_element_symbol(
                _elementAtomicNumber)

            # *** bonds
            if 'bonds' in dictSource:
                _bonds = dictSource['bonds']
                # atom1 id, atom2 id, bond types
                atomId1, atomId2, bondTypes, bondMatrix = self.__json_parser_bonds(
                    _bonds)

            # bond numbers
            bondNo = len(atomId1)

            # *** coords
            if 'coords' in dictSource:
                _coords = dictSource['coords']  # list
                # coords type, aid, xyzList
                coords_type, coords_aid, xyzList, structure_type = self.__json_parser_coords(
                    _coords)

            # *** charge
            if 'charge' in dictSource:
                _charge = dictSource['charge']  # scaler
                charge = self.__json_parser_charge(_charge)

            # *** props
            if 'props' in dictSource:
                _props = dictSource['props']  # list of dict
                propDict = self.__json_parser_props(_props)

            # *** count
            if 'count' in dictSource:
                _count = dictSource['count']  # dict
                countList = self.__json_parser_count(_count)

            # interpret
            __json_atom_position_res = self.__json_atom_position(
                atomNo, elementList, xyzList)

            # bond block
            __json_atom_bondblock_res = self.__json_atom_bondblock(
                atomNo, bondNo, elementList, bondMatrix)

            # bondBlock: used for Matview class
            # atomDetails, bondBlock, xyzCenterList = self.__json_atom_analyzer(
            #     atomNo, bondNo, elementList, xyzList, bondMatrix)

            # *** origin info
            origin_info = {
                'atom_elements': __json_atom_position_res.get('elementList'),
                'atom_atomic_number': _elementAtomicNumber,
                'atom_details': __json_atom_position_res.get('atomDetails'),
                'bond_block': __json_atom_bondblock_res.get('bondBlock'),
                'bond_list': __json_atom_bondblock_res.get('bondMatrix'),
                'xyz_list': __json_atom_position_res.get('xyzList'),
                'xyz_center_list': __json_atom_position_res.get('xyzCenterList'),
            }

            # *** define new ids for mat and update
            # *** xyzList, bondMatrix, elementlist, elementAtomicNumber
            # *** bond_list = bondMatrix
            # *** atom_block = atom_details
            # if id_sort:
            mat_position_info, atom_id_conversion, xyz_list_sorted, \
                element_list_sorted, bond_list_sorted = self.SetAtomId(
                    xyzList, elementList, bondMatrix)

            # interpret
            __json_atom_position_res2 = self.__json_atom_position(
                atomNo, element_list_sorted, xyz_list_sorted)
            # set
            atom_block_sorted = __json_atom_position_res2.get('atomDetails')
            xyz_list_sorted = __json_atom_position_res2.get('xyzList')
            xyz_center_list_sorted = __json_atom_position_res2.get(
                'xyzCenterList')

            # bond block
            __json_atom_bondblock_res2 = self.__json_atom_bondblock(
                atomNo, bondNo, element_list_sorted, bond_list_sorted)
            # set
            bond_block_sorted = __json_atom_bondblock_res2.get('bondBlock')

            # atomic number sorted
            element_atomic_number = self.arrange_prop(
                _elementAtomicNumber, atom_id_conversion[:, 0])

            # update properties
            # name
            mat_name = propDict.get('IUPAC Name')
            if mat_name is None:
                mat_name = 'NULL'

            # formula
            mat_formula = propDict.get('Molecular Formula')
            if mat_formula is None:
                # create mat formula
                mat_formula = 'NULL'

            # molecular weight
            mat_mass = propDict.get('Molecular Weight')
            if mat_mass is None:
                # mat molecular weight/mass
                mat_mass = 0
            else:
                mat_mass = float(mat_mass)

            # res
            res = {
                'header_block': '',
                'counts_line': '',
                'mat_structure': structure_type,
                'mat_cid': cid,
                'mat_name': mat_name,
                'mat_formula': mat_formula,
                'mat_mass': mat_mass,
                'atom_numbers': atomNo,
                'atom_elements': element_list_sorted,
                'atom_atomic_number': element_atomic_number,
                'atom_block': atom_block_sorted,
                'bond_numbers': bondNo,
                'bond_block': bond_block_sorted,
                'bond_list': bond_list_sorted,
                'xyz_list': xyz_list_sorted,
                'xyz_center_list': xyz_center_list_sorted,
                'compound_properties': propDict,
                'mat_info_origin': origin_info
            }

            return res

        except Exception as e:
            raise Exception(e)

    def __json_parser_id(self, data):
        '''
        return cid 
        '''
        # *** id node
        node = data['id']['cid']
        return node

    def __json_parser_atoms(self, data):
        '''
        return atom id, atom atomic number
        '''
        # aid (atom id)
        _aid = data['aid']
        # element (atomic number)
        _elementAtomicNumber = data['element']
        return _aid, _elementAtomicNumber

    def __json_parser_bonds(self, data):
        '''
        return atom1 id, atom2 is, bond type
        '''
        # aid1 (atom id)
        _aid1 = data['aid1']
        # aid2 (atom id)
        _aid2 = data['aid2']
        # order
        _order = data['order']

        # transform
        bondMatrix = np.array([_aid1, _aid2, _order])
        bondMatrix = np.transpose(bondMatrix)

        return _aid1, _aid2, _order, bondMatrix

    def __json_parser_coords(self, data):
        '''
        return coordination type, coordination id, xyzList
        '''
        data = data[0]
        _coords_type = data['type']  # list
        _coords_aid = data['aid']  # list
        _coords_conformers = data['conformers'][0]
        _x = _coords_conformers.get('x')
        # atomNo
        atomNo = len(_x)
        _y = _coords_conformers.get('y')

        _z = _coords_conformers.get('z') if _coords_conformers.get(
            'z') is not None else np.zeros(atomNo)
        # transform
        xyzList = np.array([_x, _y, _z])
        xyzList = np.transpose(xyzList)

        # 2d/3d structure
        structureType = "2d" if np.count_nonzero(_z) == 0 else '3d'

        return _coords_type, _coords_aid, xyzList, structureType

    def __json_parser_charge(self, data):
        # *** id node
        node = data
        return node

    def __json_parser_props(self, data):
        '''
        return a list of all properties
        '''
        res = {}
        for i in range(len(data)):
            _keySet = data[i]['urn']['label']
            _valueSet = list(data[i]['value'].values())
            res[str(_keySet)] = _valueSet[0]

        return res

    def __json_parser_count(self, data):
        '''
        return a list of all count
        '''
        res = {}
        for key, value in data.items():
            res[key] = value
        return res

    def __json_atom_analyzer(self, atomNo, bondNo, elementList, xyzList, bondMatrix):
        '''
        define xyz list and bond list
        '''
        # vars
        # *** detail about atoms
        atomDetails = []

        # element symbols
        atomList = elementList

        # atoms position
        for i in range(atomNo):
            _atomRow = xyzList[i, :]
            # position
            _x = float(_atomRow[0])
            _y = float(_atomRow[1])
            _z = float(_atomRow[2])
            # name
            _name = atomList[i]

            # atom info
            atom = {
                'id': i,
                'symbol': _name,
                'index': i,
                'x': _x,
                'y': _y,
                'z': _z,
                'position': {
                    'x': _x,
                    'y': _y,
                    'z': _z
                },
                'xyz': [_x, _y, _z]
            }
            # save atom info
            atomDetails.append(atom)

        # object base
        objectBaseCoordinate = Structure.CenterPoints(xyzList)
        # print(f"objectBaseCoordinate: {objectBaseCoordinate}")

        # move to the center [0,0,0]
        xyzCenterList, movingCoordinate = Structure.CenterObject(
            xyzList, objectBaseCoordinate)
        # print(f"xyzCenterList: {xyzCenterList} \n movingCoordinate: \n {movingCoordinate}")

        # bond analysis
        bondBlock = []
        atomBonds = []

        # atomNo: *** atom id in the structure ***
        for i in range(atomNo):
            # name
            _nameAtom1 = str(atomList[i])
            # index [not duplicated element]
            # _nameAtom1Index = atomList.index(_nameAtom1)
            _nameAtom1Index = i+1
            # check bond type
            for j in range(bondNo):
                _bondRow = bondMatrix[j, :]
                _nameAtomBondRow = int(_bondRow[0])
                # check bond exist
                if _nameAtom1Index == _nameAtomBondRow:
                    # atom 2 index
                    _indexAtom2 = int(_bondRow[1])
                    # atom 2 name
                    _nameAtom2 = str(atomList[_indexAtom2-1])
                    # bound type
                    _bondType = int(_bondRow[2])
                    # str bond
                    _bondName = _nameAtom1 + '-' + _nameAtom2
                    # str id bond
                    _bondId = str(_nameAtom1Index) + '-' + str(_indexAtom2)
                    # atom bond
                    atomBonds.append(
                        (_indexAtom2, _nameAtom2, _bondName, _bondType, _bondId))

            if len(atomBonds) > 0:
                # save
                _bondList = {
                    'id': _nameAtom1Index,
                    'symbol': _nameAtom1,
                    'bonds': atomBonds
                }
                bondBlock.append(_bondList)

            # reset
            atomBonds = []
            _bondList = {}

        # result
        return atomDetails, bondBlock, xyzCenterList

    def __json_atom_position(self, atomNo, elementList, xyzList):
        '''
        set mat position into the center of cartesian coordination
        '''
        try:
            # vars
            # *** detail about atoms
            atomDetails = []

            # element symbols
            atomList = elementList

            # atoms position
            for i in range(atomNo):
                _atomRow = xyzList[i, :]
                # position
                _x = float(_atomRow[0])
                _y = float(_atomRow[1])
                _z = float(_atomRow[2])
                # name
                _name = atomList[i]

                # atom info
                atom = {
                    'id': i,
                    'symbol': _name,
                    'index': i,
                    'x': _x,
                    'y': _y,
                    'z': _z,
                    'position': {
                        'x': _x,
                        'y': _y,
                        'z': _z
                    },
                    'xyz': [_x, _y, _z]
                }
                # save atom info
                atomDetails.append(atom)

            # object base
            objectBaseCoordinate = Structure.CenterPoints(xyzList)

            # move to the center [0,0,0]
            xyzCenterList, movingCoordinate = Structure.CenterObject(
                xyzList, objectBaseCoordinate)

            # res
            res = {
                'atomDetails': atomDetails,
                'elementList': elementList,
                'xyzList': xyzList,
                'xyzCenterList': xyzCenterList,
                'movingCoordinate': movingCoordinate
            }

            # res
            return res
        except Exception as e:
            raise Exception(e)

    def __json_atom_bondblock(self, atomNo, bondNo, elementList, bondMatrix):
        '''
        build bond block which is used by mat view to display mat structure
        '''
        try:
            # element symbols
            atomList = elementList

            # bond analysis
            bondBlock = []
            atomBonds = []

            # atomNo: *** atom id in the structure ***
            for i in range(atomNo):
                # name
                _nameAtom1 = str(atomList[i])
                # index [not duplicated element]
                _nameAtom1Index = int(i+1)
                # check bond type
                for j in range(bondNo):
                    # name
                    _bondRow = bondMatrix[j, :]
                    # atom id
                    _nameAtomBondRow = int(_bondRow[0])
                    # check bond exist
                    if _nameAtom1Index == _nameAtomBondRow:
                        # atom 2 index
                        _indexAtom2 = int(_bondRow[1])
                        # atom 2 name
                        _nameAtom2 = str(atomList[_indexAtom2-1])
                        # bound type
                        _bondType = int(_bondRow[2])
                        # str bond
                        _bondName = _nameAtom1 + '-' + _nameAtom2
                        # str id bond
                        _bondId = str(_nameAtom1Index) + '-' + str(_indexAtom2)
                        # atom bond
                        atomBonds.append(
                            (_indexAtom2, _nameAtom2, _bondName, _bondType, _bondId))

                if len(atomBonds) > 0:
                    # save
                    _bondList = {
                        'id': _nameAtom1Index,
                        'symbol': _nameAtom1,
                        'bonds': atomBonds
                    }
                    bondBlock.append(_bondList)

                # reset
                atomBonds = []
                _bondList = {}

            # res
            res = {
                'elementList': elementList,
                'bondMatrix': bondMatrix,
                'bondBlock': bondBlock
            }

            # res
            return res
        except Exception as e:
            Exception(e)

    def __find_element_symbol(self, atomic_numbers):
        '''
        return element symbols
        '''
        el = Element()
        # loop
        res = el.find_atom_by_property('AtomicNumber', atomic_numbers)

        # interpret
        atomList = [item['symbol'] for item in res]

        # res
        return res, atomList

    def __var_finder(data):
        '''
        find variables in a sdf file (single)
        '''
        res = re.findall(
            r"(\>\s*\<(.*)\s*\>\s*((.*\n)([^\<\>\$\$\$\$])*))", data, re.M)
        return res

    def __var_analyzer(data):
        '''
        make a dict of all properties
        '''
        res = {}
        dataSize = len(data)
        for i in range(dataSize):
            # name
            varName = data[i][1]
            # val
            varVal = (data[i][2]).split('\n')
            varVal = [item.strip() for item in varVal]
            varVal = list(filter(None, varVal))

            # check
            if len(varVal) == 1:
                _varVal = str(varVal[0])
            else:
                _varVal = varVal
            # set
            _keySet = str(varName)
            _valueSet = _varVal
            # res
            res[_keySet] = _valueSet
            # reset

        return res

    def SetAtomId(self, xyzList, elementList, bondList):
        '''
        set atom id with respect to their position in a 3d frame

        args:
            xyzList: list of element (atom) position
            elementList: list of elements
            bondList: list of bond ids and types
            robs: observer position, default [10,0,0]
        '''
        try:
            # robs position
            robs = OBS_POSITIONS
            # res
            disRes = []
            # calculate distance
            for i in range(len(xyzList)):
                _dis = np.linalg.norm(
                    np.array(robs) - np.array(xyzList[i]))
                _idOld = int(i+1)
                _symbol = elementList[i]
                disRes.append([_idOld, _dis, _symbol])

            # sort
            sortRes = sorted(disRes, key=lambda l: l[1], reverse=True)

            # reverse [old id, distance, symbol]
            sortRes.reverse()

            # save new id [old id, new id, distance, symbol]
            for j in range(len(sortRes)):
                _idNew = int(j+1)
                sortRes[j].insert(1, _idNew)

            # sorted id (old id, new id)
            idConversion = []
            for j in range(len(sortRes)):
                idConversion.append([sortRes[j][0], sortRes[j][1]])

            # build xyzList and elementlist with respect to the new ids
            xyzListSorted = []
            elementListSorted = []
            matPosition = {}
            for j in range(len(sortRes)):
                _idOld = sortRes[j][0]
                _idNew = sortRes[j][1]
                _dis = sortRes[j][2]
                _val1 = xyzList[_idOld-1]
                _val2 = elementList[_idOld-1]
                # set
                xyzListSorted.append(_val1)
                elementListSorted.append(_val2)
                # all
                matPosition[str(_idNew)] = [_idOld, _idNew, _dis, _val1, _val2]

            # set
            xyzListSorted = np.array(xyzListSorted)

            # bond id conversion
            bondList = np.array(bondList, dtype='i')
            # size
            bondRow = bondList.shape[0]
            # column 0
            bondColumn0 = np.array(bondList[:, 0])
            # column 1
            bondColumn1 = np.array(bondList[:, 1])
            # column 2
            bondColumn2 = np.array(bondList[:, 2])

            # replace
            for i in range(bondRow):
                # column 0
                _v0 = bondColumn0[i]
                # find new value
                _v1 = [item[1] for item in idConversion if item[0] == _v0][0]
                bondColumn0[i] = _v1

                # column 1
                _v2 = bondColumn1[i]
                # find new value
                _v3 = [item[1] for item in idConversion if item[0] == _v2][0]
                bondColumn1[i] = _v3

            # build bond list with new ids
            bondListSorted = np.zeros_like(bondList)
            bondListSorted[:, 0] = bondColumn0
            bondListSorted[:, 1] = bondColumn1
            bondListSorted[:, 2] = bondColumn2

            # set
            idConversion = np.array(idConversion)

            return matPosition, idConversion, xyzListSorted, elementListSorted, bondListSorted

        except Exception as e:
            Exception(e)

    def arrange_prop(self, property_value_list, atom_index_list):
        '''
        arrange a new list with respect to the index list

        args:
            property_value_list: such as atomic number
            atom_index_list: atom id between 1 and ...
        '''
        try:
            # convert atom id to list id
            index_list = np.array(atom_index_list) - 1
            property_value_list_sorted = np.take(
                property_value_list, index_list)
            # res
            return property_value_list_sorted
        except Exception as e:
            Exception(e)

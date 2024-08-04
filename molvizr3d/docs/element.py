# PUBCHEM ELEMENTS
# -----------------

# import libs
import pandas as pd
import os
import numpy as np


class Element():
    '''
    pub chem elements (periodic table of the elements)
    '''

    _elementsource = ''
    _ele = ''

    def __init__(self, atom_symbol=''):
        self.elementsource = self.__load_elements()
        self._ele = atom_symbol

    def __call__(self, ):
        print(
            'element object such as Carbon (C), get properties by calling atom_properties()')

    @property
    def elementsource(self):
        return self._elementsource

    @elementsource.setter
    def elementsource(self, value):
        self._elementsource = value

    @property
    def ele(self):
        return self._ele

    @ele.setter
    def ele(self, value):
        pass

    def __load_elements(self):
        '''
        load elements from a csv file
        '''
        try:
            # data file
            dataFile = 'PubChemElements_all.csv'
            # abs path
            pathAbs = os.path.abspath(os.path.dirname(__file__))
            # relative path to database file
            dataPathDirRel = '../data'
            # database file
            dataPath = os.path.join(pathAbs, dataPathDirRel, dataFile)

            with open(dataPath, 'rb') as f:
                df = pd.read_csv(f)

            return df

        except Exception as e:
            raise Exception(e)

    def properties(self):
        '''
        return all atom properties
        '''
        df = self.elementsource
        filt = df['Symbol'] == str(self._ele).strip()
        return df.loc[filt]

    def find_atom(self, atom_symbol):
        '''
        return the selected atom properties
        '''
        df = self.elementsource
        filt = df['Symbol'] == str(atom_symbol).strip()
        return df.loc[filt]

    def atom_properties(self, atom_symbol, atom_properties=[]):
        '''
        find desired atom properties

        args:
            atom_symbol: atom symbol
            atom_properties: a list of desired properties

        return:
            a dict of property
        '''
        try:
            # property size
            atomPropertySize = len(atom_properties)
            # res
            resDict = {}

            # check
            if atomPropertySize == 0:
                raise Exception('property list is empty.')

            df = self.elementsource
            filt = df['Symbol'] == str(atom_symbol).strip()

            # check
            if atomPropertySize == 1:
                rowRes = df.loc[filt][atom_properties[0]]
                rowRes = rowRes.to_numpy()[0]
                # add prop
                resDict[str(atom_properties[0])] = rowRes
            elif atomPropertySize > 1:
                rowRes = df.loc[filt][atom_properties]
                rowRes = rowRes.to_numpy()[0]
                # add prop
                for i in range(atomPropertySize):
                    resDict[str(atom_properties[i])] = rowRes[i]
            # res
            return resDict

        except Exception as e:
            raise Exception(e)

    def find_atom_by_property(self, atom_property_name, atom_property_value=[]):
        '''
        find desired atom properties

        args:
            atom_property_name: atom property such as atomic number
            atom_property_value: carbon = 6

        return:
            string
        '''
        try:
            # size
            atom_property_value_size = len(atom_property_value)
            atom_property_name_size = len(str(atom_property_name))

            # res
            res = []

            # check
            if atom_property_name_size == 0 or atom_property_value_size == 0:
                raise Exception('args error.')

            df = self.elementsource

            if atom_property_value_size == 1:
                filt = df[str(atom_property_name)] == atom_property_value[0]
                _rowRes = df.loc[filt][['Symbol']]
                _rowRes = _rowRes.to_numpy()[0]
                # check
                if len(_rowRes) > 0:
                    # add prop
                    _val = {'symbol': _rowRes[0], str(
                        atom_property_name): atom_property_value[0]}
                    res.append(_val)
                else:
                    raise Exception('element not found.')
            elif atom_property_value_size > 1:
                for i in range(atom_property_value_size):
                    filt = df[str(atom_property_name)
                              ] == atom_property_value[i]
                    _rowRes = df.loc[filt][['Symbol']]
                    _rowRes = _rowRes.to_numpy()[0]
                    # add
                    _val = {'symbol': _rowRes[0], str(
                        atom_property_name): atom_property_value[i]}
                    # save
                    res.append(_val)

            # res
            return res

        except Exception as e:
            raise Exception(e)

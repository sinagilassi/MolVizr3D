# NETWORK
# -------

# import packages/modules
import numpy as np
import networkx as nx


class Network():

    # Define functional groups
    functional_groups = {
        'hydroxyl': {'nodes': ['O', 'H'], 'edges': [('O', 'H')]}
    }

    def __init__(self, atomElements, atomBonds, xyzList, xyzCenterList):
        self.atomElements = atomElements
        # bond block (info)
        self.atomBonds = atomBonds
        self.xyzList = xyzList
        self.xyzCenterList = xyzCenterList

    def check_functional_groups(self, functional_groups):
        '''
        Check functional groups in a compound

        Parameters
        ----------
        functional_groups : list
            list of functional groups

        Returns
        -------
        res : dict
            a list of all count
        '''
        # create graph
        G = self.create_graph()
        # res
        return G

    def create_graph(self):
        '''
        Check functional groups in a compound

        Parameters
        ----------
        functional_groups : list
            list of functional groups

        Returns
        -------
        res : dict
            a list of all count
        '''
        # atom no
        atomNo = len(self.xyzList)
        # bond no
        bondNo = len(self.atomBonds)

        # Create a graph from atoms and bonds
        G = nx.Graph()

        # *** atom visualization
        for i in range(atomNo):
            # xyz
            _atom1X = self.xyzList[i, 0]
            _atom1Y = self.xyzList[i, 1]
            _atom1Z = self.xyzList[i, 2]
            _atom1XYZ = [_atom1X, _atom1Y, _atom1Z]

            # color
            # atom id
            _atomId = int(i+1)
            # symbol
            _atomSymbol = str(self.atomElements[i]).strip()
            # atom mark
            atomMark = str(_atomSymbol) + str(_atomId)

            # add node
            G.add_node(_atomId, symbol=_atomSymbol,
                       x=_atom1X, y=_atom1Y, z=_atom1Z, xyz=_atom1XYZ)

        # reset
        i = 0

        # *** bond visualization
        # *** using bond block
        for i in range(bondNo):
            # atom id
            _atom1Id = int(self.atomBonds[i]['id']) - 1
            # atom symbol
            _atom1Symbol = self.atomBonds[i]['symbol']
            # atom color
            _atom1Color = self.set_color(_atom1Symbol)
            # atom bond list
            _atom1BondList = self.atomBonds[i]['bonds']
            # number of bonds
            atom1BondSize = len(_atom1BondList)

            _atom1X = self.xyzList[_atom1Id, 0]
            _atom1Y = self.xyzList[_atom1Id, 1]
            _atom1Z = self.xyzList[_atom1Id, 2]
            _atom1XYZ = [_atom1X, _atom1Y, _atom1Z]

            # draw bond
            if atom1BondSize > 0:
                for j in range(atom1BondSize):
                    # atom [2] id
                    _atom2Id = int(_atom1BondList[j][0]) - 1
                    # atom [2] symbol
                    _atom2Symbol = _atom1BondList[j][1]
                    # atom [1] - atom [2] bond type
                    _bondType = int(_atom1BondList[j][3])

                    # xyz
                    _atom2X = self.xyzList[_atom2Id, 0]
                    _atom2Y = self.xyzList[_atom2Id, 1]
                    _atom2Z = self.xyzList[_atom2Id, 2]
                    _atom2XYZ = [_atom2X, _atom2Y, _atom2Z]

                    # add edges
                    G.add_edge(_atom1Id, _atom2Id)

        # res
        return G

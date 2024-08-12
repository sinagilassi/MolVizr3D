# NETWORK
# -------

# import packages/modules
import pandas as pd
import networkx as nx
from networkx.algorithms import isomorphism
# local
from .chemgraphs import ChemGraphs


class Network(ChemGraphs):

    def __init__(self, atomElements, atomBonds, xyzList, xyzCenterList, atomBonds1d):
        self.atomElements = atomElements
        # bond block (info)
        self.atomBonds = atomBonds
        self.xyzList = xyzList
        self.xyzCenterList = xyzCenterList

        # 1d vector of atom bonds
        self.atomBonds1d = atomBonds1d

        # TODO: super
        ChemGraphs.__init__(self)

        # functional group list
        self.function_group_list = {
            'hydroxyl': [self.graph_hydroxyl()],
            'carbonyl': [self.graph_carbonyl()],
            'carboxyl': [self.graph_carboxyl()],
            'ether': [self.graph_ether()],
            'alcohols': [self.graph_primary_alcohol(), self.graph_secondary_alcohol(), self.graph_tertiary_alcohol()],
            'alkane': [self.graph_alkane()],
            'alkene': [self.graph_alkene()],
            'alkyne': [self.graph_alkyne()],
            'arene': [self.graph_arene()],
            'aldehyde': [self.graph_aldehyde()],
            'ketone': [self.graph_ketone()],
            'carboxylic_acid': [self.graph_carboxylic_acid()],
            'ester': [self.graph_ester()],
            'amide': [self.graph_primary_amide(), self.graph_secondary_amide(), self.graph_tertiary_amide()],
            'primary_amide': [self.graph_primary_amide()],
            'secondary_amide': [self.graph_secondary_amide()],
            'tertiary_amide': [self.graph_tertiary_amide()],
            'amine': [self.graph_primary_amine(), self.graph_secondary_amine(), self.graph_tertiary_amine()],
            'primary_amine': [self.graph_primary_amine()],
            'secondary_amine': [self.graph_secondary_amine()],
            'tertiary_amine': [self.graph_tertiary_amine()],
            'nitrile': [self.graph_nitrile()],
            'thiol': [self.graph_thiol()],
            'alkyl_halids': [self.graph_alkyl_halide('F'), self.graph_alkyl_halide('Cl'), self.graph_alkyl_halide('Br'),
                             self.graph_alkyl_halide('I'),
                             self.graph_primary_alkyl_halide(
                                 'F'), self.graph_primary_alkyl_halide('Cl'),
                             self.graph_primary_alkyl_halide('Br'), self.graph_primary_alkyl_halide('I')]
        }

    def check_functional_groups(self, functional_groups=[]):
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
        # check functional group
        if len(functional_groups) == 0:
            functional_groups = list(self.function_group_list.keys())
        # create graph
        G = self.create_graph()
        # check functional groups
        res = self.check_functional_group(G, functional_groups)
        # res
        return res

    def check_functional_group(self, G, function_groups):
        '''
        Check a functional group exists in a compound

        Parameters
        ----------
        G : graph
            graph
        function_groups : list[str]
            functional group name like hydroxyl

        Returns
        -------
        res : dict
            a list of all count
        '''
        def node_match(n1, n2):
            '''
            Define a custom node_match function to ensure element matching
            '''
            return n1['symbol'] == n2['symbol']

        def edge_match(e1, e2):
            '''
            Define a custom edge_match function to ensure bond type matching
            '''
            return e1['type'] == e2['type']

        # res
        res_match = []

        # for each functional group
        for item in function_groups:
            # ! check functional exists in the list
            if item in self.function_group_list:
                # get a list of graphs
                function_group_graphs = self.function_group_list.get(item)

                # flag to track if functional group is found
                fg_found_any = False

                # graph function list
                for _fn in function_group_graphs:
                    # Create a GraphMatcher object for a functional group
                    fg_matcher = isomorphism.GraphMatcher(
                        G, _fn, node_match=node_match, edge_match=edge_match)

                    # Check if a functional group is in the main graph
                    fg_found = fg_matcher.subgraph_is_isomorphic()

                    if fg_found:
                        # print(f"{item} found in the molecule!")
                        fg_found_any = True
                        # res
                        res_match.append({
                            'function_group': item,
                            'result': fg_found
                        })
                        # break
                        break
                # check
                if not fg_found_any:
                    # print(f"{item} not found in the molecule.")
                    # res
                    res_match.append({
                        'function_group': item,
                        'result': False
                    })

        # res
        return res_match

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
        bondNo = len(self.atomBonds1d)

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

            # add node with cooridination
            # G.add_node(_atomId, symbol=_atomSymbol,
            #            x=_atom1X, y=_atom1Y, z=_atom1Z, xyz=_atom1XYZ)

            G.add_node(_atomId, symbol=_atomSymbol)

        # *** using bond block

        for i in range(bondNo):
            # atom 1
            _id1 = self.atomBonds1d[i]['id1']
            _symbol1 = self.atomBonds1d[i]['symbol1']
            # atom 2
            _id2 = self.atomBonds1d[i]['id2']
            _symbol2 = self.atomBonds1d[i]['symbol2']
            # bond type
            _bondType = self.atomBonds1d[i]['bond_type']
            # bond symbol
            _bondSymbol = self.atomBonds1d[i]['bond_symbol']

            # add edge
            G.add_edge(_id1, _id2, symbol=_bondSymbol, type=_bondType)

        # res
        return G

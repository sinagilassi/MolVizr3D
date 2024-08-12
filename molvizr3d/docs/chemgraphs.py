# DEFINE FUNCTIONAL GROUPS
# -----------------------

# import packages/modules
import networkx as nx


class ChemGraphs():

    def __init__(self):
        pass

    def graph_hydroxyl(self):
        '''
        Create a graph for hydroxyl
        '''
        # OH graph
        OH = nx.Graph()
        # O
        OH.add_node(1, symbol='O')
        # H
        OH.add_node(2, symbol='H')
        # bond
        OH.add_edge(1, 2, symbol='OH', type=1)

        # res
        return OH

    def graph_carbonyl(self):
        '''
        Create a graph for carbonyl
        '''
        # C=O graph
        CO = nx.Graph()
        # C
        CO.add_node(1, symbol='C')
        # O
        CO.add_node(2, symbol='O')
        # bond
        CO.add_edge(1, 2, symbol='CO', type=2)

        # res
        return CO

    def graph_carboxyl(self):
        '''
        Create a graph for carboxyl
        '''
        # C=O graph
        CO = nx.Graph()
        # C
        CO.add_node(1, symbol='C')
        # O
        CO.add_node(2, symbol='O')
        # O
        CO.add_node(3, symbol='O')
        # H
        CO.add_node(4, symbol='H')
        # bond
        CO.add_edge(1, 2, symbol='CO', type=2)
        # bond
        CO.add_edge(1, 3, symbol='CO', type=1)
        # bond
        CO.add_edge(3, 4, symbol='OH', type=1)

        # res
        return CO

    def graph_ether(self):
        '''
        Create a graph for ether
        '''
        # COC graph
        COC = nx.Graph()
        # C
        COC.add_node(1, symbol='C')
        # O
        COC.add_node(2, symbol='O')
        # C
        COC.add_node(3, symbol='C')
        # bond
        COC.add_edge(1, 2, symbol='CO', type=1)
        # bond
        COC.add_edge(2, 3, symbol='OC', type=1)

        # res
        return COC

    def graph_alkane(self):
        '''
        Create a graph for alkane
        '''
        # C-C graph
        CC = nx.Graph()
        # C
        CC.add_node(1, symbol='C')
        # C
        CC.add_node(2, symbol='C')
        # bond
        CC.add_edge(1, 2, symbol='CC', type=1)

        # res
        return CC

    def graph_alkene(self):
        '''
        Create a graph for alkene
        '''
        # C=C graph
        CC = nx.Graph()
        # C
        CC.add_node(1, symbol='C')
        # C
        CC.add_node(2, symbol='C')
        # bond
        CC.add_edge(1, 2, symbol='CC', type=2)

        # res
        return CC

    def graph_alkyne(self):
        '''
        Create a graph for alkyne
        '''
        # C≡C graph
        CC = nx.Graph()
        # C
        CC.add_node(1, symbol='C')
        # C
        CC.add_node(2, symbol='C')
        # bond
        CC.add_edge(1, 2, symbol='CC', type=3)

        # res
        return CC

    def graph_arene(self):
        '''
        Create a graph for arene
        '''
        # C6H5 graph (simplified)
        C6H5 = nx.Graph()
        # C
        C6H5.add_node(1, symbol='C')
        # C
        C6H5.add_node(2, symbol='C')
        # C
        C6H5.add_node(3, symbol='C')
        # C
        C6H5.add_node(4, symbol='C')
        # C
        C6H5.add_node(5, symbol='C')
        # C
        C6H5.add_node(6, symbol='C')
        # bond
        C6H5.add_edge(1, 2, symbol='CC', type=1)
        # bond
        C6H5.add_edge(2, 3, symbol='CC', type=2)
        # bond
        C6H5.add_edge(3, 4, symbol='CC', type=1)
        # bond
        C6H5.add_edge(4, 5, symbol='CC', type=2)
        # bond
        C6H5.add_edge(5, 6, symbol='CC', type=1)
        # bond
        C6H5.add_edge(6, 1, symbol='CC', type=2)

        # res
        return C6H5

    def graph_aldehyde(self):
        '''
        Create a graph for aldehyde
        '''
        # CHO graph
        CHO = nx.Graph()
        # C
        CHO.add_node(1, symbol='C')
        # H
        CHO.add_node(2, symbol='H')
        # O
        CHO.add_node(3, symbol='O')
        # C
        CHO.add_node(4, symbol='C')
        # bond
        CHO.add_edge(1, 2, symbol='CH', type=1)
        # bond
        CHO.add_edge(1, 3, symbol='CO', type=2)
        # bond
        CHO.add_edge(1, 4, symbol='CC', type=1)

        # res
        return CHO

    def graph_ketone(self):
        '''
        Create a graph for ketone
        '''
        # COC graph
        COC = nx.Graph()
        # C
        COC.add_node(1, symbol='C')
        # O
        COC.add_node(2, symbol='O')
        # C
        COC.add_node(3, symbol='C')
        # C
        COC.add_node(4, symbol='C')
        # bond
        COC.add_edge(1, 2, symbol='CO', type=2)
        # bond
        COC.add_edge(1, 3, symbol='CC', type=1)
        # bond
        COC.add_edge(1, 4, symbol='CC', type=1)

        # res
        return COC

    def graph_carboxylic_acid(self):
        '''
        Create a graph for carboxylic acid
        '''
        # COOH graph
        COOH = nx.Graph()
        # C
        COOH.add_node(1, symbol='C')
        # O
        COOH.add_node(2, symbol='O')
        # O
        COOH.add_node(3, symbol='O')
        # H
        COOH.add_node(4, symbol='H')
        # C
        COOH.add_node(5, symbol='C')
        # bond
        COOH.add_edge(1, 2, symbol='CO', type=2)
        # bond
        COOH.add_edge(1, 3, symbol='CO', type=1)
        # bond
        COOH.add_edge(3, 4, symbol='OH', type=1)
        # bond
        COOH.add_edge(1, 5, symbol='CC', type=1)

        # res
        return COOH

    def graph_ester(self):
        '''
        Create a graph for ester
        '''
        # COOC graph
        COOC = nx.Graph()
        # C
        COOC.add_node(1, symbol='C')
        # O
        COOC.add_node(2, symbol='O')
        # O
        COOC.add_node(3, symbol='O')
        # C
        COOC.add_node(4, symbol='C')
        # C
        COOC.add_node(5, symbol='C')
        # bond
        COOC.add_edge(1, 2, symbol='CO', type=2)
        # bond
        COOC.add_edge(1, 3, symbol='CO', type=1)
        # bond
        COOC.add_edge(3, 4, symbol='OC', type=1)
        # bond
        COOC.add_edge(1, 5, symbol='CC', type=1)

        # res
        return COOC

    def graph_amide(self):
        '''
        Create a graph for amide
        '''
        # REVIEW
        # CONH2 graph
        CONH2 = nx.Graph()
        # C
        CONH2.add_node(1, symbol='C')
        # O
        CONH2.add_node(2, symbol='O')
        # N
        CONH2.add_node(3, symbol='N')
        # H
        CONH2.add_node(4, symbol='H')
        # H
        CONH2.add_node(5, symbol='H')
        # C
        CONH2.add_node(6, symbol='C')
        # bond
        CONH2.add_edge(1, 2, symbol='CO', type=2)
        # bond
        CONH2.add_edge(1, 3, symbol='CN', type=1)
        # bond
        CONH2.add_edge(3, 4, symbol='NH', type=1)
        # bond
        CONH2.add_edge(3, 5, symbol='NH', type=1)
        # bond
        CONH2.add_edge(1, 6, symbol='CC', type=1)

        # res
        return CONH2

    def graph_primary_amide(self):
        '''
        Create a graph for a primary amide
        '''
        # CONH2 graph
        CONH2 = nx.Graph()
        # C
        CONH2.add_node(1, symbol='C')
        # O
        CONH2.add_node(2, symbol='O')
        # N
        CONH2.add_node(3, symbol='N')
        # H
        CONH2.add_node(4, symbol='H')
        # H
        CONH2.add_node(5, symbol='H')
        # C
        CONH2.add_node(6, symbol='C')
        # bond
        CONH2.add_edge(1, 2, symbol='CO', type=2)
        # bond
        CONH2.add_edge(1, 3, symbol='CN', type=1)
        # bond
        CONH2.add_edge(3, 4, symbol='NH', type=1)
        # bond
        CONH2.add_edge(3, 5, symbol='NH', type=1)
        # bond
        CONH2.add_edge(1, 6, symbol='CC', type=1)

        # res
        return CONH2

    def graph_secondary_amide(self):
        '''
        Create a graph for a secondary amide
        '''
        # CONH2 graph
        CONH2 = nx.Graph()
        # C
        CONH2.add_node(1, symbol='C')
        # O
        CONH2.add_node(2, symbol='O')
        # N
        CONH2.add_node(3, symbol='N')
        # H
        CONH2.add_node(4, symbol='H')
        # H
        CONH2.add_node(5, symbol='C')
        # C
        CONH2.add_node(6, symbol='C')
        # bond
        CONH2.add_edge(1, 2, symbol='CO', type=2)
        # bond
        CONH2.add_edge(1, 3, symbol='CN', type=1)
        # bond
        CONH2.add_edge(3, 4, symbol='NH', type=1)
        # bond
        CONH2.add_edge(3, 5, symbol='NC', type=1)
        # bond
        CONH2.add_edge(1, 6, symbol='CC', type=1)

        # res
        return CONH2

    def graph_tertiary_amide(self):
        '''
        Create a graph for a tertiary amide
        '''
        # CONH2 graph
        CONH2 = nx.Graph()
        # C
        CONH2.add_node(1, symbol='C')
        # O
        CONH2.add_node(2, symbol='O')
        # N
        CONH2.add_node(3, symbol='N')
        # H
        CONH2.add_node(4, symbol='C')
        # H
        CONH2.add_node(5, symbol='C')
        # C
        CONH2.add_node(6, symbol='C')
        # bond
        CONH2.add_edge(1, 2, symbol='CO', type=2)
        # bond
        CONH2.add_edge(1, 3, symbol='CN', type=1)
        # bond
        CONH2.add_edge(3, 4, symbol='NC', type=1)
        # bond
        CONH2.add_edge(3, 5, symbol='NC', type=1)
        # bond
        CONH2.add_edge(1, 6, symbol='CC', type=1)

        # res
        return CONH2

    def graph_primary_amine(self):
        '''
        Create a graph for primary amine
        '''
        # NH2 graph
        NH2 = nx.Graph()
        # N
        NH2.add_node(1, symbol='N')
        # H
        NH2.add_node(2, symbol='H')
        # H
        NH2.add_node(3, symbol='H')
        # C
        NH2.add_node(4, symbol='C')
        # bond
        NH2.add_edge(1, 2, symbol='NH', type=1)
        # bond
        NH2.add_edge(1, 3, symbol='NH', type=1)
        # bond
        NH2.add_edge(1, 4, symbol='NC', type=1)

        # res
        return NH2

    def graph_secondary_amine(self):
        '''
        Create a graph for secondary amine
        '''
        # NH graph
        NH = nx.Graph()
        # N
        NH.add_node(1, symbol='N')
        # H
        NH.add_node(2, symbol='H')
        # C
        NH.add_node(3, symbol='C')
        # C
        NH.add_node(4, symbol='C')
        # bond
        NH.add_edge(1, 2, symbol='NH', type=1)
        # bond
        NH.add_edge(1, 3, symbol='NC', type=1)
        # bond
        NH.add_edge(1, 4, symbol='NC', type=1)

        # res
        return NH

    def graph_tertiary_amine(self):
        '''
        Create a graph for tertiary amine
        '''
        # N graph
        N = nx.Graph()
        # N
        N.add_node(1, symbol='N')
        # C
        N.add_node(2, symbol='C')
        # C
        N.add_node(3, symbol='C')
        # C
        N.add_node(4, symbol='C')
        # bond
        N.add_edge(1, 2, symbol='NC', type=1)
        # bond
        N.add_edge(1, 3, symbol='NC', type=1)
        # bond
        N.add_edge(1, 4, symbol='NC', type=1)

        # res
        return N

    def graph_nitrile(self):
        '''
        Create a graph for nitrile
        '''
        # CN graph
        CN = nx.Graph()
        # C
        CN.add_node(1, symbol='C')
        # N
        CN.add_node(2, symbol='N')
        # C
        CN.add_node(3, symbol='C')
        # bond
        CN.add_edge(1, 2, symbol='CN', type=3)
        # bond
        CN.add_edge(1, 3, symbol='NC', type=1)

        # res
        return CN

    def graph_thiol(self):
        '''
        Create a graph for thiol
        '''
        # SH graph
        SH = nx.Graph()
        # S
        SH.add_node(1, symbol='S')
        # H
        SH.add_node(2, symbol='H')
        # C
        SH.add_node(3, symbol='C')
        # bond
        SH.add_edge(1, 2, symbol='SH', type=1)
        # bond
        SH.add_edge(1, 3, symbol='SC', type=1)

        # res
        return SH

    def graph_alkyl_halide(self, halid):
        '''
        Create a graph for primary alkyl halide
        '''
        # R-X graph (simplified)
        RX = nx.Graph()
        # C
        RX.add_node(1, symbol='C')
        # X (halogen)
        RX.add_node(2, symbol=f'{halid}')
        # bond
        RX.add_edge(1, 2, symbol=f'C{halid}', type=1)

        # res
        return RX

    def graph_primary_alkyl_halide(self, halid):
        '''
        Create a graph for primary alkyl halide
        '''
        # R-X graph (simplified)
        RX = nx.Graph()
        # C
        RX.add_node(1, symbol='C')
        # C
        RX.add_node(2, symbol='C')
        # H
        RX.add_node(3, symbol='H')
        # H
        RX.add_node(4, symbol='H')
        # H
        RX.add_node(5, symbol='H')
        # X (halogen)
        RX.add_node(6, symbol=f'{halid}')
        # bond
        RX.add_edge(1, 2, symbol='CC', type=1)
        # bond
        RX.add_edge(1, 3, symbol='CH', type=1)
        # bond
        RX.add_edge(1, 4, symbol='CH', type=1)
        # bond
        RX.add_edge(1, 5, symbol='CH', type=1)
        # bond
        RX.add_edge(2, 6, symbol=f'C{halid}', type=1)

        # res
        return RX

    def graph_secondary_alkyl_halide(self):
        '''
        Create a graph for secondary alkyl halide
        '''
        # R-X graph (simplified)
        RX = nx.Graph()
        # C
        RX.add_node(1, symbol='C')
        # C
        RX.add_node(2, symbol='C')
        # C
        RX.add_node(3, symbol='C')
        # H
        RX.add_node(4, symbol='H')
        # H
        RX.add_node(5, symbol='H')
        # X (halogen)
        RX.add_node(6, symbol='X')
        # bond
        RX.add_edge(1, 2, symbol='CC', type=1)
        # bond
        RX.add_edge(1, 3, symbol='CC', type=1)
        # bond
        RX.add_edge(1, 4, symbol='CH', type=1)
        # bond
        RX.add_edge(1, 5, symbol='CH', type=1)
        # bond
        RX.add_edge(2, 6, symbol='CX', type=1)

        # res
        return RX

    def graph_tertiary_alkyl_halide(self):
        '''
        Create a graph for tertiary alkyl halide
        '''
        # R-X graph (simplified)
        RX = nx.Graph()
        # C
        RX.add_node(1, symbol='C')
        # C
        RX.add_node(2, symbol='C')
        # C
        RX.add_node(3, symbol='C')
        # C
        RX.add_node(4, symbol='C')
        # X (halogen)
        RX.add_node(5, symbol='X')
        # bond
        RX.add_edge(1, 2, symbol='CC', type=1)
        # bond
        RX.add_edge(1, 3, symbol='CC', type=1)
        # bond
        RX.add_edge(1, 4, symbol='CC', type=1)
        # bond
        RX.add_edge(1, 5, symbol='CX', type=1)

        # res
        return RX

    def graph_primary_alcohol(self):
        '''
        Create a graph for primary alcohol
        '''
        # R-OH graph (simplified)
        ROH = nx.Graph()
        # C
        ROH.add_node(1, symbol='C')
        # H
        ROH.add_node(2, symbol='H')
        # H
        ROH.add_node(3, symbol='H')
        # O
        ROH.add_node(4, symbol='O')
        # H
        ROH.add_node(5, symbol='H')
        # bond
        ROH.add_edge(1, 2, symbol='CH', type=1)
        # bond
        ROH.add_edge(1, 3, symbol='CH', type=1)
        # bond
        ROH.add_edge(1, 4, symbol='CO', type=1)
        # bond
        ROH.add_edge(4, 5, symbol='OH', type=1)

        # res
        return ROH

    def graph_secondary_alcohol(self):
        '''
        Create a graph for secondary alcohol
        '''
        # R-OH graph (simplified)
        ROH = nx.Graph()
        # C
        ROH.add_node(1, symbol='C')
        # C
        ROH.add_node(2, symbol='C')
        # C
        ROH.add_node(3, symbol='C')
        # H
        ROH.add_node(4, symbol='H')
        # O
        ROH.add_node(5, symbol='O')
        # H
        ROH.add_node(6, symbol='H')
        # bond
        ROH.add_edge(1, 2, symbol='CC', type=1)
        # bond
        ROH.add_edge(1, 3, symbol='CC', type=1)
        # bond
        ROH.add_edge(1, 4, symbol='CH', type=1)
        # bond
        ROH.add_edge(1, 5, symbol='CO', type=1)
        # bond
        ROH.add_edge(5, 6, symbol='OH', type=1)

        # res
        return ROH

    def graph_tertiary_alcohol(self):
        '''
        Create a graph for tertiary alcohol
        '''
        # R-OH graph (simplified)
        ROH = nx.Graph()
        # C
        ROH.add_node(1, symbol='C')
        # C
        ROH.add_node(2, symbol='C')
        # C
        ROH.add_node(3, symbol='C')
        # C
        ROH.add_node(4, symbol='C')
        # O
        ROH.add_node(5, symbol='O')
        # H
        ROH.add_node(6, symbol='H')
        # bond
        ROH.add_edge(1, 2, symbol='CC', type=1)
        # bond
        ROH.add_edge(1, 3, symbol='CC', type=1)
        # bond
        ROH.add_edge(1, 4, symbol='CC', type=1)
        # bond
        ROH.add_edge(1, 5, symbol='CO', type=1)
        # bond
        ROH.add_edge(5, 6, symbol='OH', type=1)

        # res
        return ROH

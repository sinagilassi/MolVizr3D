# DISPLAY A CHEMICAL STRUCTURE
# -----------------------------

# import libs
import os
import math
import numpy as np
import matplotlib.pyplot as plt
# internal
from .observer import Observer


class Vizr3D():
    '''
    visualize chemical structure

    hint:
        xyzList is selected for visualization
    '''

    # properties
    _structure_type = ''

    def __init__(self, atomElements, atomBonds, xyzList, xyzCenterList, robs, tetaNo, phiNo, limits):
        self.atomElements = atomElements
        # bond block (info)
        self.atomBonds = atomBonds
        self.xyzList = xyzList
        self.xyzCenterList = xyzCenterList
        self.robs = robs
        self.tetaNo = tetaNo
        self.phiNo = phiNo
        self.limits = limits

        # set structure type
        structureType, perpendicularAxis, perpendicularVector, XYZ0 = self.StructureAnalyzer()
        self.structure_type = structureType

    @property
    def structure_type(self):
        return self._structure_type

    @structure_type.setter
    def structure_type(self, value):
        self._structure_type = value

    def StructureAnalyzer(self):
        '''
        check geometry structure to determine 2D/3D 

        args:
            xyzList: point coordination 
        '''
        try:
            # array
            xyzList = np.array(self.xyzList)
            # set
            X = xyzList[:, 0]
            Y = xyzList[:, 1]
            Z = xyzList[:, 2]

            XYZLen = len(X)

            # check plane structure (2D, 3D)
            X0 = np.abs(X).sum()
            Y0 = np.abs(Y).sum()
            Z0 = np.abs(Z).sum()
            XYZ0 = [X0, Y0, Z0]

            # perpendicular vector [x,y,z]
            perpendicularVector = [False, False, False]
            # set
            if X0 == 0:
                perpendicularVector[0] = True

            if Y0 == 0:
                perpendicularVector[1] = True

            if Z0 == 0:
                perpendicularVector[2] = True

            # status
            if True in perpendicularVector:
                structureType = '2D'
            else:
                structureType = '3D'

            # axis selection
            perpendicularAxis = []
            if perpendicularVector[0] is True:
                perpendicularAxis.append('X')
            if perpendicularVector[1] is True:
                perpendicularAxis.append('Y')
            if perpendicularVector[2] is True:
                perpendicularAxis.append('Z')

            return structureType, perpendicularAxis, perpendicularVector, XYZ0
        except Exception as e:
            raise

    def set_color(self, atom_symbol):
        '''
        set color for compound representation

        Ref: https://en.wikipedia.org/wiki/CPK_coloring

        '''
        colors = {
            "H": 'ffffff',
            "C": 'BCBCBC',
            "N": '0586f6',
            "O": 'f6052a',
            "F": '2dd930',
            "Cl": '2dd930',
            "Br": '950e0e',
            "I": '360e89',
            "He": '3dbaf1',
            "Ne": '3dbaf1',
            "Ar": '3dbaf1',
            "Kr": '3dbaf1',
            "Xe": '3dbaf1',
            "P": 'f1a03d',
            "S": 'f1ef3d',
            "B": 'efc867',
            "Li": '6b3ccb',
            "Na": '6b3ccb',
            "K": '6b3ccb',
            "Rb": '6b3ccb',
            "Cs": '6b3ccb',
            "Fr": '6b3ccb',
            "Be": '1c881e',
            "Mg": '1c881e',
            "Ca": '1c881e',
            "Sr": '1c881e',
            "Ba": '1c881e',
            "Ra": '1c881e',
            "Ti": '3d3e40',
            "Fe": 'a48620',
            "other": 'a729ba'
        }

        _color = colors.get(str(atom_symbol))
        if _color is None:
            return '#'+colors.get(str('other'))
        else:
            return '#'+_color

    def set_size(self, symbol, _s=4):
        '''
        set size of atom spherical shape
        '''
        _sx = 100
        _sy = 200

        sizes = {
            "H": _s*_sx,
            "C": _s*_sy,
            "N": _s*_sy,
            "O": _s*_sy,
            "F": _s*_sy,
            "Cl": _s*_sx,
            "Br": _s*_sy,
            "I": _s*_sy,
            "He": _s*_sy,
            "Ne": _s*_sy,
            "Ar": _s*_sy,
            "Kr": _s*_sy,
            "Xe": _s*_sy,
            "P": _s*_sy,
            "S": _s*_sy,
            "B": _s*_sy,
            "Li": _s*_sy,
            "Na": _s*_sy,
            "K": _s*_sy,
            "Rb": _s*_sy,
            "Cs": _s*_sy,
            "Fr": _s*_sy,
            "Be": _s*_sy,
            "Mg": _s*_sy,
            "Ca": _s*_sy,
            "Sr": _s*_sy,
            "Ba": _s*_sy,
            "Ra": _s*_sy,
            "Ti": _s*_sy,
            "Fe": _s*_sy,
            "other": _s*_sy,
        }

        _size = sizes.get(symbol)
        if _size is None:
            _sizeSet = int(sizes.get('other'))
        else:
            _sizeSet = int(_size)

        return _sizeSet

    def create_3dframe(self):
        '''
        create 3d frame dimension
        '''
        # max length
        # x
        xMin = np.min(self.xyzList[:, 0])
        xMax = np.max(self.xyzList[:, 0])
        xLen = np.abs(np.abs(xMax) - np.abs(xMin))
        # y
        yMin = np.min(self.xyzList[:, 1])
        yMax = np.max(self.xyzList[:, 1])
        yLen = np.abs(np.abs(yMax) - np.abs(yMin))
        # z
        zMin = np.min(self.xyzList[:, 2])
        zMax = np.max(self.xyzList[:, 2])
        zLen = np.abs(np.abs(zMax) - np.abs(zMin))
        # max
        xyzLenMax = np.max([xMax, yMax, zMax])
        xyzLenMin = np.min([xMin, yMin, zMin])
        xyzR = 1/xyzLenMax

        # res
        return xyzLenMax, xyzLenMin, xyzR, xLen, yLen, zLen

    def create_bond_line(self, xyz1, xyz2, bond_type, xyzL=[1, 1, 1], xyzR=0.15, perpendicularAxis=[]):
        '''
        create bond line (single, double, tipple)

        args:
            xyz1: (x,y,z) point 1
            xyz2: (x,y,z) point 2
            bond_type: bond type (1,2,3)
        '''

        bondLines = []

        # check bond_type = 2,3,
        # if bond_type == 2 or bond_type == 3:
        #     xL, yL, zL = xyzR*xyzL
        #     # check geometry (perpendicular plane)
        #     # check
        #     if len(perpendicularAxis) > 0:
        #         # *** fix/random
        #         # random.choice(perpendicularAxis)
        #         randomAxis = perpendicularAxis[0]
        #         if randomAxis == 'X':
        #             yL = 0
        #             zL = 0
        #         elif randomAxis == 'Y':
        #             xL = 0
        #             zL = 0
        #         elif randomAxis == 'Z':
        #             xL = 0
        #             yL = 0
        #     else:
        #         raise Exception('perpendicularAxis error.')

        xL, yL, zL = xyzR*np.array(xyzL)

        if bond_type == 1:
            # parallel line
            # y increase
            # line
            _l1 = [[xyz1[0], xyz2[0]], [
                xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]]]
            bondLines.append(_l1)

        elif bond_type == 2:
            # parallel line
            # y increase
            # line
            _l1 = [[xyz1[0]+xL, xyz2[0]+xL], [
                xyz1[1]+yL, xyz2[1]+yL], [xyz1[2]+zL, xyz2[2]+zL]]
            _l2 = [[xyz1[0]-xL, xyz2[0]-xL], [
                xyz1[1]-yL, xyz2[1]-yL], [xyz1[2]-zL, xyz2[2]-zL]]
            bondLines.append(_l1)
            bondLines.append(_l2)

        elif bond_type == 3:
            # parallel line
            # y increase
            # line
            _l1 = [[xyz1[0]+xL, xyz2[0]+xL], [
                xyz1[1]+yL, xyz2[1]+yL], [xyz1[2]+zL, xyz2[2]+zL]]
            _l2 = [[xyz1[0], xyz2[0]], [xyz1[1], xyz2[1]], [xyz1[2], xyz2[2]]]
            _l3 = [[xyz1[0]-xL, xyz2[0]-xL], [
                xyz1[1]-yL, xyz2[1]-yL], [xyz1[2]-zL, xyz2[2]-zL]]
            bondLines.append(_l1)
            bondLines.append(_l2)
            bondLines.append(_l3)

        return bondLines, bond_type

    def line_property(self, xyz1, xyz2):
        '''
        check line property
            which plane are engaged

        hint:
            when res contains two True, it means the False coordination contains all elements. 
        '''
        try:
            # mean value
            Xm = np.mean([xyz1[0], xyz2[0]])
            Ym = np.mean([xyz1[1], xyz2[1]])
            Zm = np.mean([xyz1[2], xyz2[2]])
            xyzMean = [Xm, Ym, Zm]

            # check plane
            isSubtractZero = np.array([False, False, False])

            # points in one line
            X = xyz1[0] - xyz2[0]
            Y = xyz1[1] - xyz2[1]
            Z = xyz1[2] - xyz2[2]
            xyzPlane = [X, Y, Z]

            # axis vector
            xyzL = np.array([0, 0, 0])

            # set
            # axis selection
            perpendicularAxis = np.array([False, False, False])
            if X == 0:
                isSubtractZero[0] = True
                perpendicularAxis[0] = True
                # xyzL = xyzL + np.array(1, 1, 0)

            if Y == 0:
                isSubtractZero[1] = True
                perpendicularAxis[1] = True
                # xyzL = xyzL + np.array(1, 1, 0)

            if Z == 0:
                isSubtractZero[2] = True
                perpendicularAxis[2] = True
                # xyzL = xyzL + np.array(1, 0, 1)

            # if X != 0 and Y != 0 and Z != 0:
            #     xyzL = np.array(1, 1, 1)

            # set xyzL
            xyzL = isSubtractZero.astype(int)

            return xyzMean, xyzPlane, isSubtractZero, xyzL, perpendicularAxis

        except Exception as e:
            raise Exception(e)

    def view3d(self, elev=None, azim=None, figSize='default', obsOption=[False, 0]):
        '''
        draw compound in the cartesian coordinate

        args:
            atomElements: atom symbol
            atomBonds: atom bonds (bond blocks)
            xyzList: atom position in the cartesian coordinate
            elev=None: plt 3d setting
            azim=None: plt 3d setting
            figSize=(10, 10): plt 3d setting
            obsOption=[False, 0]: display center point [0,0,0]
        '''
        # 3d plot
        if figSize == 'default':
            fig = plt.figure(facecolor='#000000')
        else:
            fig = plt.figure(figsize=figSize, facecolor='#000000')

        # if figSize == 'default':
        #     fig = plt.figure()
        # else:
        #     fig = plt.figure(figsize=figSize)

        ax = plt.axes(projection='3d')
        # axis display
        plt.axis('off')
        # color
        ax.set_facecolor('#000000')

        # atom no
        atomNo = len(self.xyzList)
        # bond no
        bondNo = len(self.atomBonds)

        # create 3d frame
        xyzLenMax, xyzLenMin, xyzR, xLen, yLen, zLen = self.create_3dframe()
        # set
        # scale
        # ax.set_box_aspect((np.ptp([xyzLenMin, xyzLenMax]), np.ptp(
        #     [xyzLenMin, xyzLenMax]), np.ptp([xyzLenMin, xyzLenMax])))

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
            # size
            _atomSize = self.set_size(_atomSymbol)
            # color
            _atomColor = self.set_color(_atomSymbol)

            # atom mark
            atomMark = str(_atomSymbol) + str(_atomId)
            # draw atom 1
            ax.scatter3D(_atom1X, _atom1Y, _atom1Z,
                         label=_atomSymbol, s=_atomSize, c=_atomColor)
            # ax.text(_atom1X, _atom1Y, _atom1Z, '%s' %
            #         (atomMark), size=20, zorder=1, color='r')

        # reset
        i = 0

        # *** bond visualization
        # *** using bond block
        for i in range(bondNo):
            # atom id
            _atom1Id = int(self.atomBonds[i]['id']) - 1
            # atom symbol
            _atom1Symbol = self.atomBonds[i]['symbol']
            # atom bond list
            _atom1BondList = self.atomBonds[i]['bonds']
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

                    # set color
                    lineColor = ['b', 'r', 'g']
                    lineWidth = [1, 3, 5]

                    # xyz
                    _atom2X = self.xyzList[_atom2Id, 0]
                    _atom2Y = self.xyzList[_atom2Id, 1]
                    _atom2Z = self.xyzList[_atom2Id, 2]
                    _atom2XYZ = [_atom2X, _atom2Y, _atom2Z]

                    # line property
                    xyzMean, xyzPlane, isPlane, xyzL, perpendicularAxis = self.line_property(
                        _atom1XYZ, _atom2XYZ)

                    # bond connection (points)
                    _bondConnection, _bondTypeLog = self.create_bond_line(
                        _atom1XYZ, _atom2XYZ, _bondType)

                    # size
                    _bondConnectionSize = len(_bondConnection)

                    # check
                    if _bondConnectionSize == 1:
                        _vector = _bondConnection[0]
                        ax.plot3D(_vector[0], _vector[1], _vector[2],
                                  linewidth=lineWidth[_bondType-1], c=lineColor[_bondType-1])
                    else:
                        # line color: black
                        for b in range(_bondConnectionSize):
                            _vector = _bondConnection[b]
                            ax.plot3D(_vector[0], _vector[1], _vector[2],
                                      linewidth=lineWidth[_bondType-1], c=lineColor[_bondType-1])

            # obs show
            if obsOption[0]:
                ax.scatter3D(obsOption[1], 0, 0, s=40)

        # axis setting
        ax.set_xlabel("$X$")
        ax.set_ylabel("$Y$")
        ax.set_zlabel("$Z$")
        # ax.legend()
        # set limits
        _maxVal = np.max(self.xyzList)
        ax.set_xlim(float(-_maxVal), float(_maxVal))
        ax.set_ylim(float(-_maxVal), float(_maxVal))
        ax.set_zlim(float(-_maxVal), float(_maxVal))

        # set angles/elevations
        ax.view_init(elev=elev, azim=azim)
        plt.show()

    def view3dobs(self, elev=None, azim=None, figSize=(10, 10), obsOption=[True, 0]):
        '''
        draw compound in the cartesian coordinate with observer
        '''
        # 3d plot
        fig = plt.figure(figsize=figSize)
        ax = plt.axes(projection='3d')

        # atom no
        atomNo = len(self.xyzList)
        # bond no
        bondNo = len(self.atomBonds)

        # atom visualization
        for i in range(atomNo):
            # xyz
            _atomX = self.xyzList[i, 0]
            _atomY = self.xyzList[i, 1]
            _atomZ = self.xyzList[i, 2]
            # color
            # size

            # draw atom 1
            ax.scatter3D(_atomX, _atomY, _atomZ, s=40)

        # bond visualization
        for i in range(bondNo):
            # atom id
            _atom1Id = int(self.atomBonds[i]['id']) - 1
            # atom symbol
            _atom1Symbol = self.atomBonds[i]['symbol']
            # atom bond list
            _atom1BondList = self.atomBonds[i]['bonds']
            atom1BondSize = len(_atom1BondList)

            _atom1X = self.xyzList[:, 0]
            _atom1Y = self.xyzList[:, 1]
            _atom1Z = self.xyzList[:, 2]

            # draw bond
            if atom1BondSize > 0:
                for j in range(atom1BondSize):
                    # atom [2] id
                    _atom2Id = int(_atom1BondList[j][0]) - 1
                    # atom [2] symbol
                    _atom2Symbol = _atom1BondList[j][1]
                    # atom [1] - atom [2] bond type
                    _bondType = int(_atom1BondList[j][3])

                    # set color
                    lineColor = ['k', 'b', 'c']
                    lineWidth = [1, 2, 3]

                    # xyz
                    _atom2X = self.xyzList[_atom2Id, 0]
                    _atom2Y = self.xyzList[_atom2Id, 1]
                    _atom2Z = self.xyzList[_atom2Id, 2]

                    # bond connection
                    _bondConnection = self.xyzList[[_atom1Id, _atom2Id]]
                    # draw line
                    ax.plot3D(_bondConnection[:, 0], _bondConnection[:, 1], _bondConnection[:, 2],
                              color=lineColor[_bondType-1], linewidth=lineWidth[_bondType-1])

        # obs visualization
        xyzObsList = Observer.GeneratorCircleObserver(
            self._robs, self.tetaNo, self.phiNo, self.limits['teta'])[0]
        # obs size
        xyzObsSize = len(xyzObsList)
        for i in range(xyzObsSize):
            ax.scatter3D(xyzObsList[i, :, 0],
                         xyzObsList[i, :, 1], xyzObsList[i, :, 2])
            ax.plot3D(xyzObsList[i, :, 0],
                      xyzObsList[i, :, 1], xyzObsList[i, :, 2])

        # obs show
        if obsOption[0]:
            ax.scatter3D(obsOption[1], 0, 0, s=40)

        ax.view_init(elev=elev, azim=azim)
        plt.show()

    def create_line(self, xyzList1, xyzList2, t=1):
        '''
        create a line equation and its parallel lines
        '''
        try:
            r = np.array(xyzList2) - np.array(xyzList1)

            # matrix
            rMat = np.array([xyzList1, xyzList1])

            # x
            x = r[0]*t + xyzList1[0]
            # y
            y = r[1]*t + xyzList1[1]
            # z
            z = r[2]*t + xyzList1[2]

            xL, yL, zL = 0.1*np.ones(3)
            _l1 = [[xyzList1[0]+xL, xyzList2[0]+xL], [xyzList1[1] +
                                                      yL, xyzList2[1]+yL], [xyzList1[2]+zL, xyzList2[2]+zL]]
            _l2 = [[xyzList1[0]-xL, xyzList2[0]-xL], [xyzList1[1] -
                                                      yL, xyzList2[1]-yL], [xyzList1[2]-zL, xyzList2[2]-zL]]

            # res
            return r, _l1, _l2
        except Exception as e:
            raise Exception(e)

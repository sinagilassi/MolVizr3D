# CHEMICAL STRUCTURE
# --------------------

# import libs
import numpy as np
import matplotlib.pyplot as plt


class Structure():
    '''
    chemical structure methods
    '''

    def __init__(self):
        pass

    @staticmethod
    def CenterPoints(xyzList):
        '''
        find the center coordination of an object
        '''
        # set
        xyzList = np.array(xyzList)
        # find the highest xyz
        # x
        xMax = np.max(xyzList[:, 0])
        xMin = np.min(xyzList[:, 0])
        if np.abs(xMax) != np.abs(xMin):
            xLen = np.abs(xMax - xMin)
            xCenter = xMin + (xLen/2)
        else:
            xLen = 0
            xCenter = 0 + (xLen/2)
        # print(f"X value: {xMin}, {xMax}")

        # y
        yMax = np.max(xyzList[:, 1])
        yMin = np.min(xyzList[:, 1])
        if np.abs(yMax) != np.abs(yMin):
            yLen = np.abs(yMax - yMin)
            yCenter = yMin + (yLen/2)
        else:
            yLen = 0
            yCenter = 0 + (yLen/2)
        # print(f"Y value: {yMin}, {yMax}")

        # z
        zMax = np.max(xyzList[:, 2])
        zMin = np.min(xyzList[:, 2])
        if np.abs(zMax) != np.abs(zMin):
            zLen = np.abs(zMax - zMin)
            zCenter = zMin + (zLen/2)
        else:
            zLen = 0
            zCenter = 0 + (zLen/2)
        # print(f"Z value: {zMin}, {zMax}")

        # object base
        objectBaseCoordinate = np.array([xCenter, yCenter, zCenter])

        return objectBaseCoordinate

    @staticmethod
    def CenterObject(xyzList, centerPoint):
        '''
        move an object to the center of the origin [0,0,0]
        '''
        originPoint = np.array([0, 0, 0])
        movingCoordinate = originPoint - centerPoint
        newCenterPoints = np.array(xyzList) + movingCoordinate

        return newCenterPoints, movingCoordinate

    @staticmethod
    def ObjectDislay(xyzList, xyzCenterList):
        '''
        display object in two places
        '''
        # find the center of object
        xyzCenter = Structure.CenterPoints(xyzList)
        # find the center of object after moving to the origin [0,0,0]
        xyzCenterMoving = Structure.CenterPoints(xyzCenterList)

        # 3d plot
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter3D(xyzList[:, 0], xyzList[:, 1], xyzList[:, 2])
        ax.scatter3D(xyzCenter[0], xyzCenter[1], xyzCenter[2])
        ax.scatter3D(xyzCenterList[:, 0],
                     xyzCenterList[:, 1], xyzCenterList[:, 2])
        ax.scatter3D(xyzCenterMoving[0],
                     xyzCenterMoving[1], xyzCenterMoving[2])
        # line
        ax.plot3D([xyzCenter[0], xyzCenterMoving[0]], [xyzCenter[1],
                  xyzCenterMoving[1]], [xyzCenter[2], xyzCenterMoving[2]])
        plt.show()

    @staticmethod
    def CircleCoordinate(teta, r):
        '''
        circle coordination with teta and r

        args:
        teta: angle with x-axis
        r: circle radius
        '''
        x = r*np.cos(teta)
        y = r*np.sin(teta)
        return (x, y)

    @staticmethod
    def PeriodGenerator(n=100, w=1):
        '''
        set array of 2pi period
        '''
        recordsNo = w*n
        obsAnglesRes, obsAnglesStep = np.linspace(
            0, 2, recordsNo, retstep=True)
        obsAngles = obsAnglesRes*w*np.pi

        return obsAngles, obsAnglesStep

    @staticmethod
    def PeriodLimitGenerator(n=100, w=1, limits=[0, 0]):
        '''
        set array of 2pi period
        args:
            limits: list of min, max to be removed from the period
        '''
        # print(f"limits: {limits}")
        point1 = limits[0]
        point2 = limits[1]
        recordsNo = w*n
        obsAnglesRes, obsAnglesStep = np.linspace(
            point1, point2, recordsNo, retstep=True)
        obsAngles = obsAnglesRes

        return obsAngles, obsAnglesStep

    @staticmethod
    def SphericalToCartesianCoordinate(rtpPoint):
        '''
        convert spherical to cartesian points

        args:
            rtpPoint: r, teta, phi of spherical coordinate
        '''
        # spherical
        r = rtpPoint[0]
        teta = rtpPoint[1]
        phi = rtpPoint[2]
        # cartesian
        x = r*np.sin(phi)*np.cos(teta)
        y = r*np.sin(phi)*np.sin(teta)
        z = r*np.cos(phi)
        # res
        return np.array([x, y, z])

    @staticmethod
    def LineGenerator(n=100, w=1):
        recordsNo = w*n
        obsLengthRes, obsLengthStep = np.linspace(
            -10, 10, recordsNo, retstep=True)
        obsLength = obsLengthRes*w

        return obsLength, obsLengthStep

    @staticmethod
    def CartesianToSphericalCoordinate(xyzPoint):
        '''
        convert cartesian to spherical coordinate

        args:
            xyzPoint: x, y, z of cartesian coordinate
        '''
        # cartesian
        x = xyzPoint[0]
        y = xyzPoint[1]
        z = xyzPoint[2]
        # spherical
        r = np.sqrt(x**2 + y**2 + z**2)
        teta = np.arctan(y/x)
        # np.arccos(z/np.sqrt(x**2 + y**2 + z**2))
        phi = np.arctan(np.sqrt(x**2 + y**2)/z)
        # res
        return np.array([r, teta, phi]), np.array([r, np.degrees(teta), np.rad2deg(phi)])

    @staticmethod
    def create_formula(atom_elements):
        '''
        create mat using mat elements
        '''
        try:
            # check
            if len(atom_elements) == 0:
                raise Exception('atom elements list is empty')

            my_dict = {i: atom_elements.count(i) for i in atom_elements}
            # mat name
            elList = ''
            for key, value in my_dict.items():
                if value == 1:
                    _el = str(key)
                else:
                    _el = str(key)+str(value)
                elList += _el
                elList += ''
            # transform
            elList = elList.strip()

            return elList
        except Exception as e:
            raise

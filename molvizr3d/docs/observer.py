# MATERIAL OBSERVER
# -------------------

# import libs
import numpy as np
from .structure import Structure


class Observer():

    def __init__(self):
        pass

    @staticmethod
    def GenerateCircularObserver(xyzList, obsRadius, dataNo):
        '''
        generate observer points based on circle's track
        '''
        # angular frequency
        freq = 3
        # obs angles
        obsAnglesRes = Structure.PeriodGenerator(dataNo)
        obsAngles = obsAnglesRes[0]
        obsAnglesStep = obsAnglesRes[1]

        # circle loops
        obsAnglesLoopRes = Structure.PeriodGenerator(dataNo, freq)
        obsAnglesLoop = obsAnglesLoopRes[0]
        obsAnglesLoopStep = obsAnglesLoopRes[1]

        # len
        obsAnglesLength = len(obsAngles)
        # observer coordinate
        obsCoordinate = np.zeros((freq, obsAnglesLength, 3))

        # circle 1
        for i in range(obsAnglesLength):
            _x, _y = Structure.CircleCoordinate(obsAngles[i], obsRadius)
            # outside compound
            _obsCoordinate = [_x, _y, 0]
            # inside compound
            # save
            obsCoordinate[0, i, :] = _obsCoordinate
            # reset
            _x, _y = 0, 0

        # circle 2
        for i in range(obsAnglesLength):
            _y, _z = Structure.CircleCoordinate(obsAngles[i], obsRadius)
            # outside compound
            _obsCoordinate = [0, _y, _z]
            # inside compound
            # save
            obsCoordinate[1, i, :] = _obsCoordinate
            # reset
            _y, _z = 0, 0

        # circle 3
        for i in range(obsAnglesLength):
            _x, _z = Structure.CircleCoordinate(obsAngles[i], obsRadius)
            # outside compound
            _obsCoordinate = [_x, 0, _z]
            # inside compound
            # save
            obsCoordinate[2, i, :] = _obsCoordinate
            # reset
            _x, _z = 0, 0

        # res
        return obsCoordinate, obsAngles, obsAnglesLoop, obsAnglesLoopStep

    @staticmethod
    def GeneratorLinearObserver(xyzList, obsDistance, dataNo):
        # number of observers
        obsNo = 3
        # obs length
        obsLengthRes = Structure.LineGenerator(dataNo)
        obsLength = obsLengthRes[0]
        obsLengthStep = obsLengthRes[1]

        # obs total length
        obsLengthLoopRes = Structure.LineGenerator(dataNo, w=obsNo)
        obsLengthLoop = obsLengthLoopRes[0]
        obsLengthLoopStep = obsLengthLoopRes[1]

        # len
        obsLengthSize = len(obsLength)
        # observer cordinate
        obsCoordinate = np.zeros((obsNo, obsLengthSize, 3))

        # line 1 (parallel x)
        for i in range(obsLengthSize):
            _x, _y, _z = obsLength[i], obsDistance, 0
            # outside compound
            _obsCoordinate = [_x, _y, _z]
            # inside compound
            # save
            obsCoordinate[0, i, :] = _obsCoordinate
            # reset
            _x, _y, _z = 0, 0, 0

        # line 2 (parallel y)
        for i in range(obsLengthSize):
            _x, _y, _z = obsDistance, obsLength[i], 0
            # outside compound
            _obsCoordinate = [_x, _y, _z]
            # inside compound
            # save
            obsCoordinate[1, i, :] = _obsCoordinate
            # reset
            _x, _y, _z = 0, 0, 0

        # line 1 (parallel z)
        for i in range(obsLengthSize):
            _x, _y, _z = obsLength[i], 0, obsDistance
            # outside compound
            _obsCoordinate = [_x, _y, _z]
            # inside compound
            # save
            obsCoordinate[2, i, :] = _obsCoordinate
            # reset
            _x, _y, _z = 0, 0, 0

        # res
        return obsCoordinate, obsLength, obsLengthLoop, obsLengthLoopStep

    @staticmethod
    def GeneratorCircleObserver(r, tetaNo, phiNo, limits=[]):
        '''
        Generate circle xyz points in cartesian coordinate
        '''
        # rad [rad]
        tetaRes = Structure.PeriodGenerator(n=tetaNo)
        teta = tetaRes[0]
        tetaStep = tetaRes[1]
        # print(f"teta: {teta.shape}")
        tetaNo = len(teta)
        # phi [rad]
        if len(limits) > 0:
            phiRes = Structure.PeriodLimitGenerator(n=phiNo, limits=limits)
            phi = phiRes[0]
            phiStep = phiRes[1]
        else:
            phiRes = Structure.PeriodGenerator(n=phiNo)
            phi = phiRes[0]
            phiStep = phiRes[1]
        # print(f"phi: {phi.shape}")
        phiNo = len(phi)

        # total rotation
        angFreq = Structure.PeriodGenerator(n=tetaNo, w=tetaNo)[0]
        # print(f"angFreq: {angFreq}")
        # frequency
        freq = np.max(angFreq)/(2*np.pi)
        # print(f"freq: {freq}")
        # period
        period = 1/freq
        # print(f"period: {period}")
        # sample no [in a second]
        sampleNo = len(angFreq)
        # print(f"sampleNo: {sampleNo}")
        # sampling rate [number of samples in a second]
        samplingRate = sampleNo
        # sampling interval
        samplingInterval = 1/samplingRate
        # time span
        timeSpan = np.arange(0, 1, samplingInterval)

        # [number of circles, number of points in a circle, [x,y,z] points]
        xyzPoints = np.zeros((tetaNo, phiNo, 3))

        for i in range(tetaNo):
            for j in range(phiNo):
                _rtpPoint = np.array([r, teta[i], phi[j]])
                _xyzPoint = Structure.SphericalToCartesianCoordinate(_rtpPoint)
                # save
                xyzPoints[i, j, :] = _xyzPoint

        # print(f"xyzPoints: {xyzPoints.shape}")
        # res
        return xyzPoints, teta, phi, tetaStep, phiStep, angFreq, freq, period, sampleNo, samplingRate, samplingInterval, timeSpan

    @staticmethod
    def ObsWatchPathGenerator(r, tetaNo, phiNo, tetaLimits=[], phiLimit=[]):
        '''
        generate a circle path containing xyz points in the cartesian coordinate

        args:
            r: distance between observer points and element points
            tetaNo: number of circle paths
            phiNo: number of observer points in a circle path
            tetaLimits=[]: angles define the limit of observer, default: [0, pi]
            phiLimit=[]: angles define the limit of observer, default: [0, 2pi]

        hints:
            the last obs point is the mirror of angle 0, thus it is ignored.
        '''
        # rad [rad]
        if len(tetaLimits) > 0:
            tetaRes = Structure.PeriodLimitGenerator(
                n=tetaNo+1, limits=tetaLimits)
            teta = tetaRes[0]
            tetaStep = tetaRes[1]
        else:
            tetaRes = Structure.PeriodGenerator(n=tetaNo+1)
            teta = tetaRes[0]
            tetaStep = tetaRes[1]

        # phi [rad]
        if len(phiLimit) > 0:
            phiRes = Structure.PeriodLimitGenerator(
                n=phiNo+1, limits=phiLimit)
            phi = phiRes[0]
            phiStep = phiRes[1]
        else:
            phiRes = Structure.PeriodGenerator(n=phiNo+1)
            phi = phiRes[0]
            phiStep = phiRes[1]

        # time span [for one loop]
        timeLoop = 1
        # sampling no [1 second]
        sampleNo = phiNo
        # sampling rate [number of samples in a second]
        samplingRate = sampleNo
        # sampling interval
        samplingInterval = 1/samplingRate
        # time span
        timeSpan = np.arange(0, timeLoop, samplingInterval)
        # frequency set
        freq = 1
        # total rotation [1 second]
        angFreq = 2*np.pi*freq
        # period
        period = 1/freq

        # total sampling time [s]
        samplingTimeTotal = tetaNo
        # total number of samples (total time)
        sampleNoTotal = samplingTimeTotal*sampleNo
        # sampling rate with respect to total time
        samplingIntervalTotal = samplingTimeTotal/sampleNoTotal
        # time span
        timeSpanTotal = np.arange(0, samplingTimeTotal, samplingIntervalTotal)

        # [number of circles, number of points in a circle, [x,y,z] points]
        xyzPoints = np.zeros((tetaNo, phiNo, 3))

        for i in range(tetaNo):
            for j in range(phiNo):
                _rtpPoint = np.array([r, teta[i], phi[j]])
                _xyzPoint = Structure.SphericalToCartesianCoordinate(_rtpPoint)
                # save
                xyzPoints[i, j, :] = _xyzPoint

        # res
        return xyzPoints, teta, phi, tetaStep, phiStep, angFreq, freq, period, sampleNo, samplingRate, samplingInterval, timeSpan, timeSpanTotal

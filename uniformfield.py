# x full width 32cm
# y full width 77.6mm
# z full width 1m

import numpy as _np
import argparse
import pybdsim


def generate(magnetic_field_strength, length):

    B = magnetic_field_strength #T
    use_length = length*100     #convert metres to centimetres

    # in cm
    xFullWidth = 32
    yFullWidth = 7.76
    zFullWidth = 100

    # pole position offsets in cm
    poleposx = 13
    poleposy = 0
    poleposz = -0.5*use_length + 150.5  # converts from length of spectrometer
                                        # to z position of pole

    fullWidths = _np.array([xFullWidth, yFullWidth, zFullWidth])
    positivePoint = 0.25 * fullWidths

    x = [poleposx - positivePoint[0], poleposx + positivePoint[0]]
    y = [poleposy - positivePoint[1], poleposy + positivePoint[1]]
    z = [poleposz - positivePoint[2], poleposz + positivePoint[2]]

    data = []
    # loop over and build up 3d lists of lists of lists
    for xi in x:
        u = []
        for yi in y:
            v = []
            for zi in z:
                v.append([xi,yi,zi,0,B,0])
            u.append(v)
        data.append(u)

    # convert to numpy array
    data = _np.array(data)

    # construct a BDSIM format field object and write it out
    f = pybdsim.Field.Field3D(data)
    outputfilename = '/unix/pdpwa/jchappell/data/fields/BDSIM_uniform_field__' \
                     + str(B) + 'T__length_' + str(length) + 'm.dat'
    f.Write(outputfilename)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""
    This script generates a uniform magnetic field corresponding to the 
    dipole magnet within the BDSIM Awake Spectrometer class.""",
                                     formatter_class=
                                     argparse.RawTextHelpFormatter)

    parser.add_argument('--fieldstrength', dest='fieldStrength',
                        default=None, help = """
                        This is the value of the desired uniform field 
                        strength measured in Tesla.""")

    parser.add_argument('--length', dest='length', default=None, help = """
    This is the value of the length supplied to the gmad file that 
    initialises the Awake spectrometer class. Measured in metres.""")

    arguments = parser.parse_args()

    _fieldStrength = float(arguments.fieldStrength)

    _length = float(arguments.length)

    generate(_fieldStrength, _length)

    print "Field generated with strength: " + str(_fieldStrength) + 'T.'
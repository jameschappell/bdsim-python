import pybdsim
import matplotlib.pyplot as plt
import numpy as np


class FourDData(object):
    def __init__(self, filename, xind=0, yind=1, zind=2, tind=3):
        d = pybdsim.Field.Load(filename)

        # '...' fills in unknown number of dimensions with ':' meaning
        # all of that dimension
        if (xind >= 0):
            self.x = d[..., xind].flatten()
        if (yind >= 0):
            self.y = d[..., yind].flatten()
        if (zind >= 0):
            self.z = d[..., zind].flatten()
        if (tind >= 0):
            self.t = d[..., tind].flatten()

        # index from end as we don't know the dimensionality
        self.fx = d[..., -3].flatten()
        self.fy = d[..., -2].flatten()
        self.fz = d[..., -1].flatten()

        self.mag = _np.sqrt(self.fx ** 2 + self.fy ** 2 + self.fz ** 2)


class ThreeDData(FourDData):
    def __init__(self, filename):
        FourDData.__init__(self, filename, tind=-1)


def Niceties(xlabel, ylabel):
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar()
    plt.tight_layout()
    plt.axes().set_aspect('equal', 'datalim')


def Plot3DXYZ(filename, scale=None):
    """
    Plots (B_x, B_y, B_z) as function of x, y and z.
    """

    d = ThreeDData(filename)
    plt.figure()
    plt.quiver(d.x, d.y, d.z, d.fx, d.fy, d.fz, cmap=plt.cm.magma,
               pivot='mid', scale=scale)
    Niceties('X (cm)', 'Y (cm)')


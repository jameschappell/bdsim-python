import pybdsim
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


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

        self.mag = np.sqrt(self.fx ** 2 + self.fy ** 2 + self.fz ** 2)


class ThreeDData(FourDData):
    def __init__(self, filename):
        FourDData.__init__(self, filename, tind=-1)


def Niceties(xlabel, ylabel, zlabel):
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #plt.zlabel(zlabel)
    #plt.colorbar()
    plt.tight_layout()
    plt.axes().set_aspect('equal', 'datalim')


def Plot3DXYZ(filename, scale=None):
    """
    Plots (B_x, B_y, B_z) as function of x, y and z.
    """

    points = np.array([[-3.0, -3.88, -99.5],
                       [-3.0, 3.88, -99.5],
                       [29.0, -3.88, -99.5],
                       [29.0, 3.88, -99.5],
                       [-3.0, -3.88, 0.5],
                       [-3.0, 3.88, 0.5],
                       [29.0, -3.88, 0.5],
                       [29.0, 3.88, 0.5]])

    d = ThreeDData(filename)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.quiver(d.x, d.y, d.z, d.fx, d.fy, d.mag, cmap=plt.cm.magma,
               pivot='middle', length=0.2)

    ax.scatter3D(points[:, 0], points[:, 1], points[:, 2], c='r')

    ax.set_xlim(-10, 40)
    ax.set_ylim(-5, 5)
    ax.set_zlim(-120, 20)
    ax.set_xlabel('X (cm)')
    ax.set_ylabel('Y (cm)')
    ax.set_zlabel('Z (cm)')
    #Niceties('X (cm)', 'Y (cm)', 'Z (cm)')


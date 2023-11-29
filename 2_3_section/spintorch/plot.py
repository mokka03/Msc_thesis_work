import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, CenteredNorm
from matplotlib.ticker import MaxNLocator
from .geom import WaveGeometryMs, WaveGeometry
from .solver import MMSolver

def plot_geom(plotdir,geom,epoch=0):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    h1 = ax.imshow(geom,origin="lower", cmap=plt.cm.summer)
    plt.colorbar(h1, ax=ax, label='Saturation magnetization (A/m)')
    fig.savefig(plotdir+"geom_%d.png" %(epoch))
    plt.close(fig)


def plot_m(plotdir,m,epoch=0):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    h1 = ax.imshow(m,origin="lower", cmap=plt.cm.winter)
    plt.colorbar(h1, ax=ax, label='Magnetic field (a.u.)')
    fig.savefig(plotdir+"m_%d.png" %(epoch))
    plt.close(fig)


def plot_loss(loss_iter, plotdir):
    fig = plt.figure()
    plt.plot(loss_iter, 'o-')
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    fig.savefig(plotdir+'loss.png')
    plt.close(fig)
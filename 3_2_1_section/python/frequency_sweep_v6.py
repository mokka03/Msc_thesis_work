"""Optimize a focusing model"""
import torch
import os
import numpy as np
from scipy.io import savemat, loadmat

import spintorch
from spintorch.utils import tic, toc, stat_cuda
from spintorch.plot import wave_integrated, wave_snapshot

import warnings
warnings.filterwarnings("ignore", message=".*Casting complex values to real.*")

'''Directories'''
basedir = 'demux_v6_frequency_sweep/'
plotdir = 'plots/' + basedir
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)
savedir = 'models/' + basedir
if not os.path.isdir(savedir):
    os.makedirs(savedir)    

"""Parameters"""
dx = 100e-9      # discretization (m)
dy = 100e-9      # discretization (m)
dz = 99e-9      # discretization (m)
nx = 700        # size x    (cells)
ny = 740        # size y    (cells)
damping_with = 60
desired_output = 4

Ms0 = 116.23e3      # saturation magnetization (A/m)
Ms_max = 118.72e3      # saturation magnetization (A/m)
Ms_min = 116.23e3      # saturation magnetization (A/m)
B0 = 214.1e-3     # bias field (T)
Bt = 1e-3       # excitation field amplitude (T)

dt = 20e-12     # timestep (s)
f1 = 2.05e9        # source frequency 1 (Hz)
f2 = 2.1e9    # source frequency 2 (Hz)
timesteps = 7200 # number of timesteps for wave propagation

dev = torch.device('cuda')  # 'cuda' or 'cpu'

'''Geometry, sources, probes, model definitions'''
geom = spintorch.WaveGeometryMs((nx, ny), (dx, dy, dz), (Ms0,Ms_max,Ms_min), B0,damping_with)
src = spintorch.WaveLineSource(60, 0, 60, ny-1, dim=1)
probes = []
Np = 16  # number of probes
probe_radius = 15
for p in range(Np):
    probes.append(spintorch.WaveIntensityProbeDisk(nx-240, int(145+probe_radius*2*p), probe_radius))
model = spintorch.MMSolver(geom, dt, [src], probes)

print('Running on', dev)
model.to(dev)   # sending model to GPU/CPU

## Load trained Msat
Msat = loadmat('Msat_v6.mat').get('Msat')
Msat = torch.tensor(Msat.transpose())
model.geom.Msat0[:,:] = Msat

'''Define the source signal'''
## 2 input frequencies
t = torch.arange(0, timesteps*dt, dt, device=dev).unsqueeze(0).unsqueeze(2) # time vector
    
'''Plot spin-wave propagation'''
model.retain_history = False
with torch.no_grad():
    
    ## Save Msat
    spintorch.plot.geometry(model, epoch=0, plotdir=plotdir)
    Msat = model.geom.Msat.detach()
    savemat(savedir + "Msat.mat", {"Msat": Msat.to(torch.device("cpu")).numpy().transpose()})
    
    for f_MHz in range(1950,2205,5):
        f = f_MHz*1e6
        INPUTS = Bt*torch.sin(2*np.pi*f*t)
        
        ## Calculate wave propagation
        tic()
        u = model(INPUTS).sum(dim=1)
        stat_cuda("after wave propagation")
        toc()

        ## Save my
        my = model.m_last[0,1,]-model.m0[0,1,].unsqueeze(0).cpu()
        wave_snapshot(model, my, (plotdir+'snapshot_time%d_%dMHz.png' % (timesteps,f_MHz)),r"$m_z$")
        wave_integrated(model, my, (plotdir+'integrated_%dMHz.png' % (f_MHz)))
        savemat(savedir+"my_%dMHz.mat" % (f_MHz), {"my": my.to(torch.device("cpu")).numpy().transpose()})

        my = model.m_half[0,1,]-model.m0[0,1,].unsqueeze(0).cpu()
        wave_snapshot(model, my, (plotdir+'snapshot_time%d_%dMHz.png' % (timesteps/2,f_MHz)),r"$m_z$")

        ## Save input
        # savemat(savedir+"src.mat", {"src": INPUTS.to(torch.device("cpu")).numpy()})

        ## Save output
        savemat(savedir+"outputs_%dMHz.mat" % (f_MHz), {"outputs": u[0,].to(torch.device("cpu")).numpy()})

    # Plot and save alpha
    import matplotlib.pyplot as plt
    alpha = model.Alpha(False).to(torch.device("cpu")).numpy().transpose()
    alpha2D = np.flip(alpha[:,:,0,0],0)
    plt.imshow(alpha2D)
    plt.savefig(plotdir+"alpha.png")
    savemat(savedir+"alpha.mat", {"alpha": alpha})
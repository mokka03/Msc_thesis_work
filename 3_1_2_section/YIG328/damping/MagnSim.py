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
basedir = 'MagnSim/'
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
ny = 1000        # size y    (cells)
damping_with = 60
desired_output = 4

Ms0 = 116.23e3      # saturation magnetization (A/m)
Ms_max = 118.72e3      # saturation magnetization (A/m)
Ms_min = 116.23e3      # saturation magnetization (A/m)
B0 = 214.3e-3     # bias field (T)
Bt = 1e-3       # excitation field amplitude (T)

dt = 20e-12     # timestep (s)
f1 = 2.05e9        # source frequency 1 (Hz)
# f2 = 2.15e9    # source frequency 2 (Hz)
timesteps = 7000 # number of timesteps for wave propagation

dev = torch.device('cuda')  # 'cuda' or 'cpu'

for b in range(2141,2146,1):
    B0 = b*1e-4
    '''Geometry, sources, probes, model definitions'''
    geom = spintorch.WaveGeometryMs((nx, ny), (dx, dy, dz), (Ms0,Ms_max,Ms_min), B0,damping_with)
    src = spintorch.WaveLineSource(60, 0, 60, ny-1, dim=1)
    probes = []
    Np = 16  # number of probes
    probe_radius = 15
    for p in range(Np):
        probes.append(spintorch.WaveIntensityProbeDisk(nx-140, int(145+probe_radius*2*p), probe_radius))
    model = spintorch.MMSolver(geom, dt, [src], probes)

    print('Running on', dev)
    model.to(dev)   # sending model to GPU/CPU

    '''Define the source signal'''
    ## 2 input frequencies
    t = torch.arange(0, timesteps*dt, dt, device=dev).unsqueeze(0).unsqueeze(2) # time vector
    X1 = Bt*torch.sin(2*np.pi*f1*t)
    # X2 = Bt*torch.sin(2*np.pi*f2*t)


    epoch = 0
    '''Plot spin-wave propagation'''
    model.retain_history = False
    with torch.no_grad():        
        for i in range(1):
            if i == 0:
                INPUTS = X1
            # else:
            #     INPUTS = X2
            
            ## Calculate wave propagation
            tic()
            u = model(INPUTS).sum(dim=1)
            stat_cuda("after wave propagation")
            toc()

            ## Save my
            my = model.m_last[0,1,]-model.m0[0,1,].unsqueeze(0).cpu()
            savemat(savedir+"my_%d_%d.mat" % (b,i), {"my": my.to(torch.device("cpu")).numpy().transpose()})
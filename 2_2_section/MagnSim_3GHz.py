"""Optimize a focusing model"""
import torch
import os
import numpy as np
from scipy.io import savemat, loadmat

import spintorch
from spintorch.utils import tic, toc, stat_cuda
from spintorch.OOMMFio import OOMMF2torch

import warnings
warnings.filterwarnings("ignore", message=".*Casting complex values to real.*")

'''Directories'''
basedir = 'MagnSim_3GHz/'
plotdir = 'plots/' + basedir
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)
savedir = 'models/' + basedir
if not os.path.isdir(savedir):
    os.makedirs(savedir)    

"""Parameters"""
dx = 1000e-9      # discretization (m)
dy = 1000e-9      # discretization (m)
dz = 540e-9      # discretization (m)
nx = 500        # size x    (cells)
ny = 400        # size y    (cells)
damping_with = 50
desired_output = 4

Ms0 = 140e3      # saturation magnetization (A/m)
Ms_max = 140e3      # saturation magnetization (A/m)
Ms_min = 140.4e3      # saturation magnetization (A/m)
B0 = 278e-3     # bias field (T)
Bt = 1e-3       # excitation field amplitude (T)

dt = 10e-12     # timestep (s)
f1 = 3e9        # source frequency 
timesteps = 12000 # number of timesteps for wave propagation


'''Geometry, sources, probes, model definitions'''
### Geometry
geom = spintorch.WaveGeometryMs((nx, ny), (dx, dy, dz), (Ms0,Ms_max,Ms_min), B0,damping_with)
### Sources
source_field = OOMMF2torch('CPW_DL_Real_3GHz.ohf')
source_field = source_field.transpose(1,2)
source_field = source_field/torch.max(source_field[0,])
src = spintorch.CPWSource(source_field)
### Probes
probes = []
Np = 19  # number of probes
for p in range(Np):
    probes.append(spintorch.WaveIntensityProbeDisk(nx-15, int(ny*(p+1)/(Np+1)), 2))
model = spintorch.MMSolver(geom, dt, timesteps, [src], probes)

dev = torch.device('cuda')  # 'cuda' or 'cpu'
print('Running on', dev)
model.to(dev)   # sending model to GPU/CPU

'''Define the source signal and output goal'''
t = torch.arange(0, timesteps*dt, dt, device=dev).unsqueeze(0).unsqueeze(2) # time vector
X = Bt*torch.sin(2*np.pi*f1*t)  # sinusoid signal at f1 frequency, Bt amplitude

INPUTS = X  # here we could cat multiple inputs
OUTPUTS = torch.tensor([int(Np/2)]).to(dev) # desired output

'''Plot spin-wave propagation'''
epoch = 0
with torch.no_grad():
    model.retain_history = True
    ## Calculate wave propagation
    tic()
    u = model(INPUTS).sum(dim=1)
    stat_cuda("after wave propagation")
    toc()


    ## Calculate the induced voltage in CPW
    m = model.m_lasts

    # m = torch.stack(model.m_lasts, 1)[0,]-model.m0[0,].unsqueeze(0).cpu()
    # m = m.to(dev)
    nt = m.shape[0]
    cpw_xpos = torch.tensor([349,350,351,352,  357,358,359,360,  365,366,367,368]) # position of CPW in x dimension
    cpw_zpos = torch.tensor(420e-9) # position of CPW in z dimension
    outputCPW = spintorch.PickupCPW(cpw_xpos, cpw_zpos, (nx,ny,nt),(dx,dy,dz), dt, timesteps)
    outputCPW.to(dev)   # sending model to GPU/CPU
    print(m.shape)
    tic()
    avg = outputCPW(m, model.geom.Msat)
    toc()
    savemat(savedir + "avg.mat", {"avg": avg.cpu().numpy()})
    print(avg)

    ## Save Msat
    Msat = model.geom.Msat.detach().cpu().numpy().transpose()
    spintorch.plot.plot_geom(plotdir,Msat)
    savemat(savedir + "Msat.mat", {"Msat": Msat})

    ## Save my
    my = (model.m_last[0,1,]-model.m0[0,1,].unsqueeze(0).cpu()).numpy().transpose()
    spintorch.plot.plot_m(plotdir,my)
    savemat(savedir+"my.mat", {"my": my})

    # ## Save m
    # m = torch.stack(model.m_lasts, 1)[0,]-model.m0[0,].unsqueeze(0).cpu()
    # savemat(savedir+"m.mat", {"m": m.to(torch.device("cpu")).numpy().transpose()})
    # print(m.shape)

    # # Plot and save alpha
    # import matplotlib.pyplot as plt
    # alpha = model.Alpha(False).to(torch.device("cpu")).numpy().transpose()
    # alpha2D = np.flip(alpha[:,:,0,0],0)
    # plt.imshow(alpha2D)
    # plt.savefig(plotdir+"alpha.png")
    # savemat(savedir+"alpha.mat", {"alpha": alpha})

    # # Plot my
    # plt.imshow(my[0,].transpose(0,1))
    # plt.savefig(plotdir+"my.png")
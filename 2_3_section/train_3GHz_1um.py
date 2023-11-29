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
basedir = 'train_3GHz_1um/'
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
nx = 250        # size x    (cells)
ny = 200        # size y    (cells)
damping_with = 50
desired_output = 4

Ms0 = 140e3      # saturation magnetization (A/m)
Ms_max = 139.5e3      # saturation magnetization (A/m)
Ms_min = 141.5e3      # saturation magnetization (A/m)
B0 = 274e-3     # bias field (T)

dt = 20e-12     # timestep (s)
f1 = 3e9        # source frequency 
timesteps = 6000 # number of timesteps for wave propagation


'''Geometry, sources, probes, model definitions'''
### Geometry
geom = spintorch.WaveGeometryMs((nx, ny), (dx, dy, dz), (Ms0,Ms_max,Ms_min), B0,damping_with)
### Sources
source_field = loadmat('CPW_DL_Real_3GHz_1um.mat').get('excitation_mask')
source_field = torch.tensor(source_field).float()
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
X = torch.sin(2*np.pi*f1*t)  # sinusoid signal at f1 frequency, Bt amplitude

INPUTS = X  # here we could cat multiple inputs
OUTPUTS = torch.tensor([int(Np/2)]).to(dev) # desired output

'''Define optimizer and lossfunction'''
optimizer = torch.optim.Adam(model.parameters(), lr=0.5)

'''Load checkpoint'''
epoch = epoch_init = 19 # select previous checkpoint (-1 = don't use checkpoint)
if epoch_init>=0:
    checkpoint = torch.load(savedir + 'model_e%d.pt' % (epoch_init))
    epoch = checkpoint['epoch']
    loss_iter = checkpoint['loss_iter']
    model.load_state_dict(checkpoint['model_state_dict'])
    optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
else:
    loss_iter = []

nt = model.timesteps-model.tlasts
cpw_xpos = torch.tensor([175,176,177,178,  183,184,185,186,  191,192,193,194]) # position of CPW in x dimension
cpw_zpos = torch.tensor(420e-9) # position of CPW in z dimension

outputCPW = spintorch.PickupCPW(cpw_xpos, cpw_zpos, (nx,ny,nt),(dx,dy,dz), dt, timesteps)
outputCPW.to(dev)   # sending model to GPU/CPU

'''Train the network'''
tic()
model.retain_history = True
for epoch in range(epoch_init+1, 40):
    optimizer.zero_grad()
    u = model(INPUTS).sum(dim=1)

    avg = outputCPW(model.m_lasts, model.geom.Msat)
    print(avg)
    loss = -avg*1e7
    loss_iter.append(loss.item())  # store loss values
    spintorch.plot.plot_loss(loss_iter, plotdir)
    stat_cuda('after forward')
    loss.backward()
    optimizer.step()
    stat_cuda('after backward')
    print("Epoch finished: %d -- Loss: %.6f" % (epoch, loss))
    toc()   

    ## Save
    spintorch.plot.plot_geom(plotdir,model.geom.Msat.detach().cpu().numpy().transpose(),epoch=epoch)
    savemat(savedir+"Msat%d.mat" %(epoch), {"Msat": model.geom.Msat.detach().cpu().numpy().transpose()})
    spintorch.plot.plot_m(plotdir,(model.m_last[0,1,]-model.m0[0,1,].unsqueeze(0).cpu()).numpy().transpose(),epoch=epoch)
    savemat(savedir+"m%d.mat" %(epoch), {"m": (model.m_lasts).detach().cpu().numpy().transpose()})
    '''Save model checkpoint'''
    torch.save({
                'epoch': epoch,
                'loss_iter': loss_iter,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict()
                }, savedir + 'model_e%d.pt' % (epoch))




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
    # m = torch.stack(model.m_lasts, 1)[0,]-model.m0[0,].unsqueeze(0).cpu()
    # m = m.to(dev)
    # nt = m.shape[0]
    # outputCPW = spintorch.PickupCPW(torch.arange(349,361), 420e-9, (nx,ny,nt),(dx,dy,dz), dt)
    # outputCPW.to(dev)   # sending model to GPU/CPU
    # tic()
    # avg = outputCPW(m, model.geom.Msat)
    # toc()
    # savemat(savedir + "deltaV.mat", {"deltaV": deltaV.cpu().numpy()})

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
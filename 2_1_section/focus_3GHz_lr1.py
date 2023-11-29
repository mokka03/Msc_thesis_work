"""Optimize a focusing model"""
import torch
import os
import numpy as np
from scipy.io import savemat, loadmat
import matplotlib.pyplot as plt

import spintorch
from spintorch.utils import tic, toc, stat_cuda
from spintorch.OOMMFio import OOMMF2torch
from spintorch.plot import wave_integrated, wave_snapshot

import warnings
warnings.filterwarnings("ignore", message=".*Casting complex values to real.*")

'''Directories'''
basedir = 'focus_3GHz_lr1/'
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
Ms_max = 140.3e3      # saturation magnetization (A/m)
Ms_min = 139.7e3      # saturation magnetization (A/m)
B0 = 278e-3     # bias field (T)
Bt = 1e-3       # excitation field amplitude (T)

dt = 20e-12     # timestep (s)
f1 = 3e9        # source frequency 
timesteps = 6000 # number of timesteps for wave propagation

dev = torch.device('cuda')  # 'cuda' or 'cpu'

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
Np = 5  # number of probes
for p in range(Np):
    probes.append(spintorch.WaveIntensityProbeDisk(nx-100, int(ny*(p+1)/(Np+1)), 10))
model = spintorch.MMSolver(geom, dt, [src], probes)

dev = torch.device('cuda')  # 'cuda' or 'cpu'
print('Running on', dev)
model.to(dev)   # sending model to GPU/CPU

'''Define the source signal and output goal'''
t = torch.arange(0, timesteps*dt, dt, device=dev).unsqueeze(0).unsqueeze(2) # time vector
X = 0.5*Bt*torch.sin(2*np.pi*f1*t)  # sinusoid signal at f1 frequency, Bt amplitude

INPUTS = X  # here we could cat multiple inputs
OUTPUTS = torch.tensor([int(Np/2)]).to(dev) # desired output

'''Define optimizer and lossfunction'''
optimizer = torch.optim.Adam(model.parameters(), lr=1)

def my_loss(output, target_index):
    target_value = output[:,target_index]
    loss = output.sum(dim=1)/target_value-1
    return (loss.sum()/loss.size()[0]).log10()

'''Load checkpoint'''
epoch = epoch_init = -1 # select previous checkpoint (-1 = don't use checkpoint)
if epoch_init>=0:
    checkpoint = torch.load(savedir + 'model_e%d.pt' % (epoch_init))
    epoch = checkpoint['epoch']
    loss_iter = checkpoint['loss_iter']
    model.load_state_dict(checkpoint['model_state_dict'])
    optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
else:
    loss_iter = []


'''Train the network'''
tic()
model.retain_history = True
for epoch in range(epoch_init+1, 10):
    optimizer.zero_grad()
    u = model(INPUTS).sum(dim=1)
    spintorch.plot.plot_output(u[0,], OUTPUTS[0]+1, epoch, plotdir)
    loss = my_loss(u,OUTPUTS)
    loss_iter.append(loss.item())  # store loss values
    spintorch.plot.plot_loss(loss_iter, plotdir)
    stat_cuda('after forward')
    loss.backward()
    optimizer.step()
    stat_cuda('after backward')
    print("Epoch finished: %d -- Loss: %.6f" % (epoch, loss))
    toc()   

    '''Save model checkpoint'''
    torch.save({
                'epoch': epoch,
                'loss_iter': loss_iter,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict()
                }, savedir + 'model_e%d.pt' % (epoch))

'''Save model checkpoint'''
torch.save({
            'epoch': epoch,
            'loss_iter': loss_iter,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict()
            }, savedir + 'model_e%d.pt' % (epoch))


'''Plot spin-wave propagation'''
with torch.no_grad():

    ## Save Msat
    # spintorch.plot.geometry(model, epoch=epoch, plotdir=plotdir)
    Msat = model.geom.Msat.detach()
    savemat(savedir + "Msat.mat", {"Msat": Msat.to(torch.device("cpu")).numpy().transpose()})

    ## Save my
    my = model.m_last[0,1,]-model.m0[0,1,].unsqueeze(0).cpu()
    # wave_snapshot(model, my, (plotdir+'snapshot_time%d_epoch%d.png' % (timesteps,epoch)),r"$m_z$")
    # wave_integrated(model, my, (plotdir+'integrated_epoch%d.png' % (epoch)))
    savemat(savedir+"my.mat", {"my": my.to(torch.device("cpu")).numpy().transpose()})

    # Plot my
    plt.imshow(my[0,].transpose(0,1))
    plt.savefig(plotdir+"my.png")
    
    my = model.m_half[0,1,]-model.m0[0,1,].unsqueeze(0).cpu()
    # wave_snapshot(model, my, (plotdir+'snapshot_time%d_epoch%d.png' % (timesteps/2,epoch)),r"$m_z$")
    ## Save input
    savemat(savedir+"src.mat", {"src": INPUTS.to(torch.device("cpu")).numpy()})

    # Plot and save alpha
    alpha = model.Alpha(False).to(torch.device("cpu")).numpy().transpose()
    alpha2D = np.flip(alpha[:,:,0,0],0)
    plt.imshow(alpha2D)
    plt.savefig(plotdir+"alpha.png")
    savemat(savedir+"alpha.mat", {"alpha": alpha})

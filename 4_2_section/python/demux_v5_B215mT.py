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
basedir = 'demux_v5_B215mT/'
plotdir = 'plots/' + basedir
if not os.path.isdir(plotdir):
    os.makedirs(plotdir)
savedir = 'models/' + basedir
if not os.path.isdir(savedir):
    os.makedirs(savedir)    

"""Parameters"""
dx = 100e-9      # discretization (m)
dy = 100e-9      # discretization (m)
dz = 100e-9      # discretization (m)
nx = 700        # size x    (cells)
ny = 740        # size y    (cells)
damping_with = 60
desired_output = 4

Ms0 = 115.4e3      # saturation magnetization (A/m)
Ms_max = 117.39e3      # saturation magnetization (A/m)
Ms_min = 115.4e3      # saturation magnetization (A/m)
B0 = 215e-3     # bias field (T)
Bt = 1e-3       # excitation field amplitude (T)

dt = 20e-12     # timestep (s)
f1 = 2.1e9        # source frequency 1 (Hz)
f2 = 2.15e9    # source frequency 2 (Hz)
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
Msat = loadmat('Msat_walls.mat').get('Msat')
Msat = torch.tensor(Msat)
print(Msat)
model.geom.Msat0[:,:] = torch.transpose(Msat,0,1)
# print(model.geom.Msat0)

'''Define the source signal'''
## 2 input frequencies
t = torch.arange(0, timesteps*dt, dt, device=dev).unsqueeze(0).unsqueeze(2) # time vector
X1 = Bt*torch.sin(2*np.pi*f1*t)
X2 = Bt*torch.sin(2*np.pi*f2*t)

# INPUTS = torch.cat((X1, X2))  # here we could cat multiple inputs
# OUTPUTS = torch.tensor([8,Np-1-8]).to(dev) # desired output

'''Define optimizer and lossfunction'''
optimizer = torch.optim.Adam(model.parameters(), lr=1) # 0.01

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
model.retain_history = False
for epoch in range(epoch_init+1, 60):
    for i in range(2):
        '''Define input signal and it's output'''
        if i == 0:
            INPUTS = X1
            OUTPUTS = torch.tensor([int(desired_output)]).to(dev) # desired output
        else:
            INPUTS = X2
            OUTPUTS = torch.tensor([int(Np-1-desired_output)]).to(dev) # desired output
        

        '''Train the network'''
        tic()
        optimizer.zero_grad()
        u = model(INPUTS).sum(dim=1)

        # Plot
        spintorch.plot.plot_output(u[0,], OUTPUTS[0]+1, epoch, plotdir)
        spintorch.plot.geometry(model, epoch=-1, plotdir=plotdir)
        my = model.m_last[0,1,]-model.m0[0,1,].unsqueeze(0).cpu()
        wave_snapshot(model, my, (plotdir+'snapshot_last.png'),r"$m_z$")
        
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

    
'''Plot spin-wave propagation'''
model.retain_history = False
with torch.no_grad():
    
    ## Save Msat
    spintorch.plot.geometry(model, epoch=epoch, plotdir=plotdir)
    Msat = model.geom.Msat.detach()
    savemat(savedir + "Msat.mat", {"Msat": Msat.to(torch.device("cpu")).numpy().transpose()})
    
    for i in range(2):
        if i == 0:
            INPUTS = X1
        else:
            INPUTS = X2
        
        ## Calculate wave propagation
        tic()
        u = model(INPUTS).sum(dim=1)
        stat_cuda("after wave propagation")
        toc()

        ## Save my
        my = model.m_last[0,1,]-model.m0[0,1,].unsqueeze(0).cpu()
        wave_snapshot(model, my, (plotdir+'snapshot_time%d_epoch%d_input%d.png' % (timesteps,epoch,i)),r"$m_z$")
        wave_integrated(model, my, (plotdir+'integrated_epoch%d_input%d.png' % (epoch,i)))
        savemat(savedir+"my_input%d.mat" % (i), {"my": my.to(torch.device("cpu")).numpy().transpose()})

        my = model.m_half[0,1,]-model.m0[0,1,].unsqueeze(0).cpu()
        wave_snapshot(model, my, (plotdir+'snapshot_time%d_epoch%d_input%d.png' % (timesteps/2,epoch,i)),r"$m_z$")
        ## Save input
        savemat(savedir+"src.mat", {"src": INPUTS.to(torch.device("cpu")).numpy()})

        # Plot and save alpha
        import matplotlib.pyplot as plt
        alpha = model.Alpha(False).to(torch.device("cpu")).numpy().transpose()
        alpha2D = np.flip(alpha[:,:,0,0],0)
        plt.imshow(alpha2D)
        plt.savefig(plotdir+"alpha.png")
        savemat(savedir+"alpha.mat", {"alpha": alpha})
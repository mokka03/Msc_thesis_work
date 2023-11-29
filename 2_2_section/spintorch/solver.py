"""Micromagnetic solver with backpropagation"""

import torch 
from torch import nn, cat, cross, tensor, zeros, empty
from torch.utils.checkpoint import checkpoint
import numpy as np
from scipy.io import savemat ##

from .geom import WaveGeometryMs 
from .demag import Demag 
from .exch import Exchange
from .damping import Damping
from .sot import SOT

# torch.backends.cudnn.deterministic = True
# torch.backends.cudnn.benchmark = False


class MMSolver(nn.Module):
    gamma_LL = 1.7595e11    # gyromagnetic ratio (rad/Ts)
    relax_timesteps = 100
    retain_history = False

    def __init__(self, geometry, dt: float, timesteps, sources=[], probes=[]):
        super().__init__()

        self.register_buffer("dt", tensor(dt))               # timestep (s)

        self.geom = geometry
        self.sources = nn.ModuleList(sources)
        self.probes = nn.ModuleList(probes)
        self.demag_2D = Demag(self.geom.dim, self.geom.d)
        self.exch_2D = Exchange(self.geom.d)
        self.Alpha = Damping(self.geom.dim,self.geom.damping_with)
        self.torque_SOT = SOT(self.geom.dim)
        SOT.gamma_LL = self.gamma_LL
        # self.torque_SOT.active=False

        m0 = zeros((1, 3,) + self.geom.dim)
        # m0[:, 1,] =  1     # set initial magnetization in y direction
        m0[:, 2,] =  1     ## OOP magnetization
        self.register_buffer("m0", m0)
        self.m_history = []
        self.timesteps = timesteps
        self.tlasts = int(self.timesteps - self.timesteps/100)
        self.register_buffer("m_lasts", torch.zeros((self.timesteps-self.tlasts,3,self.geom.dim[0],self.geom.dim[1])))
        self.fwd = False  # differentiates main fwd pass and checpointing runs

    def forward(self, signal):
        self.fwd = True
        if isinstance(self.geom, WaveGeometryMs):
            Msat = self.geom()
            B_ext = self.geom.B
        else:
            B_ext = self.geom()
            Msat = self.geom.Ms

        self.relax(B_ext, Msat) # relax magnetization and store in m0 (no gradients)
        outputs = self.run(self.m0, B_ext, Msat, signal) # run the simulation
        self.fwd = False
        return cat(outputs, dim=1)

    def run(self, m, B_ext, Msat, signal):
        """Run the simulation in multiple stages for checkpointing"""
        outputs = []
        N = int(np.sqrt(signal.size()[1])) # number of stages 
        B_ext_0 = B_ext ##
        self.tstep = 1  ##
        for stage, sig in enumerate(signal.chunk(N, dim=1)):
            output, m = checkpoint(self.run_stage, m, B_ext, B_ext_0, Msat, sig)    ## B_ext_0
            outputs.append(output)
            # self.ii+=1   ##
        return outputs
        
    def run_stage(self, m, B_ext, B_ext_0, Msat, signal):   ## B_ext_0
        """Run a subset of timesteps (needed for 2nd level checkpointing)"""
        outputs = empty(0,device=self.dt.device)
        # Loop through the signal 
        # jj = 1  ##
        for sig in signal.split(1,dim=1):
            B_ext = self.inject_sources(B_ext, B_ext_0, sig)    ## B_ext_0
            # Propagate the fields (with checkpointing to save memory)
            m = checkpoint(self.rk4_step_LLG, m, B_ext, Msat)
            # Measure the outputs (checkpointing helps here as well)
            outputs = checkpoint(self.measure_probes, m, Msat, outputs)
            if self.fwd:    ##
                if self.tstep == int(self.timesteps/2):    ##
                    self.m_half = m.detach().cpu()    ##
                elif self.tstep == self.timesteps:    ##
                    self.m_last = m.detach().cpu()    ##
                if self.tstep > self.tlasts:    ##
                    self.m_lasts[self.tlasts-self.tstep,:,:,:] = m[0,]  ##
                    # self.m_lasts.append(m.detach().cpu())    ##
                if self.retain_history:    ##
                    self.m_history.append(m.detach().cpu())    ##
            self.tstep +=1  ##
        return outputs, m
        
    def inject_sources(self, B_ext, B_ext_0, sig):   ## B_ext_0
        """Add the excitation signal components to B_ext"""
        for i, src in enumerate(self.sources):
            B_ext = src(B_ext, B_ext_0, sig[0,0,i])  ## B_ext_0
        return B_ext

    def measure_probes(self, m, Msat, outputs): 
        """Extract outputs and concatenate to previous values"""
        probe_values = []
        for probe in self.probes:
            probe_values.append(probe((m-self.m0)*Msat))
        outputs = cat([outputs, cat(probe_values).unsqueeze(0).unsqueeze(0)], dim=1)
        return outputs
        
    def relax(self, B_ext, Msat): 
        """Run the solver with high damping to relax magnetization"""
        with torch.no_grad():
            for n in range(self.relax_timesteps):
                self.m0 = self.rk4_step_LLG(self.m0, B_ext, Msat, relax=True)
                
    def B_eff(self, m, B_ext, Msat): 
        """Sum the field components to return the effective field"""
        return B_ext + self.exch_2D(m,Msat) + self.demag_2D(m,Msat) 
    
    def rk4_step_LLG(self, m, B_ext, Msat, relax=False):
        """Implement a 4th-order Runge-Kutta solver"""
        h = self.gamma_LL * self.dt  # this time unit is closer to 1
        k1 = self.torque_LLG(m, B_ext, Msat, relax)
        k2 = self.torque_LLG(m + h*k1/2, B_ext, Msat, relax)
        k3 = self.torque_LLG(m + h*k2/2, B_ext, Msat, relax)
        k4 = self.torque_LLG(m + h*k3, B_ext, Msat, relax)
        return (m + h/6 * ((k1 + 2*k2) + (2*k3 + k4)))
    
    def torque_LLG(self, m, B_ext, Msat, relax=False):
        """Calculate Landau-Lifshitz-Gilbert torque"""
        m_x_Beff = cross(m, self.B_eff(m, B_ext, Msat),1)
        return (-(1/(1+self.Alpha(relax)**2) * (m_x_Beff + self.Alpha(relax)*cross(m, m_x_Beff,1))) 
                + self.torque_SOT(m, Msat))
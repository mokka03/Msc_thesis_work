import torch
import skimage
from torch.utils.checkpoint import checkpoint
from .utils import curl


class WaveProbe(torch.nn.Module):
	def __init__(self, x, y):
		super().__init__()

		self.register_buffer('x', torch.tensor(x, dtype=torch.int64))
		self.register_buffer('y', torch.tensor(y, dtype=torch.int64))

	def forward(self, m):
		return m[:,0, self.x, self.y]

	def coordinates(self):
		return self.x.cpu().numpy(), self.y.cpu().numpy()

class WaveIntensityProbe(WaveProbe):
	def __init__(self, x, y):
		super().__init__(x, y)

	def forward(self, m):
		return super().forward(m).pow(2)

class WaveIntensityProbeDisk(WaveProbe):
	def __init__(self, x, y, r):
		x, y = skimage.draw.disk((x, y), r)
		super().__init__(x, y)

	def forward(self, m):
		return super().forward(m).sum().pow(2).unsqueeze(0)
	
class PickupCPW(torch.nn.Module):
	def __init__(self, cpw_xpos, cpw_zpos, N, dr, dt, timesteps):
		super().__init__()
		self.register_buffer('cpw_xpos', torch.tensor(cpw_xpos))
		self.register_buffer('cpw_zpos', torch.tensor(cpw_zpos))
		self.register_buffer('dx', torch.tensor(dr[0]))
		self.register_buffer('dy', torch.tensor(dr[1]))
		self.register_buffer('dz', torch.tensor(dr[2]))
		self.register_buffer('dt', torch.tensor(dt))
		self.register_buffer('nx', torch.tensor(N[0]))
		self.register_buffer('ny', torch.tensor(N[1]))
		self.register_buffer('nt', torch.tensor(N[2]))
		self.mu0 = 4.0*torch.pi*1.0E-7

		self.register_buffer('time', torch.arange(timesteps-self.nt+1,timesteps) * self.dt) ## Ezt még átírni, a timesteps nem kell bele
		self.register_buffer('X', torch.zeros([self.nx,self.ny]))
		self.register_buffer('Y', torch.zeros([self.nx,self.ny]))
		x = torch.arange(0,self.nx)*self.dx
		y = torch.arange(0,self.ny)*self.dy
		self.X, self.Y = torch.meshgrid(x, y, indexing='ij')

		# self.register_buffer('distance_curl_M', torch.zeros([self.nx,self.ny]))

		# self.register_buffer('deltaV', torch.zeros((self.nt-1,self.cpw_xpos.shape[0])))

		# self.register_buffer('amplitude', torch.zeros(len(self.cpw_xpos)))
		# self.register_buffer('phase', torch.zeros(len(self.cpw_xpos)))

		self.register_buffer('CPW_signal', torch.tensor([0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0])/4)
		self.register_buffer('CPW_ground', torch.tensor([1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1])/4)
		self.register_buffer('f', torch.tensor(6e9))
		



	def forward(self, m, Msat):
		deltaV = torch.zeros((self.nt-1,self.cpw_xpos.shape[0])).to(torch.device('cuda'))
		amplitude =  torch.zeros(len(self.cpw_xpos)).to(torch.device('cuda'))
		phase = torch.zeros(len(self.cpw_xpos)).to(torch.device('cuda'))
		for tt in range(self.nt):
			print('Working on %d out of %d' %(tt+1,self.nt))
			self.m_ = m[tt,]*Msat
			self.m_ = self.m_[:,:,:,None]
			
			### Computing the curl of magnetization (M)
			curl_M_datax, curl_M_datay, curl_M_dataz = curl(torch.cat((self.m_[0,],self.m_[0,]), 2), # M_datax
												   torch.cat((self.m_[1,],self.m_[1,]), 2), # M_datay
												   torch.cat((self.m_[2,],self.m_[2,]), 2), # M_dataz
												   self.dx,self.dy,self.dz)

			### Get vectorpotential
			# A_x, A_y, A_z = self.get_vector_porential(curl_M_datax, curl_M_datay, curl_M_dataz)
			A_y = checkpoint(self.get_vector_porential, curl_M_datax, curl_M_datay, curl_M_dataz)

			if tt > 0:
				###### Calculating electric field
				# E_x = - (A_x - A_x_prev)/self.dt
				E_y = - (A_y - A_y_prev)/self.dt
				# E_z = - (A_z - A_z_prev)/self.dt

				for filament  in range(self.cpw_xpos.shape[0]):
					deltaV[tt-1,filament] = -torch.sum(E_y[filament,:],0)*self.dy

			# A_x_prev = A_x
			A_y_prev = A_y
			# A_z_prev = A_z


		ffit_save = torch.zeros_like(deltaV)
		for filament  in range(self.cpw_xpos.shape[0]):
			# fit Fourier curve;
			sin_vec = torch.sin(2*torch.pi*self.f*self.time)
			cos_vec = torch.cos(2*torch.pi*self.f*self.time)
			re = torch.sum(cos_vec*deltaV[:,filament])/torch.pi/(cos_vec.shape[0]*self.dt*self.f)*self.dt*self.f*2*torch.pi
			im = torch.sum(sin_vec*deltaV[:,filament])/torch.pi/(cos_vec.shape[0]*self.dt*self.f)*self.dt*self.f*2*torch.pi
			ffit = re*cos_vec+im*sin_vec
			ffit_save[:,filament] = ffit
			amplitude[filament] = torch.abs(re+1j*im)
			phase[filament] = torch.angle(re+1j*im)
		avg = torch.sum(amplitude*torch.exp(1j*phase)*self.CPW_signal) - torch.sum(amplitude*torch.exp(1j*phase)*self.CPW_ground)

		# return torch.abs(avg)
		return deltaV, ffit_save
	

	def get_vector_porential(self, curl_M_datax, curl_M_datay, curl_M_dataz):
		###### Calculating vector potential
		A_datax = torch.zeros_like(self.X)
		A_datay = torch.zeros_like(self.X)
		A_dataz = torch.zeros_like(self.X)

		for j in range(self.ny):
			for i in self.cpw_xpos:

				# distance_curl_M = torch.sqrt(torch.square(self.X-self.X[i,j]) + torch.square(self.Y-self.Y[i,j]) + torch.square(self.cpw_zpos))

				#### Volume
				# fractionx = curl_M_datax / torch.sqrt(torch.square(self.X-self.X[i,j]) + torch.square(self.Y-self.Y[i,j]) + torch.square(self.cpw_zpos))
				fractiony = curl_M_datay / torch.sqrt(torch.square(self.X-self.X[i,j]) + torch.square(self.Y-self.Y[i,j]) + torch.square(self.cpw_zpos))
				# fractionz = curl_M_dataz / torch.sqrt(torch.square(self.X-self.X[i,j]) + torch.square(self.Y-self.Y[i,j]) + torch.square(self.cpw_zpos))

				# A_datax[i,j] = self.mu0/(4*torch.pi)*self.dx*self.dy*self.dz*torch.sum(torch.nan_to_num(fractionx))
				A_datay[i,j] = self.mu0/(4*torch.pi)*self.dx*self.dy*self.dz*torch.sum(torch.nan_to_num(fractiony))
				# A_dataz[i,j] = self.mu0/(4*torch.pi)*self.dx*self.dy*self.dz*torch.sum(torch.nan_to_num(fractionz))

		# return A_datax[self.cpw_xpos,:], A_datay[self.cpw_xpos,:], A_dataz[self.cpw_xpos,:]
		return A_datay[self.cpw_xpos,:]

import torch
import skimage
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

		self.register_buffer('A_x', torch.zeros((len(self.cpw_xpos), self.ny, self.nt)))
		self.register_buffer('A_y', torch.zeros((len(self.cpw_xpos), self.ny, self.nt)))
		self.register_buffer('A_z', torch.zeros((len(self.cpw_xpos), self.ny, self.nt)))
		# self.register_buffer('E_x', torch.zeros((len(self.cpw_xpos), self.ny, self.nt)))
		# self.register_buffer('E_y', torch.zeros((len(self.cpw_xpos), self.ny, self.nt)))
		# self.register_buffer('E_z', torch.zeros((len(self.cpw_xpos), self.ny, self.nt)))

		self.register_buffer('deltaV', torch.zeros(self.cpw_xpos.shape[0]+1))

		# self.register_buffer('m_t', torch.zeros((3,self.nx,self.ny)))
		self.register_buffer('m_', torch.zeros((3,self.nx,self.ny,2)))

		self.register_buffer('n_hat_z', torch.zeros((3,self.nx,self.ny)))
		self.register_buffer('n_hat_x', torch.zeros_like(self.n_hat_z[:,1,:]))
		self.register_buffer('n_hat_y', torch.zeros_like(self.n_hat_z[:,:,1]))
		self.register_buffer('C_00x', torch.zeros_like(self.n_hat_z))
		self.register_buffer('C_100', torch.zeros_like(self.n_hat_z[:,1,:]))
		self.register_buffer('C_n00', torch.zeros_like(self.n_hat_z[:,1,:]))
		self.register_buffer('C_010', torch.zeros_like(self.n_hat_z[:,:,1]))
		self.register_buffer('C_0n0', torch.zeros_like(self.n_hat_z[:,:,1]))

		self.register_buffer('amplitude', torch.zeros(self.A_x.shape[0]))
		self.register_buffer('phase', torch.zeros(self.A_x.shape[0]))

		self.register_buffer('CPW_signal', torch.tensor([0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0])/4)
		self.register_buffer('CPW_ground', torch.tensor([1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1])/4)
		self.register_buffer('f', torch.tensor(6e9))
		



	def forward(self, m, Msat):
		for tt in range(self.nt):
			print('Working on %d out of %d' %(tt+1,self.nt))
			self.m_t = m[tt,]
			self.m_[:,:,:,0] = self.m_t*Msat
			self.m_[:,:,:,1] = self.m_t*Msat
			M_datax = self.m_[0,]
			M_datay = self.m_[1,]
			M_dataz = self.m_[2,]
			### Computing the curl of magnetization (M)
			curl_M_datax, curl_M_datay, curl_M_dataz = curl(M_datax,M_datay,M_dataz,self.dx,self.dy,self.dz)
			curl_M_datax = curl_M_datax[:,:,0]
			curl_M_datay = curl_M_datay[:,:,0]
			curl_M_dataz = curl_M_dataz[:,:,0]

			###### Calculating vector potential
			A_datax = torch.zeros_like(self.X)
			A_datay = torch.zeros_like(self.X)
			A_dataz = torch.zeros_like(self.X)

			#### Cross product for surface intergal
			### Top
			self.n_hat_z[2,] = 1
			self.C_00x = torch.cross(self.m_t,self.n_hat_z,dim=0)
			### Bottom
			self.n_hat_z[2,] = -1
			self.C_00x = self.C_00x + torch.cross(self.m_t,self.n_hat_z,dim=0)
			### Right
			self.n_hat_x[0,] = 1
			self.C_100 = torch.cross(self.m_t[:,self.nx-1,:],self.n_hat_x,dim=0)
			### Left
			self.n_hat_x[0,] = -1
			self.C_n00 = torch.cross(self.m_t[:,0,:],self.n_hat_x,dim=0)
			### Up
			self.n_hat_y[1,] = 1
			self.C_010 = torch.cross(self.m_t[:,:,self.ny-1],self.n_hat_y,dim=0)
			### Down
			self.n_hat_y[1,] = -1
			self.C_0n0 = torch.cross(self.m_t[:,:,0],self.n_hat_y,dim=0)

			for j in range(self.ny):
				for i in self.cpw_xpos:

					distance_curl_M = torch.sqrt(torch.square(self.X-self.X[i,j]) + torch.square(self.Y-self.Y[i,j]) + torch.square(self.cpw_zpos))

					#### Volume
					fractionx = curl_M_datax / distance_curl_M
					fractiony = curl_M_datay / distance_curl_M
					fractionz = curl_M_dataz / distance_curl_M

					A_datax[i,j] = self.mu0/(4*torch.pi)*self.dx*self.dy*self.dz*torch.sum(torch.nan_to_num(fractionx))
					A_datay[i,j] = self.mu0/(4*torch.pi)*self.dx*self.dy*self.dz*torch.sum(torch.nan_to_num(fractiony))
					A_dataz[i,j] = self.mu0/(4*torch.pi)*self.dx*self.dy*self.dz*torch.sum(torch.nan_to_num(fractionz))

					#### Surface
					### Top and bottom
					surf_frac = self.C_00x / distance_curl_M
					A_s_00x = self.mu0/(4.0*torch.pi)*self.dx*self.dy*torch.sum(torch.nan_to_num(surf_frac),(1,2))
					### Right
					surf_frac = self.C_100 / distance_curl_M[self.nx-1,:]
					A_s_100 = self.mu0/(4.0*torch.pi)*self.dy*self.dz*torch.sum(torch.nan_to_num(surf_frac),1)
					### Left
					surf_frac = self.C_n00 / distance_curl_M[1,:]
					A_s_n00 = self.mu0/(4.0*torch.pi)*self.dy*self.dz*torch.sum(torch.nan_to_num(surf_frac),1)
					### Up
					surf_frac = self.C_010 / distance_curl_M[:,self.ny-1]
					A_s_010 = self.mu0/(4.0*torch.pi)*self.dx*self.dz*torch.sum(torch.nan_to_num(surf_frac),1)
					### Down
					surf_frac = self.C_0n0 / distance_curl_M[:,0]
					A_s_0n0 = self.mu0/(4.0*torch.pi)*self.dx*self.dz*torch.sum(torch.nan_to_num(surf_frac),1)

					A_datax[i,j] = A_datax[i,j] + A_s_00x[0] + A_s_100[0] + A_s_n00[0] + A_s_010[0] + A_s_0n0[0]
					A_datay[i,j] = A_datay[i,j] + A_s_00x[1] + A_s_100[1] + A_s_n00[1] + A_s_010[1] + A_s_0n0[1]
					A_dataz[i,j] = A_dataz[i,j] + A_s_00x[2] + A_s_100[2] + A_s_n00[2] + A_s_010[2] + A_s_0n0[2]

			self.A_x[:,:,tt] = A_datax[self.cpw_xpos,:]
			self.A_y[:,:,tt] = A_datay[self.cpw_xpos,:]
			self.A_z[:,:,tt] = A_dataz[self.cpw_xpos,:]

			###### Calculating electric field
			self.E_x = - (self.A_x[:,:,1:self.nt] - self.A_x[:,:,0:self.nt-1])/self.dt
			self.E_y = - (self.A_y[:,:,1:self.nt] - self.A_y[:,:,0:self.nt-1])/self.dt
			self.E_z = - (self.A_z[:,:,1:self.nt] - self.A_z[:,:,0:self.nt-1])/self.dt


		
		for filament  in range(self.E_x.shape[0]):
			self.deltaV = -torch.sum(self.E_y[filament,:,:],0)*self.dy

			# fit Fourier curve;
			c = (self.dt*self.f)/0.16 # correction term
			sin_vec = torch.sin(2*torch.pi*self.f*self.time)
			cos_vec = torch.cos(2*torch.pi*self.f*self.time)
			re = torch.sum(cos_vec*self.deltaV)/torch.pi/(cos_vec.shape[0]*self.dt*self.f)*c
			im = torch.sum(sin_vec*self.deltaV)/torch.pi/(cos_vec.shape[0]*self.dt*self.f)*c
			ffit = re*cos_vec+im*sin_vec

			self.amplitude[filament] = torch.abs(re+1j*im)
			self.phase[filament] = torch.angle(re+1j*im)
		avg = torch.sum(self.amplitude*torch.exp(1j*self.phase)*self.CPW_signal) - torch.sum(self.amplitude*torch.exp(1j*self.phase)*self.CPW_ground)

		return torch.abs(avg)